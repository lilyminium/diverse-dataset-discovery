import datetime
import functools
import json
import pathlib
import multiprocessing
import textwrap

import tqdm

from openff.toolkit import Molecule, ForceField
from discoset.data import SMARTS
from discoset.templates.executable import EXECUTABLE_TEMPLATE
from yammbs.checkmol import ChemicalEnvironment, analyze_functional_groups

def replace_v4_with_v3(forcefield_str: str):
    REPLACEMENTS = {
        r'<vdW version="0.4" potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1.0" cutoff="9.0 * angstrom ** 1" switch_width="1.0 * angstrom ** 1" periodic_method="cutoff" nonperiodic_method="no-cutoff">': r'<vdW version="0.3" potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1.0" cutoff="9.0 * angstrom" switch_width="1.0 * angstrom" method="cutoff">',
        r'<Electrostatics version="0.4" scale12="0.0" scale13="0.0" scale14="0.8333333333" scale15="1.0" cutoff="9.0 * angstrom ** 1" switch_width="0.0 * angstrom ** 1" periodic_potential="Ewald3D-ConductingBoundary" nonperiodic_potential="Coulomb" exception_potential="Coulomb">': r'<Electrostatics version="0.3" scale12="0.0" scale13="0.0" scale14="0.8333333333" scale15="1.0" cutoff="9.0 * angstrom" switch_width="0.0 * angstrom" method="PME">',
    }
    for old, new in REPLACEMENTS.items():
        forcefield_str = forcefield_str.replace(old, new)
    return forcefield_str


def single_label_mol_forcefield(mol: Molecule, forcefield, empty_entry: dict[str, bool]):
    entry = dict(empty_entry)
    labels = forcefield.label_molecules(mol.to_topology())[0]
    for parameters in labels.values():
        for parameter in parameters.values():
            label = parameter.id if parameter.id else parameter.name
            if label in entry:
                entry[label] = True
    return entry


def single_label_smiles_checkmol(smiles: str, empty_entry: dict[str, bool]):
    entry = dict(empty_entry)
    groups = analyze_functional_groups(smiles)
    if groups:
        for group in groups:
            entry[group.value] = True
    return entry


def single_label_inchi(
    inchi: str,
    forcefield: ForceField = None,
    empty_checkmol_entry: dict[str, bool] = {},
    empty_forcefield_entry: dict[str, bool] = {},
):
    try:
        mol = Molecule.from_inchi(inchi, allow_undefined_stereo=True)
        smiles = mol.to_smiles()
    except Exception as e:
        return
    checkmol_entry = single_label_smiles_checkmol(smiles, empty_checkmol_entry)
    forcefield_entry = single_label_mol_forcefield(mol, forcefield, empty_forcefield_entry)
    return {
        "inchi": inchi,
        "smiles": smiles,
        "checkmol": checkmol_entry,
        "forcefield": forcefield_entry,
    }


def get_all_inchis(
    verbose: bool = False,
) -> set[str]:
    """
    Download all available inchi keys from the QCArchive datasets.
    """
    from qcportal import PortalClient

    client = PortalClient("https://api.qcarchive.molssi.org:443", cache_dir=".")
    datasets = client.list_datasets()

    all_inchis = set()

    if verbose:
        datasets = tqdm.tqdm(datasets, desc="Downloading inchi keys from datasets")

    for dataset_data in datasets:
        dataset_id = dataset_data["id"]
        dataset = client.get_dataset_by_id(dataset_id)
        n_entries_no_standard_inchi = 0
        for entry in dataset.iterate_entries():
            if "standard_inchi" in entry.attributes:
                all_inchis.add(entry.attributes["standard_inchi"])
            else:
                n_entries_no_standard_inchi += 1
        print(
            f"{n_entries_no_standard_inchi} in dataset {dataset_id} do not have standard_inchi"
        )
    
    if verbose:
        print(f"Found {len(all_inchis)} unique InChI keys")

    return all_inchis




def generate(
    forcefield_path: str,
    output_directory: str,
    verbose: bool = False,
    n_processes: int = 1,
    forcefield_threshold: int = 10,
    checkmol_threshold: int = 10,
    version_number: int = 1,
    output_checkmol_file: str = None,
    output_forcefield_file: str = None,
):
    # get all inchi keys
    all_inchis = sorted(get_all_inchis(verbose=verbose))

    # load forcefield
    forcefield = ForceField(forcefield_path)

    # set up empty FF entry
    empty_forcefield_entry = {}
    handlers_to_remove = ["Constraints", "LibraryCharges", "ToolkitAM1BCC", "Electrostatics", "vdW"]
    for parameter_handler in forcefield._parameter_handlers:
        if parameter_handler in handlers_to_remove:
            continue
        handler = forcefield.get_parameter_handler(parameter_handler)
        for parameter in handler.parameters:
            label = parameter.id if parameter.id else parameter.name
            empty_forcefield_entry[label] = False

    # set up empty checkmol entry
    checkmol_columns = sorted([
        col.value for col in ChemicalEnvironment
    ])
    empty_checkmol_entry = dict.fromkeys(checkmol_columns, False)

    labeller = functools.partial(
        single_label_inchi,
        forcefield=forcefield,
        empty_checkmol_entry=empty_checkmol_entry,
        empty_forcefield_entry=empty_forcefield_entry,
    )

    with multiprocessing.Pool(n_processes) as pool:
        if verbose:
            results = list(
                tqdm.tqdm(
                    pool.imap(labeller, all_inchis),
                    total=len(all_inchis),
                    desc="Labelling molecules",
                )
            )
        else:
            results = list(pool.imap(labeller, all_inchis))
    
    results = [x for x in results if x is not None]

    # collate results and sum
    checkmol_results = dict.fromkeys(checkmol_columns, 0)
    forcefield_results = dict.fromkeys(empty_forcefield_entry.keys(), 0)
    for result in results:
        for key, value in result["checkmol"].items():
            if value:
                checkmol_results[key] += 1
        for key, value in result["forcefield"].items():
            if value:
                forcefield_results[key] += 1

    if output_checkmol_file:
        with open(output_checkmol_file, "w") as f:
            json.dump(checkmol_results, f)
    
    if output_forcefield_file:
        with open(output_forcefield_file, "w") as f:
            json.dump(forcefield_results, f)

    lowest_checkmol_groups = sorted(
        checkmol_results.items(), key=lambda x: x[1],
    )
    
    lowest_forcefield_groups = sorted(
        forcefield_results.items(), key=lambda x: x[1],
    )

    needed_checkmol_groups = []
    count = -10
    while lowest_checkmol_groups and count < checkmol_threshold:
        group, count = lowest_checkmol_groups.pop(0)
        if count < checkmol_threshold:
            needed_checkmol_groups.append(group)
    
    needed_forcefield_groups = {"b": [], "a": [], "t": [], "i": [], "n": []}
    count = -10
    while lowest_forcefield_groups and count < forcefield_threshold:
        group, count = lowest_forcefield_groups.pop(0)
        if count < forcefield_threshold:
            prefix = group[0]
            needed_forcefield_groups[prefix].append(group)


    # format strings
    checkmol_groups_lines = [
        f'"{group}": "{SMARTS.get(group, "[*:1]")}",  # {checkmol_results[group]} matches'
        for group in needed_checkmol_groups
    ]
    checkmol_groups_string = textwrap.indent(
        "\n".join(checkmol_groups_lines),
        " " * 4,
    )

    low_coverage_parameters_lines = []
    for prefix, groups in needed_forcefield_groups.items():
        for group in groups:
            low_coverage_parameters_lines.append(
                f'"{group}",  # {forcefield_results[group]} matches'
            )
        low_coverage_parameters_lines.append("")
    low_coverage_parameters_string = textwrap.indent(
        "\n".join(low_coverage_parameters_lines),
        " " * 4,
    )

    date = datetime.datetime.now().strftime("%Y-%m-%d")
    forcefield_name = pathlib.Path(forcefield_path).name
    forcefield_xml_string = replace_v4_with_v3(forcefield.to_string())

    executable_template = EXECUTABLE_TEMPLATE.format(
        date=date,
        forcefield_name=forcefield_name,
        checkmol_groups_string=checkmol_groups_string,
        low_coverage_parameters_string=low_coverage_parameters_string,
        forcefield_xml_string=forcefield_xml_string,
    )

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(exist_ok=True, parents=True)
    suffix = forcefield_path.split("-", 1)[1].split(".offxml")[0]

    output_file = output_directory / f"select-interesting-molecules_{suffix}_v{version_number}.py"
    with open(output_file, "w") as file:
        file.write(executable_template)
    print(f"Script written to {output_file}")
