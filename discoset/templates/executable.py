EXECUTABLE_TEMPLATE = '''#!/usr/bin/env python3
"""
Date generated: {date}

This one-file script takes an input of molecules in SMILES format
and sorts them by the number of matches to rare environments.

The rare environments include a list of functional groups and {forcefield_name} parameters
for which there is low coverage in our already available data.

"""

import argparse
import contextlib
import functools
import multiprocessing
import tempfile

from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
import pandas as pd

parser = argparse.ArgumentParser(
    description=(
        "Sort a list of SMILES strings by the number of matches to rare environments. "
        "This takes in a list of SMILES strings and outputs a CSV file with the SMILES "
        "strings and the number of matches to rare environments.\\n"
        "The rare environments include a list of functional groups and {forcefield_name} parameters "
        "for which there is low coverage in our already available data.\\n"
        "The output CSV file will contain columns for SMILES strings (column: 'SMILES', type: str) "
        "and the number of matches to rare environments (column: 'Count', type: int). "
        "If the flag --no-write-all is set, only these two columns will be written. "
        "If the flag is not set, all the rare environments that were searched for will be written "
        "as well, with boolean values to indicate if they are present or not."
    ),
    formatter_class=argparse.RawDescriptionHelpFormatter  # preserve newlines
)
parser.add_argument(
    "-i", "--input",
    type=str,
    help="Path to a file containing SMILES strings, with one on each line.",
)
parser.add_argument(
    "-o", "--output",
    type=str,
    help="Path to the output CSV file.",
)
parser.add_argument(
    "-np", "--nproc",
    type=int,
    default=1,
    help="Number of processes to use.",
)
parser.add_argument(
    "--no-write-all",
    action="store_false",
    help=(
        "Write all columns to the CSV file. "
        "If False, only columns for SMILES and the number "
        "matches to rare groups will be written. "
        "If True, all the rare groups that were searched for will be written."
    )
)
parser.add_argument(
    "-c", "--count-threshold",
    type=int,
    default=1,
    help=(
        "Threshold for considering a parameter 'rare'. "
        "Only groups with a count greater than or equal to this threshold will be written."
    )
)
parser.add_argument(
    "-n", "--only-top-n",
    type=int,
    default=-1,
    help="Only write the top N entries to the output file.",
)



def main():
    args = parser.parse_args()
    with open(args.input) as f:
        smiles = [line.strip() for line in f]
    search_all_smiles(
        smiles,
        args.output,
        nprocs=args.nproc,
        write_all=args.no_write_all,
        count_threshold=args.count_threshold,
        only_top_n=args.only_top_n,
    )


def draw_checkmol():
    """
    Draw the checkmol SMARTS for verification.
    
    This is not used in the main script, but can be used to generate a
    visual representation of the checkmol groups.
    """
    from rdkit import Chem
    from rdkit.Chem import Draw

    rdmols = [
        Chem.MolFromSmarts(smirks)
        for smirks in CHECKMOL_GROUPS.values()
    ]
    img = Draw.MolsToGridImage(
        rdmols,
        molsPerRow=4,
        subImgSize=(300, 300),
        legends=list(CHECKMOL_GROUPS.keys()),
    )
    img.save("checkmol.png")


def search_all_smiles(
    smiles,
    output_file: str,
    nprocs: int = 1,
    write_all: bool = True,
    count_threshold: int = 1,
    only_top_n: int = -1,
):
    """
    Search all SMILES strings for rare environments and write to a CSV file.

    Parameters
    ----------
    smiles : list[str]
        List of SMILES strings to search.
    output_file : str
        Path to the output CSV file.
    nprocs : int, optional
        Number of processors to use, by default 1.
    count_threshold : int, optional
        Threshold for considering a parameter 'rare', by default 1.
    only_top_n : int, optional
        Only write the top N entries to the output file. Default -1.
    write_all : bool, optional
        Write all columns to the CSV file. If False, only columns for SMILES
        and the number of diverse environments will be written.
        If True, all the rare groups that were searched for will be written
        as well. Default True.
    """
    forcefield = load_forcefield()
    empty_entry = {{ group: False for group in CHECKMOL_GROUPS }}
    for parameter_id in LOW_COVERAGE_PARAMETERS:
        empty_entry[parameter_id] = False
    all_entries = []

    labeller = functools.partial(
        label_single_smiles,
        forcefield=forcefield,
        empty_entry=empty_entry,
    )
    with multiprocessing.Pool(nprocs) as pool:
        all_entries = list(
            _progress_bar(
                pool.imap(labeller, smiles),
                desc="Labeling SMILES",
                total=len(smiles),
            )
        )
    all_entries = [x for x in all_entries if x is not None]

    all_entries = sorted(all_entries, key=lambda x: x["Count"], reverse=True)
    df = pd.DataFrame(all_entries)

    if not len(df):
        print(f"No valid matches found -- skipping writing to {{output_file}}")

    keys = ["SMILES", "Count"]
    if write_all:
        keys += list(empty_entry.keys())
    
    df = df[keys]
    df = df[df["Count"] >= count_threshold]
    if only_top_n > 0:
        df = df.head(only_top_n)

    df.to_csv(output_file, index=False)
    print(f"Wrote {{len(df)}} entries to {{output_file}}")



def _progress_bar(iterable_, **kwargs):
    """Try to use tqdm if it is available, otherwise return the iterable."""
    try:
        import tqdm
        return tqdm.tqdm(iterable_, **kwargs)
    except ImportError:
        return iterable_

def load_forcefield():
    """Load the OpenFF 2.2.0 force field from the string below."""
    # write force field to temp file and load
    with tempfile.NamedTemporaryFile("w") as f:
        f.write(FORCEFIELD)
        forcefield = ForceField(f.name)
    return forcefield

def label_single_smiles(
    smi: str,
    forcefield: ForceField,
    empty_entry,
):
    """Label a single SMILES string with rare environments."""
    # ignore warnings about stereo
    with capture_toolkit_warnings():
        try:
            mol = Molecule.from_smiles(smi, allow_undefined_stereo=True)
        except Exception as e:
            return None

    entry = dict(empty_entry)

    # does it match any checkmol groups?
    for group, smirks in CHECKMOL_GROUPS.items():
        matches = mol.chemical_environment_matches(smirks)
        if len(matches):
            entry[group] = True

    labels = forcefield.label_molecules(mol.to_topology())[0]
    for parameters in labels.values():
        for parameter in parameters.values():
            label = parameter.id if parameter.id else parameter.name
            if label in entry:
                entry[label] = True

    entry["Count"] = sum(entry.values())
    entry["SMILES"] = smi
    return entry

@contextlib.contextmanager
def capture_toolkit_warnings(run: bool = True):  # pragma: no cover
    """A convenience method to capture and discard any warning produced by external
    cheminformatics toolkits excluding the OpenFF toolkit. This should be used with
    extreme caution and is only really intended for use when processing tens of
    thousands of molecules at once."""

    import logging
    import warnings

    if not run:
        yield
        return

    warnings.filterwarnings("ignore")

    toolkit_logger = logging.getLogger("openff.toolkit")
    openff_logger_level = toolkit_logger.getEffectiveLevel()
    toolkit_logger.setLevel(logging.ERROR)

    yield

    toolkit_logger.setLevel(openff_logger_level)






#  ____                                                           
# |  _ \\__ _ _ __ ___                                            
# | |_) / _` | '__/ _ \\                                           
# |  _ < (_| | | |  __/                                           
# |_| \\_\\__,_|_|  \\___|                                   _       
#   ___ _ ____   _(_)_ __ ___  _ __  _ __ ___   ___ _ __ | |_ ___ 
#  / _ \\ '_ \\ \\ / / | '__/ _ \\| '_ \\| '_ ` _ \\ / _ \\ '_ \\| __/ __|
# |  __/ | | \\ V /| | | | (_) | | | | | | | | |  __/ | | | |_\\__ \\
#  \\___|_| |_|\\_/ |_|_|  \\___/|_| |_|_| |_| |_|\\___|_| |_|\\__|___/



# If updating this script, this section should be updated as well.

# checkmol groups:
# these were identified using the checkmol software.
CHECKMOL_GROUPS = {{
{checkmol_groups_string}
}}

LOW_COVERAGE_PARAMETERS = [
{low_coverage_parameters_string}
]

FORCEFIELD = """\\
{forcefield_xml_string}
"""


if __name__ == "__main__":
    main()
'''