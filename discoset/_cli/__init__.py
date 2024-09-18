import click

from discoset.generate_new_patterns import generate

@click.command("generate")
@click.option(
    "-ff", "--forcefield-path",
    required=True,
    type=str,
    help="Path to the forcefield file",
)
@click.option(
    "-o", "--output-directory",
    required=True,
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    help="Path to the output directory to write the script to",
)
@click.option(
    "-v", "--verbose",
    is_flag=True,
    help="Print out more information",
)
@click.option(
    "-np", "--n-processes",
    default=1,
    type=int,
    help="Number of processes to use",
)
@click.option(
    "-fft", "--forcefield-threshold",
    default=10,
    type=int,
    help="Threshold for considering a parameter 'rare'",
)
@click.option(
    "-fgt", "--checkmol-threshold",
    default=10,
    type=int,
    help="Threshold for considering a parameter 'rare'",
)
@click.option(
    "-vn", "--version-number",
    default=1,
    type=int,
    help="Version number for the script",
)
@click.option(
    "-ocf", "--output-checkmol-file",
    default=None,
    type=click.Path(exists=False, file_okay=True, dir_okay=False),
    help=(
        "Output JSON file for checkmol pattern counts. "
        "If not provided, will not be written."
    ),
)
@click.option(
    "-off", "--output-forcefield-file",
    default=None,
    type=click.Path(exists=False, file_okay=True, dir_okay=False),
    help=(
        "Output JSON file for forcefield pattern counts. "
        "If not provided, will not be written."
    ),
)
def generate_cli(
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
    generate(
        forcefield_path=forcefield_path,
        output_directory=output_directory,
        verbose=verbose,
        n_processes=n_processes,
        forcefield_threshold=forcefield_threshold,
        checkmol_threshold=checkmol_threshold,
        version_number=version_number,
        output_checkmol_file=output_checkmol_file,
        output_forcefield_file=output_forcefield_file,
    )

@click.group()
def cli():
    pass

cli.add_command(generate_cli)



__all__ = ["cli"]