import click

import fishsim.src.run_merfish
import fishsim.src.subdivide


@click.command(context_settings={"max_content_width": 120}, name="simulate")
@click.option(
    "-c",
    "--config-file",
    required=True,
    type=click.Path(exists=True),
    help="Path to the configuration file in YAML format.",
)
@click.option(
    "-o",
    "--output-dir-name",
    required=True,
    type=click.Path(),
    help="Define the name of the synethic data directory",
)
def simulate(**kwargs):
    """Generate synthetic data according to merFISH barcoding scheme."""
    fishsim.src.run_merfish(**kwargs)


@click.group(name="fishsim")
def main():
    pass


main.add_command(run_merfish)
