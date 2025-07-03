import click

import fishsim.src.run_merfish
import fishsim.src.subdivide

@click.command(
    context_settings={"max_content_width": 120},
    name="run_merfish"
)
@click.option('--config-file', required=True, type=click.Path(exists=True), help='Path to the configuration file (config.yml).')
@click.option('--output-dir-name', required=True, type=str, help='Define the name of the synethic data directory')
def run_merfish(**kwargs):
    """ Generate synthetic data according to merFISH barcoding scheme.
    """
    fishsim.src.run_merfish(**kwargs)


@click.group(name="fishsim")
def main():
    pass


main.add_command(run_merfish)
