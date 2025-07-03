import click

import anglerfish_simulator.src.run_merfish
import anglerfish_simulator.src.subdivide

@click.command(
    context_settings={"max_content_width": 120},
    name="run_merfish"
)
@click.option('--config-file', required=True, type=click.Path(exists=True), help='Path to the configuration file (config.yml).')
@click.option('--output-dir-name', required=True, type=str, help='Define the name of the synethic data directory')
def run_merfish(**kwargs):
    """ Generate synthetic data according to merFISH barcoding scheme.
    """
    anglerfish_simulator.src.run_merfish(**kwargs)



@click.command(
    context_settings={"max_content_width": 120},
    name="subdivide"
)
@click.option('--dataset-name', required=True, type=str, help='Define the name of the synethic data directory for simulated images, which will be used to create the chopped directory')
@click.option('--total-fov', required=True, type=int, help='total number of simulated tiles/fovs')
@click.option('--img-size', default=40, type=int, help='image size')
@click.option('--num-rounds', default=8, type=int, help='total number of imaging rounds per channel')
@click.option('--wv-beads', default=473, type=int, help='wavelength channel for beads')
@click.option('--wv-nuclei', default=100, type=int, help='wavelength channel for nuclei')
@click.option('--wv-spots1', default=650, type=int, help='wavelength channel 1 for spots')
@click.option('--wv-spots2', default=561, type=int, help='wavelength channel 2 for spots')
def subdivide(**kwargs):
    """ Subdivide synthetic data into 40x40 images for AnglerFISH model training
    """
    anglerfish_simulator.src.subdivide(**kwargs)



@click.group(name="anglerfish_simulator")
def main():
    pass


main.add_command(run_merfish)
main.add_command(subdivide)
