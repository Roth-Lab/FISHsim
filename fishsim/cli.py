import click

import fishsim.run


@click.command(context_settings={"max_content_width": 120}, name="simulate")
@click.option(
    "-b",
    "--codebook-file",
    required=True,
    type=click.Path(exists=True),
    help="Path to the codebook file in TSV format.",
)
@click.option(
    "-c",
    "--config-file",
    required=True,
    type=click.Path(exists=True),
    help="Path to the simulation configuration file in YAML format.",
)
@click.option(
    "-d",
    "--data-org-file",
    required=True,
    type=click.Path(exists=True),
    help="Path to the data organisation file in TSV format.",
)
@click.option(
    "-e",
    "--emitter-file",
    required=True,
    type=click.Path(),
    help="Path to file to save the emitter file information to.",
)
@click.option(
    "-i",
    "--img-file",
    required=True,
    type=click.Path(),
    help="Path to file to save the stacked images in.",
)
@click.option(
    "--dist-file",
    type=click.Path(),
    help="Path to file to specifying distribution of gene expression.",
)
@click.option(
    "--seed",
    default=None,
    type=int,
    help="Random seed",
)
def simulate(**kwargs):
    """Generate synthetic data according to merFISH barcoding scheme."""
    fishsim.run.simulate(**kwargs)


@click.group(name="fishsim")
def main():
    pass


main.add_command(simulate)
