# FISHsim (Multiplexed Spots Simulator)

This is a simulation software that generates photorealistic images of flourescent probes.

## Installation
We recommend using conda to set up a virtual environment. 

```
conda env create -f environments.yml
conda activate fishsim
pip install -e .
```


## Program Usage

## Running as an installed package
If you have installed fishsim as a package, then

1. Run in command line to simulate images.
```
fishsim simulate --config-file examples/config_clean.yaml --codebook-file examples/codebook.tsv --data-org-file examples/data_org.tsv --emitter-file outputs/clean/emitters.tsv --img-file outputs/clean/images.tif --seed 42
```

> **Note:** If the distribution column is provided in the specified codebook, the number_emitter specified in the CLI or the config file is overridden.
