# FISHsim (Multiplexed Spots Simulator)

This is a simulation software that generates photorealistic images of flourescent probes.

## Installation
We recommend using conda to set up a virtual environment. 

```
conda env create -f environments.yml
```



## Configuring the Program Parameters

Update `src/config.yml` to set the parameters for the simulator.

Some parameters can also be set using the Command Line Interface. This overrides the settings defined in the `config.yml` file for those specific parameters


## Program Usage

## Running as an installed package
If you have installed fishsim as a package, then

1. Update the parameters in resources\configs\config.yml

2. Run in command line to simulate images. For dataset-name, no need to include 'merfish_' in the beginning
```
fishsim run_merfish --config-file .\resources\config_serval_experiment_default.yml --output-dir-name serval_experiment_default
```
3. If you wish to subdivide. For dataset-name, no need to include 'merfish_' in the beginning
```
fishsim subdivide --dataset-name '' --total-fov 10
```



> **Note:** If the distribution column is provided in the specified codebook, the number_emitter specified in the CLI or the config file is overridden.


## Running the Unit Tests

Run *pytest* in the `src` directory the run the unit tests associated with the program
