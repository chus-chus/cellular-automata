# Cellular automata simulator

### Simulator for the study of genetically stable and mutator cell populations in a cancer setting

<img src="https://github.com/chus-chus/cellular-automata/raw/master/example_figures/header_system_evolution.gif"/>

## Installation

Just clone or download the repo! Create a new `Python` environment and / or make sure you have
`numpy` and `matplotlib` installed.

(Tested on `Python 3.8`)

## Usage

Navigate to the downloaded folder. Go to the `src` folder. From here, you can run the simulations
with the `main.py` script:

```
python3 main.py
```

Get help on the parameters by calling

``` 
python3 main.py --help
```

You can change any parameter in either the command line or specify parameters for particular 
experiments in `.txt` files. For example, changing parameters from the command line:

```
python3 main.py --epochs 150 --mutatorRR 0.4
```

A handy way of launching experiments is by specifying the parameters
in `.txt` files. These need to be in the `experiment_templates` folder. They
have to follow this format (changing the same params as in the prev. example):

``` 
--epochs
150
--mutatorRR
0.4
```

Now, you can launch experiments specifying the files as argument:

``` 
python3 main.py @../experiment_templates/file.txt
```

In both ways of specifying parameters, the ones that are not specified will be
set to a default. Check the help to see default values.

### Simulation output

For each simulation, a folder in the directory `experiment_results` will be created (this one will also be created if
it does not exist). The folder can be identified by the datetime of the experiment.
It contains figures regarding the evolution of cell types and their states, as well as an animation of the 
population evolution (only if the param `createAnimation` is `True`). It also contains a `.txt` file with 
the params. used to run the experiment (with the correct format), making replication even easier.

## Example

We run a simulation with the following params:

```
--epochs
100
--mutatorRR
5
--stableRR
1
--stableRepairProb
0.3
--mutatorRepairProb
0.2
--createAnimation
True
--randomSeed
1
```

That is, a simulation with 100 epochs, where the replication rates and repair probability is different for mutator
and stable cells. Other parameters are set to default. This outputs the following figures in the corresponding folder 
of `experiment_results/`:

<img src="https://github.com/chus-chus/cellular-automata/raw/master/example_figures/system_evolution.gif" />
<img src="https://github.com/chus-chus/cellular-automata/raw/master/example_figures/total_population_evolution.png" />

<img src="https://github.com/chus-chus/cellular-automata/raw/master/example_figures/healthy_population_evolution.png" />
<img src="https://github.com/chus-chus/cellular-automata/raw/master/example_figures/damaged_population_evolution.png" />
<img src="https://github.com/chus-chus/cellular-automata/raw/master/example_figures/mutated_population_evolution.png" />

## Acknowledgements
This work has been created for the subject *Mathematical Models in Biology* of the Master of Mathematics @ UPC. The corresponding report
detailing the simulation models can be found in this repo.