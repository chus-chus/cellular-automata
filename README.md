# Cellular automata simulator

##### Simulator for the study of genetically stable and mutator cells in a cancer setting

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

A handy way of launching experiments is by speciying the parameters
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
set to a default.