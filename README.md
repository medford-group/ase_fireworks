# ase_fireworks

This is a package implmenting fireworks for an arbitrary ASE calculator. For information on how to set up fireworks please see the [fireworks documentation](https://materialsproject.github.io/fireworks/). Currently implemented fireworks include ASE based optimization and simple "press start" calculations for any calculator. This has not been extensively tested on all ASE calculators. Lattice optimization fireworks have been written in rough form, but are not currently functional.

## Installation

### pip:

You can use pip to directly install this repo:

`pip install git+https://github.com/medford-group/ase_fireworks/tree/master/ase_fireworks`

### Installing manually:

simply clone the repo and run the setup file:

 ```
git clone https://github.com/medford-group/ase_fireworks.git
python setup.py install
```

## usage

To use this package, you will want to make sure you have a fireworks server up and running. Once you have ensured that, you will likely want to use the `get_ase_wflows` function to stack up an arbirary number of runs. This function can take in a single atoms object (for a single run) or a list of atoms objects. Similarly, you can pass in a singel set of input parameters for the calculators or a list. If a single set of input parameters is passed in but a list of atoms is given, it is assumed that you want to use the same parameters for all calculations. You need to pass in the name of the calculator class you want to use and where it needs to be imported from in the `calculator` and `calculator_module` arguments. The input argument `optimizer` is also a switch to turn on optimization. You may choose any optimization method [available in ASE](https://wiki.fysik.dtu.dk/ase/ase/optimize.html). If it is left blank or set to `None`, no optimization will be run. When `optimizer` is not passed in, the calculator is simply run as is with the input parameters provided (i.e. `calculator.calculate()`)

Here is an example:
```
from ase_fireworks import get_ase_wflows
from fireworks import LaunchPad
from ase.build import molecule, bulk

atoms = bulk('Al')

# enter launchpad info
lpad = LaunchPad(
    host= "localhost",
    port= 27017,
    name= "atoms",
    username = "admin",
    password = "admin"

                )

# a .json file containing mongodb information if you intend to write results
# to a mongodb
db_file = './db.json' # a .json file containing mongodb information

wf = get_ase_wflows(
                   # atoms can be a single object or a list of objects
                   atoms,
                   # can be a single set of input parameters or a list
                   parameters = {'rgkmax': 5.0,
                                 'gmaxvr': 0.0,
                                 'mixtype':3
                                },
                   # the name of the calculator class
                   calculator = 'ELK',
                   to_db = True,
                   db_file = db_file,
                   # if optimizer is set, it will run an optimization
                   #optimizer = 'GPMin',
                   #fmax = 0.05,
                   # where to import the calculator from
                   calculator_module = 'ase.calculators.elk',
                   identifiers = None)

lpad.add_wf(wf)
```
