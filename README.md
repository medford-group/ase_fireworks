# ase_fireworks

This is a package implmenting fireworks for an arbitrary ASE calculator. For information on how to set up fireworks please see the [fireworks documentation](https://materialsproject.github.io/fireworks/). Currently implemented fireworks include ASE based optimization and simple "press start" calculations for any calculator. This has not been extensively tested on all ASE calculators.

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

wf = get_ase_wflows(atoms,
                   parameters = {'rgkmax': 5.0,
                                 'gmaxvr': 0.0,
                                 'mixtype':3
                                },
                   calculator = 'ELK',
                   to_db = True,
                   db_file = db_file,
                   #optimizer = 'GPMin',
                   #fmax = 0.05,
                   calculator_module = 'ase.calculators.elk',
                   identifiers = None)

lpad.add_wf(wf)
```
