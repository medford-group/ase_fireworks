import os
import datetime
from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from fireworks.core.fworker import FWorker
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen import Structure
from ase_fireworks import get_ase_wflows
import yaml
import json




base_dir = '/gpfs/pace1/project/chbe-medford/medford-share/users/bcomer3/pseudopotential_generation/'

proc_nums = {
             'dcdft':10
            }

pkg_info = {
            'sparc':{'calculator':'SPARC','calculator_module':'sparc.sparc'},
            'abinit':{'calculator':'Abinit','calculator_module':'ase.calculators.abinit'},
            'espresso':{'calculator':'Espresso','calculator_module':'espresso.espresso'},
            'elk':{'calculator':'ELK','calculator_module':'ase.calculators.elk'},
            }

software_packages = ['abinit']


systems = ['dcdft']

def populate_launchpad(software, systems, optimizer = None):
    """
    A simple function to fill a workflow with a set of systems
    """
    # load in fireworks
    launch_pad = yaml.load(open('../config/my_launchpad.yaml', 'r'))

    # this is messy, but it has to be done
    del launch_pad['ssl_ca_file']
    del launch_pad['strm_lvl']
    del launch_pad['user_indices']
    del launch_pad['wf_user_indices']

    lpad = LaunchPad(**launch_pad)

    # set up Abinit's input settings

    db_file = os.getcwd() + '/../config/db.json'
    for system_class in systems:
        # load in the json file
        systems = json.load(open('{}{}.json'.format(base_dir+'/staging/structures/', system_class), 'rb'))
        parameters = json.load(open('{}{}.json'.format(base_dir+'/staging/parameters/', system_class),'rb'))
        # reformat into lists
        ids = []
        systems_list = []
        parameters_list = []
        for id_, system in systems.items():
            systems_list.append(system)
            parameters_list.append(parameters[id_])
            ids.append(id_)
        # convert from pymatgen structures to ase atoms objects
        systems_list = [AseAtomsAdaptor.get_atoms(Structure.from_dict(a)) for a in systems_list]

        wf = get_ase_wflows(systems_list,
                   parameters = parameters,
                   calculator = pkg_info[software]['calculator'],
                   to_db = True,
                   db_file = db_file,
                   optimizer = optimizer,
                   calculator_module = pkg_info[software]['calculator_module'],
                   identifiers = None)
        # add the workflow
        lpad.add_wf(wf)


os.chdir(base_dir)
for package in software_packages:
    for procs, launches in proc_nums.items():
        os.chdir('{}/{}/runs'.format(procs,package))
        os.system('lpad -c ../config/ reset --password {}'\
                  .format(datetime.datetime.now().strftime('%Y-%m-%d')))

        populate_launchpad(package, systems)

        #os.system('qlaunch -c ../config rapidfire --nlaunches {}'.format(launches))
        os.chdir(base_dir)


