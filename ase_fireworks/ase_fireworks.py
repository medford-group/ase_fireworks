from .mongo import atoms_dict, dict_atoms, json
import os
from fireworks import explicit_serialize, FiretaskBase, FWAction, Firework, Workflow
from ase.atoms import Atoms
from ase.atom import Atom
from .mongo import mongo_atoms_doc, mongo_doc_atoms, MongoDatabase
from pymatgen.io.ase import AseAtomsAdaptor
from collections import OrderedDict
from ase.io.jsonio import encode


@explicit_serialize
class RunCalculatorASE(FiretaskBase):
    """
    Runs the given ASE calculator based on some input dictionary

    Args:
        atoms (Atoms Object): An ase atoms object to run SPARC on
        parameter_dict (Dict): A dictionary of input parameters to run SPARC with
        calculator (str): the name of the calculator to be imported
        to_db (bool): If True, the firework will attmept to connect to a mongodb
                      and store the results there from the calculator.todict
                      function and the atoms object.
        identifier (str/list/dict): a tag to be placed in the mongodb entry to 
                      identify this run from others, can be any type that can go
                      into a mongodb.
        db_file (str): the path to a mongodb json file containing information on
                       how to access the database. Details on this file (and an
                       example) are in the exmaple directory.
        calculator_module (str): the location of the calculator you want to use.
                       The default should be 'ase.calculators'. Otherwise
                       something like 'espresso.espresso' might be what you want
        
    """
    required_params = ['atoms', 'parameter_dict','calculator',
                       'to_db','identifier','db_file','calculator_module']
    optional_params = []
    _fw_name = 'Run ASE calculator'
    def run_task(self,fw_spec):
        # get the calculator
        exec(r'from {} import {} as calculator'.format(self['calculator_module'],
                                                       self['calculator']))
        # set up the class
        parameter_dict = self['parameter_dict']
        atoms = dict_atoms(self['atoms'])
        calc = calculator(**parameter_dict)
        atoms.set_calculator(calc)
        # run the calculation
        calc.calculate(atoms = atoms)
        # store in mongodb
        if self['to_db'] == True:
            from json import load
            db_info = load(open(self['db_file'],'r'))
            id_ = calc.calc_to_database(**db_info )
            if self['identifier'] is not None:
                db = MongoDatabase(**db_info)
                db.modify(id_,{'identifier':self['identifier']})
                print(id_)


class ASE_Run_FW(Firework):
    def __init__(self, atoms, 
                 parameters = {},  
                 calculator = 'VASP',
                 name = 'ASE Run',
                 db_file = None,
                 parents = None,
                 to_db = False,
                 identifier = None,
                 calculator_module = 'ase.calculators',
                 **kwargs
                 ):  
        """
        Runs an SCF calculation in SPARC on the given structure

        Args:
            atoms (Atoms Object or Pymatgen Structure): Input structure
        
            name (str): Name for the Firework
        
            parameters (Dict): the input parameters for SPARC
        
            sparc_command (str): optional, a command used by ase to run SPARC.
                This can also be set with the $ASE_SPARC_COMMAND environment
                variable.

            psuedo_potentials_path (str): optional, a path to the pseudopotentials
                you'd like used. This can also be set using the $PSP_PATH environment
                Variable.

            db_file (str): Not implemented

            parents ([Firework]): Parents of this Firework
        """
    
        t = []
        try:
            translator = AseAtomsAdaptor()
            atoms = translator.get_atoms(atoms)
        except:
            pass
        t.append(RunCalculatorASE(atoms = atoms, parameter_dict = parameters,
             calculator = calculator,
             to_db = to_db,
             db_file = db_file,
             identifier = identifier,
             calculator_module = calculator_module))
        super(ASE_Run_FW, self).__init__(t, parents=parents, name="{}-{}".
                                         format(
                                             dict_atoms(atoms).get_chemical_formula(), name),
                                         **kwargs)

def get_ase_wflows(structures, 
                   parameters = {},             
                   calculator = 'VASP',
                   to_db = True,
                   db_file = None,
                   optimizer = None,
                   fmax = None,
                   identifiers = None,
                   calculator_module = 'ase.calculators',):
    """
    A function to generate an arbitrary number of DFT calculations in a single workflow.
    This is designed to be very simple, no stringing workflows together, just throw in
    the structures and parameters and go. You may pass in a list of structures or a single
    structure. Likewise, you may pass in a single dictionary of parameters of a list 
    corresponding to the list of structures. If a list of structures is passed in with only
    a single parameter dict it is assumed that you want to use the same parameters for all
    the calculations

    inputs:
        structures (ASE Atoms object/list): a single ASE atoms object or a list of atoms
                        objects for the structure(s) you'd like to calculate
        parameters (dict/list): a dictionary of list of dictionaries containing the input
                        arguments for the ASE calculators you're going to run. The list
                        of dictionaries must correspond to the list of structures. If 
                        only one dictionary is provided the same dictionary is used for
                        all structures.
        calculator (str): the name of the calculator to be imported
        to_db (bool): If True, the firework will attmept to connect to a mongodb
                        and store the results there from the calculator.todict
                        function and the atoms object.
        identifier (str/list/dict): a tag to be placed in the mongodb entry to 
                        identify this run from others, can be any type that can go
                        into a mongodb.
        db_file (str): the path to a mongodb json file containing information on
                         how to access the database. Details on this file (and an
                         example) are in the exmaple directory.
        optimizer (str): if this variable is left as None no optimization is performed.
                         Otherwise, input a string containing the ASE optimizer you'd 
                         like to use.
                         (https://wiki.fysik.dtu.dk/ase/ase/optimize.html)
        calculator_module (str): the location of the calculator you want to use.
                         The default should be 'ase.calculators'. Otherwise
                         something like 'espresso.espresso' might be what you want


    returns:
        workflows(list): a list of fireworks workflows
    """
    fws = []

    # check inputs
    if type(structures) != list:
        structures = [structures]
    if type(parameters) != list:  # If no list of parameters is given, use the same for all
        parameters = [parameters] * len(structures) 
    if type(identifiers) != list:
        identifiers = [identifiers]
    if len(parameters) != len(structures):
        raise Exception('The number of parameter dictionaries did not match the number of strucutures')

    # build the workflow from individual fireworks
    if optimizer is not None:
        if fmax is None:
            fmax = 0.05
        # for ASE optimization runs
        for struct, param, identifier in zip(structures, parameters, identifiers):
            name = struct.get_chemical_formula()
            fws.append(ASE_Optimize_FW(atoms_dict(struct),param,
                       calculator = calculator,
                       to_db = to_db,
                       db_file = db_file,
                       identifier = identifier,
                       calculator_module = calculator_module,))
    else:
        # for simple runs
        if fmax is not None:
            Warning('fmax was set, but an optimizer was not chosen, thus no optimization will be performed. To run an optimization, pass in the optimizer argument')
        for struct, param, identifier in zip(structures, parameters, identifiers):
                       fws.append(ASE_Optimize_FW(atoms_dict(struct),param,
                       calculator = calculator,
                       to_db = to_db,
                       db_file = db_file,
                       optimizer = optimizer,
                       identifier = identifier,
                       calculator_module = calculator_module,))

    return Workflow(fws, name="{} calculations wf, e.g.,".format(len(fws)))


@explicit_serialize
class Optimize_Lattice_ASE_FW(FiretaskBase):
    """
    Runs SPARC using ASE based on some input dictionary

    Args:
        atoms (Atoms Object): An ase atoms object to run SPARC on
        parameter_dict (Dict): A dictionary of input parameters to run SPARC with
        
    """
    required_params = ['atoms', 'parameter_dict','calculator',
                       'to_db','identifier','calculator_module']
    optional_params = []

    _fw_name = 'Optimize Lattice ASE'
    def run_task(self,fw_spec):
        # get the calculator
        exec(r'from {} import {} as calculator'.format(self['calculator_module'],
                                                       self['calculator']))
        def Eng(abc, calcargs):
            abc = np.array(abc)
            orig = dict_atoms(self['atoms'].copy())
            new = copy(orig)
            lattice_scale =  np.array([max(abs(a)) for a in atoms.cell])
            new.set_cell(np.multiply(abc/lattice_scale,new.cell.T).T,
                         scale_atoms=True)
            calc = calculator(**calcargs)
            new.set_calculator(calc)
            E = new.get_potential_energy()
            del new, orig
            return E
        
        import numpy as np
        from scipy.optimize import fmin, minimize 
        from copy import copy
        parameter_dict = self['parameter_dict']
        atoms = dict_atoms(self['atoms'])
        x0 = [max(abs(a)) for a in atoms.cell]
        xopt = minimize(Eng, x0, args = (parameter_dict),method='Nelder-Mead',
                        options= {
                                    #'maxiter':150,
                                    'fatol':0.01,'xatol':0.01,
                                    #'adaptive': True
                                    }
                                        )
        x0 = [max(abs(a)) for a in atoms.cell]
        atoms.cell = np.multiply(xopt.x/x0,atoms.cell.T).T
        atoms.set_cell(np.multiply(xopt.x/x0,atoms.cell.T).T,
                       scale_atoms=True)
        calc = SPARC(**parameter_dict)
        atoms.set_calculator(calc)
        calc.calculate()

        if self['to_db'] == True:
            id_ = calc.calc_to_mongo(
                           )
            if self['identifier'] is not None:
                db = MongoDatabase(
                           ) 
                db.modify(id_,{'identifier':self['identifier']})
                print(id_)

class OptimizeLatticeASE(Firework):
    def __init__(self, atoms,
                 parameters = {},
                 name = 'SPARC SCF',
                 sparc_command = None,
                 psuedo_potentials_path = None,
                 db_file = None,
                 parents = None,
                 to_db = False,
                 identifier = None,
                 **kwargs
                 ):

        t = []
        try:
            translator = AseAtomsAdaptor()
            atoms = translator.get_atoms(atoms)
        except:
            pass
        t.append(OptimizeLattice(atoms = atoms, parameter_dict = parameters,
             sparc_command = sparc_command,
             psuedo_potentials_path = psuedo_potentials_path,
             identifier = identifier,
             to_db = to_db))
        super(OptimizeLatticeSPARC, self).__init__(t, parents=parents, name="{}-{}".
                                         format(
                                             dict_atoms(atoms).get_chemical_formula(), name),
                                         **kwargs)

def get_sparc_lattice_optimizations(structures, parameters = {},
                                sparc_command = None,
                                to_db = True,
                                psuedo_potentials_path = None,
                                identifiers = None):
    fws = []
    if type(parameters) != list:  # If no list of parameters is given, use the same for all
        parameters = [parameters] * len(structures)
    if type(identifiers) != list and identifiers is not None:  # If no list of parameters is given, use the same for all
        identifiers = [identifiers] * len(structures)

    for struct, param, identifier in zip(structures, parameters,identifiers):
        name = struct.get_chemical_formula()
        fws.append(OptimizeLatticeSPARC(atoms_dict(struct),param,
                    sparc_command = sparc_command,
                    psuedo_potentials_path = psuedo_potentials_path,
                    identifier = identifier,
                    to_db = to_db))
    return Workflow(fws, name="{} tests wf, e.g.,".format(len(fws)))


@explicit_serialize
class OptimizeASE(FiretaskBase):
    """
    Runs an ASE based optimization using the ASE calculator class chosen

    Args:
        atoms (Atoms Object): An ase atoms object to run SPARC on
        parameter_dict (Dict): A dictionary of input parameters to run SPARC with
        calculator (str): the name of the calculator to be imported
        to_db (bool): If True, the firework will attmept to connect to a mongodb
                      and store the results there from the calculator.todict
                      function and the atoms object.
        identifier (str/list/dict): a tag to be placed in the mongodb entry to 
                      identify this run from others, can be any type that can go
                      into a mongodb.
        db_file (str): the path to a mongodb json file containing information on
                       how to access the database. Details on this file (and an
                       example) are in the exmaple directory.
        optimizer (str): Which optimizer you want to use. You can use any of the
                       optimizers implemented in ASE 
                       (https://wiki.fysik.dtu.dk/ase/ase/optimize.html)
        calculator_module (str): the location of the calculator you want to use.
                       The default should be 'ase.calculators'. Otherwise
                       something like 'espresso.espresso' might be what you want
 
    """
    required_params = ['atoms', 'parameter_dict','calculator',
                       'to_db','identifier','db_file','fmax',
                       'optimizer','calculator_module']
    optional_params = []

    _fw_name = 'Optimize with {}'.format('calculator')
    def run_task(self,fw_spec):
        if self['optimizer'] == None:
            optimizer = 'QuasiNewton'
        else:
            optimizer = self['optimizer']
        exec(r'from ase.optimize import {} as opt'.format(optimizer), globals())
        exec(r'from {} import {} as calculator'.format(self['calculator_module'],
                                                       self['calculator']), globals())

        parameter_dict = self['parameter_dict']
        atoms = dict_atoms(self['atoms'])
        calc = calculator(**parameter_dict)
        atoms.set_calculator(calc)
        relax = opt(atoms,logfile='opt.log',
                    trajectory='opt.traj',
                    restart='opt.pckl')
        relax.run(self['fmax'])
        if self['to_db'] == True:
            from json import load
            db_info = load(open(self['db_file'],'r'))
            id_ = calc.calc_to_database(**db_info )
            if self['identifier'] is not None:
                db = MongoDatabase(**db_info)
                db.modify(id_,{'identifier':self['identifier']})
                print(id_)


class ASE_Optimize_FW(Firework):
    def __init__(self, atoms,
                 parameters = {},
                 calculator = 'VASP',
                 name = 'ASE Optimize',
                 db_file = None,
                 parents = None,
                 to_db = False,
                 identifier = None,
                 optimizer = 'QuasiNewton',
                 fmax = 0.02,
                 calculator_module = 'ase.calculators',
                 **kwargs,
                 ):
        """
        Runs an SCF calculation in SPARC on the given structure

        Args:
            atoms (Atoms Object or Pymatgen Structure): Input structure
        
            name (str): Name for the Firework
        
            parameters (Dict): the input parameters for SPARC
        
            sparc_command (str): optional, a command used by ase to run SPARC.
                This can also be set with the $ASE_SPARC_COMMAND environment
                variable.

            psuedo_potentials_path (str): optional, a path to the pseudopotentials
                you'd like used. This can also be set using the $PSP_PATH environment
                Variable.

            db_file (str): Not implemented

            parents ([Firework]): Parents of this Firework
        """

        t = []
        try:
            translator = AseAtomsAdaptor()
            atoms = translator.get_atoms(atoms)
        except:
            pass
        t.append(OptimizeASE(atoms = atoms, parameter_dict = parameters,
             calculator = calculator,
             to_db = to_db,
             db_file = db_file,
             identifier = identifier,
             optimizer = optimizer,
             fmax = fmax,
             calculator_module = calculator_module))
        super(ASE_Optimize_FW, self).__init__(t, parents=parents, name="{}-{}".
                                         format(
                                             dict_atoms(atoms).get_chemical_formula(), name),
                                         **kwargs)





