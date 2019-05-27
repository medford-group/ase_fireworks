import os
from tempfile import mkstemp
from shutil import move
from os import remove, close


dbs = ['dcdft']
envs = {'sparc':'module load anaconda3/4.2.0;source activate atm',
        'espresso':'module load use.own;module load anaconda/2-4.2.0;source activate dev_esp',
        'abinit':'module load anaconda3/4.2.0;source activate atm',
}

n_p = {
    'dcdft':[1,8],
      }

def replace(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    with open(abs_path,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    close(fh)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)

cur_dir = os.getcwd()

for db in dbs:
    os.system('cp -r example ' + db)
    os.chdir(db)
    base_dir = os.getcwd()
    print(base_dir)
    for software in ['espresso', 'sparc', 'abinit']:
        replace(software + '/config/db.json', 'COLLECTION', db)
        replace(software + '/config/my_launchpad.yaml', 'NAME', software + '_' + db) 
        replace(software + '/config/my_launchpad.yaml', 'LOGPATH',
                base_dir + '/' + software)
        replace(software + '/config/my_qadapter.yaml', 'SOFTPATH',
                base_dir + '/' + software)
        replace(software + '/config/my_qadapter.yaml','ENV_SETUP', envs[software])
        replace(software + '/config/FW_config.yaml','CONFIGPATH', '{}/{}/config'.format(base_dir,software))

        replace(software + '/config/my_qadapter.yaml',\
                    'NODES', str(n_p[db][0]))
        replace(software + '/config/my_qadapter.yaml',\
                    'PROCS', str(n_p[db][1]))

                    
                    
    os.chdir(cur_dir)
