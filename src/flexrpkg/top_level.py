import os
import time
import pandas as pd
import numpy as np
import sys

def intro_messages():
    ## Intro messages
    print('******************************************************************************')
    print(' ')
    print('Welcome to FLEXR!')
    print('A program for automated multi-conformer model building using unbiased electron')
    print('density map sampling.')
    time.sleep(1)
    print(' ')
    print('Brought to you by the Fischer Lab at St. Jude Children\'s Research Hospital')
    print('Copyright reserved')
    print('Please cite: ')
    print('')
    print('Stachowski, T. R. & Fischer, M.')
    print('FLEXR: automated multi-conformer model building using electron-density map sampling.')
    print('2023. Acta Cryst. D79.')
    print('https://doi.org/10.1107/S2059798323002498')
    print(' ')
    time.sleep(1)
    print('Type -h for info about options')
    print(' ')
    print('******************************************************************************')
    time.sleep(2)
    print('Let\'s get started')
    time.sleep(1)

def get_coot_loc():
    try:
        cootexe = '/opt/homebrew/bin/coot'
        version = os.popen('/opt/homebrew/bin/coot --version').read()
        cootversion = version.split()[3]
        pythonversion = version.split()[8][:4]
        cootloc = '/opt/homebrew/Cellar/coot/'+cootversion+'/lib/python'+pythonversion+'/site-packages/coot/'
        libraryloc = '/opt/homebrew/Cellar/coot/'+cootversion+'/lib/python'+pythonversion+'/site-packages/coot/library/rotamer_library_coot.csv'
    except:
        try:
            cootexe = '/usr/local/bin/coot'
            version = os.popen('/usr/local/bin/coot --version').read()
            cootversion = version.split()[3]
            pythonversion = version.split()[8][:4]
            cootloc = '/usr/local/Cellar/coot/'+cootversion+'/lib/python'+pythonversion+'/site-packages/coot/'
            libraryloc = '/usr/local/Cellar/coot/'+cootversion+'/lib/python'+pythonversion+'/site-packages/coot/library/rotamer_library_coot.csv'
        except:
            cootloc = 'NULL'
            libraryloc = 'NULL'
            cootexe = 'NULL'
    print('FLEXR thinks things are located here:')
    print(libraryloc,cootloc,cootexe)
    return str(libraryloc),str(cootloc),str(cootexe)


def check_library():
    #disable printing pandas warnings
    pd.options.mode.chained_assignment = None

    ## check if rotamer library can be found
    print('Checking dependencies...')
    try:
        ## define location of rotamer library
        print('Loading ideal rotamer library...')
        #library = pd.read_csv('./library/rotamer_library_coot.csv',header=0)
        #library = pd.read_csv('/opt/homebrew/Cellar/coot/1.1.01/lib/python3.11/site-packages/library/rotamer_library_coot.csv',header=0)
        libraryloc,cootloc,cootexe = get_coot_loc()
        library = pd.read_csv(libraryloc,header=0)
        chi_labels = ['chi1_mean','chi2_mean','chi3_mean','chi4_mean']
        for i in chi_labels:
            library[i] = library[i].apply(lambda x : x+360 if x<0 else x)
        return library
    except FileNotFoundError:
        print('Cannot find library.')
        print('Exiting.')
        #sys.exit()
        raise FileNotFoundError()


def args():

    from src.tools import arguments
    from src.tools.arguments import create_parser

    CLI = create_parser()
    try:
        ARGS = CLI.parse_args()
    except:
        ARGS = CLI.parse_args([])

    return ARGS

def test_input_file(filename):

    if filename is None:
        print('No input file defined.')
        print('Exiting...')
        print(' ')
        sys.exit()

def create_log(ARGS):
    with open('log','w') as f:
      for arg in vars(ARGS):
          f.write(str(arg)+' '+str(getattr(ARGS, arg)))
          f.write('\n')
