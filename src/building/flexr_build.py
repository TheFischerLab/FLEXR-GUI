"""

Part of FLEXR
Script for building in alternative conformations with Coot

"""
import os
import sys
import pandas as pd
from glob import glob


try:
    import coot
    import coot_utils
except ImportError:
    print('Cannot find Coot.')
    print('Not building...')
    print('Exiting...')
    sys.exit()


def gen_atom_list():

    # Possible amino acid residue atoms
    atom_list = ['N','CA','C','O','CB','CG','SD','CE','CD1','CD2','OG1','CG2','OD1','ND2',\
    'ND1','CE1','NE2','CG1','OD2','CD','OG','CE2','CZ','OH','NE','NH1','NH2','OE1','OE2',\
    'NE1','CE3','CZ2','CZ3','CH2','SG','NZ']

    return atom_list

def get_alt_locs():
    # List of possible alt locs - more than you should need
    alt_loc = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P']

    return alt_loc

def build_check(buildopt):
    if buildopt:

        # avoid large backup files
        coot.turn_off_backup(molnum)
        return True
    else:
        print('Building option set to false. Exiting.')
        print('Done.')

def check_files(filein,build_list):
    # If no input file defined: exit
    if (filein is None) or (build_list is None):
        print('Input file not defined...')
        print('Exiting...')
        print('Done.')
        sys.exit()

def building(build_list,filein,branchopt,molnum,exitopt,atom_list,alt_loc):

        #load model
        print('Building alternative conformers...')

        #need to double check
        if branchopt == 'CA':
            branchopt = 0
        else:
            branchopt = 1

        build_list_file = build_list[:-4]

        threshold = build_list.split('_')[-2]
        build_list = pd.read_csv(build_list,header=0)

        molnum = int(molnum)

        if filein == 'None':
            flexrmolnum = coot.copy_molecule(molnum)
            filein = ''
            # doesn't work
            coot.set_molecule_name(flexrmolnum,build_list_file+'_flexr.pdb')
            # delete existing alt confs
            #
        else:
            flexrmolnum = coot.read_pdb(filein)

        for p in build_list[['chain','res_n']].drop_duplicates().values:
            print(p)
            chain=p[0]
            resno=int(p[1])
            for k in range(0,len(build_list[(build_list['res_n']==resno)&\
                                (build_list['chain']==chain)])):
                print('Building alt number: ',k+1)
                rotamer = build_list[build_list['res_n']==resno]['rotamer'].values[k]
                print('For: ',resno,rotamer,chain)
                coot.set_go_to_atom_chain_residue_atom_name(chain,resno,'CA')
                if k == 0:
                    coot.set_residue_to_rotamer_name(flexrmolnum,chain,resno,'','',rotamer)
                if k == 1:
                    coot.altconf()
                    coot.set_add_alt_conf_split_type_number(branchopt)
                    coot_utils.with_auto_accept(\
                    [coot.add_alt_conf_py,flexrmolnum,chain,resno,'','',1])
                    coot.set_residue_to_rotamer_name(flexrmolnum,chain,resno,'','B',rotamer)
                if k > 1:
                    coot.altconf()
                    coot.set_add_alt_conf_split_type_number(branchopt)
                    coot_utils.with_auto_accept(\
                    [coot.add_alt_conf_py,flexrmolnum,chain,resno,'','A',1])
                    for atom in atom_list:
                        coot.set_atom_string_attribute(\
                        flexrmolnum,chain,resno,'',atom,'','alt-conf',alt_loc[k])
                    coot.set_residue_to_rotamer_name(flexrmolnum,chain,resno,'',alt_loc[k],rotamer)
        coot.set_go_to_atom_chain_residue_atom_name(chain,resno+1,'CA')
        #output model
        coot.write_pdb_file(flexrmolnum,build_list_file+'_flexr.pdb')
        coot.write_cif_file(flexrmolnum,build_list_file+'_flexr.cif')

        print('Building finished.')
        print('')

        return flexrmolnum


#def refinement():
        #print('List of residues with altconfsz:')
        #print('')
        #coot_utils.residues_with_alt_confs(molnum)
        #coot_utils.refine_residues_with_alt_conf()

        print('')
def exit(exitopt):

    if (exitopt == 'True'):
        coot.coot_no_state_real_exit(0)
        print('Exiting Coot...')
        print('Done.')

def building_run(build_list,filein,branchopt,molnum,exitopt):

        check_files(filein,build_list)

        atom_list = gen_atom_list()

        alt_loc = get_alt_locs()

        flexrmolnum = building(build_list,filein,branchopt,molnum,exitopt,atom_list,alt_loc)

        exit(exitopt)

        return flexrmolnum

if __name__ == '__main__':
    print(sys.argv[:])
    building_run(sys.argv[3], sys.argv[4], sys.argv[5],sys.argv[6],sys.argv[7])
