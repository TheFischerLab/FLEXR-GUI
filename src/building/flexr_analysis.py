"""

Part of FLEXR
Analysis script using Coot

"""
import os
import sys
import pandas as pd
import numpy as np

try:
    import coot
    import coot_utils
except ImportError:
    print('Cannot find Coot.')
    print('Exiting...')
    sys.exit()

def output_summaries(altsfile,imol,flexrmolnum):
    try:
        import matplotlib.pyplot as plt
        from matplotlib_venn import venn2,venn2_circles
    except:
        print('Plotting packages missing.')
        #print('Cannot produce output summary figures.')

    #try:
    # residues with alts
    og = coot_utils.residues_with_alt_confs(imol)
    flex = coot_utils.residues_with_alt_confs(flexrmolnum)

    #format
    og = [x[1] + str(x[2]) for x in og]
    flex = [x[1] + str(x[2]) for x in flex]

    # sets
    og = [str(x) for x in og]
    flex = [str(x) for x in flex]
    og = set(og)
    flex = set(flex)
    subsets = [og,flex]

    # create figure
    figure, (ax1, ax2) = plt.subplots(1, 2,figsize=(10, 4))

    # 1. venn diagram
    v1 = venn2(subsets,('Original','FLEXR'),ax=ax1)
    circles = venn2_circles(subsets=subsets, linestyle='solid',color='black',ax=ax1)
    ax1.title.set_text("Residues modeled \n with alternative conformations")

    alts = pd.read_csv(altsfile).value_counts()
    alts = alts.reset_index()[['res_n','chain']].value_counts().reset_index()
    alts.columns = [['res_n','chain','nalts']]
    alts = alts['nalts'].value_counts()
    alts = alts.reset_index()
    alts.columns = ['nalts','count']

    # 2. bar plot
    ax2.bar(x=alts['nalts'].astype(str),height=alts['count'].astype(int),edgecolor='black',linewidth=1.5)
    #ax2.set_yticks(np.arange(0,max(alts['count']+1),3))
    ax2.set_xlabel('Alternative conformations per residue (n) \n in FLEXR model')
    ax2.set_ylabel('Residues (n)')

    figure.tight_layout()
    figure.savefig('summary.png',dpi=500)

    # 3. array

    matches = list(set(flex).intersection(set(og)))
    flexonly = list(set(flex)-set(og))
    ogonly = list(set(og)-set(flex))

    sorted(matches).reverse()
    sorted(flexonly).reverse()
    sorted(ogonly).reverse()

    #matches = pd.Series(matches,name='Reproduced')
    #flexonly = pd.Series(flexonly,name='FLEXR only')
    #ogonly = pd.Series(ogonly,name='Original only')

    #print(matches)
    #print(flexonly)
    #print(ogonly)

    with open('alts-original-flexr-comparison.txt', 'w+') as f:
        f.write('Reproduced: '+str(matches)+'\n')
        f.write('Original only: '+str(ogonly)+'\n')
        f.write('FLEXR only: '+str(flexonly)+'\n')


    #except:
    #    print('Cannot produce output summary figures.')















