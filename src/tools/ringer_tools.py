# -*- coding: utf-8 -*-
#
#    Tools for analyzing Ringer measurements.
#    Authors: Tim Stachowski & Marcus Fischer
#    Email: tim.stachowski@stjude.org
#    Copyright 2022 St. Jude Children's Research Hospital, Memphis, TN.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
import argparse
import time
import os
import subprocess
import sys
import numpy as np
import pandas as pd

def ringer_parser(file,step):
    """ Organize Ringer outputs from mmtbx """

    step=step
    angles = np.arange(0,360,step=step)

    ## figure out if we're dealing with cryo-em or x-ray ringer output
    dataframe = pd.read_csv(file,header=None)
    if file.endswith('emringer.csv'):
        print('FLEXR Cryo-EM')
        dataframe.columns = ['res','chi','peak']+[*angles]
    else:
        print('FLEXR X-ray')
        dataframe.columns = ['res','map','chi','peak']+[*angles]
        ## map type, future user input option
        dataframe = dataframe[dataframe['map']=='2mFo-DFc']

    ## extract residue number
    res_ns = []
    for j in dataframe['res']:
        res_n = ''.join(i for i in j if i.isdigit())
        res_ns.append(int(res_n))
    dataframe['res_n'] = res_ns

    ## extract chain
    chain = []
    for j in dataframe['res']:
        chain_n = ''.join(i for i in j if not i.isdigit())
        chain_n = chain_n.strip()[-1]
        chain.append(chain_n)
    dataframe['chain'] = chain

    ## extract res_type
    dataframe['res_type'] = dataframe['res'].str[:3]

    return dataframe,angles

def peak_find(file,dataframe,chi,sigma_threshold,plot,height,prominence,width,distance,angles):
    """ Find peaks in Ringer plots """

    try:
        from scipy.signal import find_peaks
    except ImportError:
        print('Please install SciPy')
        sys.exit()

    if plot:
        print('Plotting...')
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            D = False
            print('Matplotlib not found... not plotting')
            sys.exit()
        try:
            os.makedirs(file[:-4]+'_'+chi+'_'+str(sigma_threshold)+'_peaks', exist_ok=True)
            #subprocess.check_output(['mkdir',file[:-4]+'_'+chi+'_'+str(sigma_threshold)+'_peaks'])
        except FileExistsError as e:
            pass

    output = []

    ## loop through residues
    dataframe = dataframe[dataframe['chi']==chi]
    for index,i in dataframe.iterrows():
        res_all = i['res']
        res_n = i['res_n']
        chain = i['chain']
        res = i['res_type']
        chi = chi
        sigma_raw = i[angles].values
        sigma_raw = [float(x) for x in sigma_raw]
        sigma_raw = np.array(sigma_raw)
        dat = np.column_stack((angles,sigma_raw))
        dat = [(x,y) for (x,y) in dat if y >= sigma_threshold]
        sigma_above_threshold = [y for (x,y) in dat]
        angles_cut = [x for (x,y) in dat]
        peaks, properties = find_peaks(sigma_above_threshold,
                                       height=height,
                                       prominence=prominence,
                                       width=width,
                                       distance=distance)
        peak_n = len(peaks)
        area = []
        lefts = []
        rights = []
        for peak in range(peak_n):
            left = properties['left_bases'][peak]
            lefts.append(left)
            right = properties['right_bases'][peak]
            rights.append(right)
            data_to_int = np.column_stack([angles_cut[left:right],
                                           sigma_above_threshold[left:right]])
            angles_int = [j for (j,y) in data_to_int]
            sigma_int = [y for (j,y) in data_to_int]
            area.append(np.trapz(sigma_int,angles_int))
        areas_norm = []
        if len(area) > 1:
            if (lefts[0] == lefts[-1]) or (rights[0] == rights[-1]):
                if area[0] > area[-1]:
                    area[0] = area[0]-area[-1]
                    area[0] = np.absolute(area[0])
                else:
                    area[-1] = area[-1] - area[0]
                    area[-1] = np.absolute(area [-1])
                sum_area = sum(area)
                for areai in area:
                    areas_norm.append(areai/sum_area)
            else:
                sum_area = sum(area)
                for areai in area:
                    areas_norm.append(areai/sum_area)
        if len(area) == 1:
            areas_norm = [1.0]
        if len(area) == 0:
            areas_norm = [0]
        if plot:
            plt.figure(figsize=(5,5))
            plt.plot(angles,sigma_raw,color='black')
            plt.plot(angles_cut,sigma_above_threshold,'red',linestyle='dotted')
            plt.plot([angles_cut[x] for x in peaks],
                     [sigma_threshold]*len(peaks), "x")
            plt.axhline(y=sigma_threshold, linestyle="--", color="gray")
            plt.legend(['raw data',
                        'data used for peak \n finding and integration',
                        'peak','sigma cutoff'],
                        bbox_to_anchor=(1.05, 1),
                        loc='upper left')
            plt.xlabel(chi+'  angle (˚)')
            plt.ylabel('sigma')
            plt.title(res)
            if len(peaks)>0:
                plt.ylim(0-sigma_threshold,np.max(sigma_raw)+sigma_threshold)
            plt.tight_layout()

            plotname = str(res_n)+'_'+chain+'_'+chi+'_'+str(sigma_threshold)+'_peaks.png'
            plt.savefig(plotname,dpi=100)
            plt.close()
            try:
                subprocess.check_output(['mv',plotname,file[:-4]+'_'+chi+'_'+str(sigma_threshold)+'_peaks'])
            except subprocess.CalledProcessError as e:
                continue
        peaks = [angles_cut[x] for x in peaks]
        output.append((res_all,res_n,res,chi,peaks,peak_n,areas_norm,chain))
    output = pd.DataFrame(output)
    output.columns = ['res','res_n','res_type','chi','peak_angles','peaks_n','areas_norm','chain']
    return output


def find_flips(peaks1,peaks2,chi):
    """Find changes (flips) in major/minor conformation of a residue between
       two datafiles based on peak integration """
    dataframe_f1 = peaks1
    dataframe_f2 = peaks2

    if chi == 'chi1':
        good_residues = "Ser, Gln, Asn, Glu, Asp, Arg, Lys, Met, Cys, Leu"
    if chi == 'chi2':
        good_residues = "Gln, Glu, Arg, Lys, Met, Ile"
    if chi == 'chi3':
        good_residues = "Lys, Arg, Met"
    if chi == 'chi4':
        good_residues = "Lys, Arg"
    good_residues = good_residues.split(', ')
    good_residues = [x.upper() for x in good_residues]

    # only consider ringer plots with 2 peaks
    dataframe_f1 = dataframe_f1[dataframe_f1['peaks_n']==2]
    dataframe_f2 = dataframe_f2[dataframe_f2['peaks_n']==2]

    # only consider one chi angle
    dataframe_f1 = dataframe_f1[dataframe_f1['chi']==chi]
    dataframe_f2 = dataframe_f2[dataframe_f2['chi']==chi]

    # only consider residues in allowable chi list
    dataframe_f1['res'] = dataframe_f1['res'].str[:3]
    dataframe_f1 = dataframe_f1[dataframe_f1['res'].isin(good_residues)]
    dataframe_f2['res'] = dataframe_f2['res'].str[:3]
    dataframe_f2 = dataframe_f2[dataframe_f2['res'].isin(good_residues)]

    # format
    dataframe_f1['areas_norm'] = dataframe_f1['areas_norm'].astype(str).str[1:-1]
    dataframe_f1['areas_norm'] = [x.split(',') for x in dataframe_f1['areas_norm']]
    dataframe_f1['areas_norm'] = \
    [(float(x),float(y)) for x,y in dataframe_f1['areas_norm']]
    dataframe_f2['areas_norm'] = dataframe_f2['areas_norm'].astype(str).str[1:-1]
    dataframe_f2['areas_norm'] = [x.split(',') for x in dataframe_f2['areas_norm']]
    dataframe_f2['areas_norm'] = \
    [(float(x),float(y)) for x,y in dataframe_f2['areas_norm']]

    #format
    dataframe_f1['res_n'] = dataframe_f1['res_n'].astype(int)
    dataframe_f2['res_n'] = dataframe_f2['res_n'].astype(int)

    flips_pair = []

    for j in range(len(dataframe_f1)):
        res_n = int(dataframe_f1.iloc[j,0])
        if res_n in dataframe_f2['res_n'].values:
            f1_p1 = dataframe_f1.iloc[j,5][0]
            f2_p1 = dataframe_f2[dataframe_f2['res_n']==res_n]['areas_norm'].values[0][0]
            if f1_p1 > 0.5 > f2_p1:
                print('FLIP! at residue',res_n,chi)
                flips_pair.append((res_n,chi))
            if f1_p1 < 0.5 < f2_p1:
                print('FLIP! at residue',res_n,chi)
                flips_pair.append((res_n,chi))
    flips = pd.DataFrame(flips_pair,columns=['res_n',chi])
    return flips

def find_gainloss(peaks1,peaks2,chi):
    """ Calculate the changes in number of peaks
    (conformations) between two datafiles
    """
    dataframe_f1 = peaks1
    dataframe_f2 = peaks2

    if chi == 'chi1':
        good_residues = "Ser, Gln, Asn, Glu, Asp, Arg, Lys, Met, Cys, Leu"
    if chi == 'chi2':
        good_residues = "Gln, Glu, Arg, Lys, Met, Ile"
    if chi == 'chi3':
        good_residues = "Lys, Arg, Met"
    if chi == 'chi4':
        good_residues = "Lys, Arg"
    good_residues = good_residues.split(', ')
    good_residues = [x.upper() for x in good_residues]

    # only consider one chi angle
    dataframe_f1 = dataframe_f1[dataframe_f1['chi']==chi]
    dataframe_f2 = dataframe_f2[dataframe_f2['chi']==chi]

    # only consider residues in allowable chi list
    dataframe_f1['res'] = dataframe_f1['res'].str[:3]
    dataframe_f1 = dataframe_f1[dataframe_f1['res'].isin(good_residues)]
    dataframe_f2['res'] = dataframe_f2['res'].str[:3]
    dataframe_f2 = dataframe_f2[dataframe_f2['res'].isin(good_residues)]

    dataframe_f1 = dataframe_f1[['res_n','peaks_n','res']]
    dataframe_f1.columns = ['res_n','peak_cryo','res']
    dataframe_f2 = dataframe_f2[['res_n','peaks_n','res']]
    dataframe_f2.columns = ['res_n','peak_rt','res']

    dataframe = pd.merge(dataframe_f1,dataframe_f2,on='res_n',how='left')
    dataframe = dataframe.fillna(0)

    differences = dataframe.peak_cryo.values - dataframe.peak_rt.values
    residues = dataframe['res_n'].values[:]

    dataframe = np.column_stack((residues,differences))

    dataframe = pd.DataFrame(dataframe,columns = ['res','peak_gain_loss'])
    dataframe['res'] = dataframe['res'].astype(int)
    dataframe['chi'] = chi

    return dataframe

def match_and_build(alt,res,library_tmp,rotamer_geometry):
    """ This is for matching peaks found from find_peaks to rotamer geometries """
    chi_labels = ['chi1_mean','chi2_mean','chi3_mean','chi4_mean']
    labels = chi_labels[:len(alt)]
    conf = np.column_stack((labels,alt))
    conf = pd.DataFrame(conf.T)
    conf.columns = conf.iloc[0]
    conf = conf.drop(conf.index[0])
    res_lib = library_tmp[conf.columns]
    match = np.isclose(conf.values.astype(float).astype(int),\
    res_lib.values.astype(int),atol=rotamer_geometry)
    matches = [x.all() for x in match]
    matches = [i for i, x in enumerate(matches) if x]
    if len(matches) != 0:
        ## if multiple matches to that combination of angles, take the most common one
        rot_index = library_tmp.iloc[matches,:].sort_values(by='frequency%')
        #print(matches,rot_index)
        res2build = rot_index['rotamer'].values[0]
        res2buildfreq = rot_index['frequency%'].values[0]
        #print(res,res2build)
    if len(matches) == 0:
        res2build = np.nan
        res2buildfreq = np.nan
    return res2build,res2buildfreq

def pearson_correlation(final_df,ringers,chi):
    """ Produces matrices and a CSV values of pairwise Pearson
    CC calculations per residue.
    Needs an alignment from Ringer-delta
    """

    try:
        import seaborn as sns
    except ImportError:
        print('Please install Seaborn')

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('Please install matplotlib')

    print(' ')
    print('Calculating Pearson CC values....')
    time.sleep(1)

    try:
        subprocess.check_output(['mkdir',chi+'_cc'])
    except subprocess.CalledProcessError as e:
        print('Please delete the directory ',chi+'_cc',' and re-run.')
        sys.exit(1)

    files = pd.read_csv('alignment_new_index.csv',header=0)
    files = files.columns
    output = []
    for i in set(final_df['new_resseq']):
        res_df = []
        LABELS = []
        for j in range(len(ringers)):
            file = files[j]
            if i in ringers[j]['new_resseq'].values:
                LABEL = [x.strip() for x in ringers[j][ringers[j]['new_resseq']==i]\
                .iloc[0,[0,-2,-3]].values]
                LABEL = ' '.join(LABEL)
                LABELS.append(file+','+LABEL)
                tmp = ringers[j][ringers[j]['new_resseq']==i].iloc[0,3:-4]
                res_df.append(tmp)
        res_df = pd.DataFrame(res_df,columns=None)
        res_df = res_df.T
        res_df.columns = LABELS
        corrs = []
        for m in range(len(LABELS)):
            x = res_df.iloc[:,m]
            for n in range(len(LABELS)):
                y = res_df.iloc [:,n]
                corr = x.corr(y,method='pearson')
                corrs.append(corr)
                output.append((LABELS[m].split(',')[0],LABELS[m].split(',')[1],\
                LABELS[n].split(',')[0],LABELS[n].split(',')[1],corr))
        corrs = np.array(corrs)
        matrix = corrs.reshape((len(LABELS),len(LABELS)))
        ax = sns.heatmap(matrix,xticklabels=LABELS,yticklabels=LABELS,vmax=1,\
        annot=True,linewidths=.5)
        ax.set(title='Pearson Correlation Coefficient:   '+str(i))
        plt.tight_layout()
        plt.savefig(os.getcwd()+'/'+chi+'_cc/'+'cc_'+str(i)+'.png',dpi=300)
        plt.close()
    output = pd.DataFrame(output)
    output.columns = ['file1','res1','file2','res2','cc']
    output.to_csv('cc_'+chi+'.csv',header=True,index=False)
    return output
    print('Values saved to: cc_'+chi+'.csv')
    print('Matrix plots saved to: ./'+chi+'_cc')
    print(' ')

def parse_peak_find(filename,sigmathreshold,plot,peakheight,peakprominence,peakwidth,\
                    peakdistance,step,mode):
    ## load data, find peaks, and build single dataframe
    print(filename)
    if mode == 'FLEXRSCORE':
        chis = ['chi0','chi1','chi2','chi3','chi4']
    else:
        chis = ['chi1','chi2','chi3','chi4']
    df = pd.DataFrame([])
    parsed,rot_angles = ringer_parser(filename,step)
    try:
        for chi in chis:
            tmp = peak_find(filename,parsed,chi,sigmathreshold,plot,peakheight,peakprominence,peakwidth,peakdistance,rot_angles)
            df = pd.concat([df,tmp])
        df = df.sort_values(by=['chain','res_n','chi'])
        df = df.drop_duplicates(subset=['res','chi'])
        if '/' in filename:
            file_location = filename.split('/')
            file_out = file_location[-1]
            file_out = '/'.join(file_location[:-1])+'/peak_finder_output_'+str(sigmathreshold)+'_'+file_out
            df.to_csv(file_out,header=True,index=False)
        else:
            df.to_csv('peak_finder_output_'+str(sigmathreshold)+'_'+filename,header=True,index=False)
        return df
    except ValueError as e:
        print('')
        print('Cannot parse Ringer CSV')
        print('')
        if mode == 'FLEXRSCORE':
            print('Ringer CSV missing alanine measurements.')
            print('Re-run EMRinger on your hydrogenated model.')
        else:
            print('Make sure you do not have any alts in your model')
            print('and re-run Ringer.')
        sys.exit(1)

def assemble_matches(df,library,rotamer_geometry):
    print(' ')
    print('Matching Ringer peaks to ideal rotamers...')
    ## loop through all residues an extract info
    ## compare chi angles from Ringer to ideal angles in library
    from itertools import product

    ## assign what residues will always have peaks at certain chi angles
    chi1_allowed = ['SER','ARG','LYS','MET','CYS']
    chi1_branched = ['THR','VAL']
    chi2_branched = ['LEU','ASN','ASP']
    chi2_allowed = ['ILE']
    chi3_branched = ['GLU','GLN']
    rings = ['PHE','TYR','TRP','HIS']

    alt_confs = []
    if len(df)>0:
        for resn,restype,chain in df[['res_n','res_type','chain']].drop_duplicates().values:

            angles = df[(df['res_n']==resn)&(df['res_type']==restype)&(df['chain']==chain)]['peak_angles']
            peakareas = df[(df['res_n']==resn)&(df['res_type']==restype)&(df['chain']==chain)]['areas_norm']
            res = restype
            library_tmp = library[library['res_type']==res]

            angles = [x[:] for x in np.array(angles)]
            angles = [[x for x in sublist] for sublist in angles]

            peakareas = [x[:] for x in np.array(peakareas)]
            peakareas = [[x for x in sublist] for sublist in peakareas]

            if len(angles) > 1:
                alts = list(product(*angles))
                areas = list(product(*peakareas))
            if len(angles) == 1:
                alts = list(angles)
                areas = list(peakareas)
            if len(alts) > 0:
                print(res,angles,resn,chain)
                ## test which residues have certain number of peaks
                ## based on branching vs unbrached vs ring
                if res in chi1_allowed:
                    if len(alts) > 1:
                        for count,k in enumerate(alts):
                            area = sum(areas[count])
                            if len(k) > 1:
                                res2build,res2buildfreq = match_and_build(k,res,library_tmp,rotamer_geometry)
                                alt_confs.append((resn,restype,chain,res2build,res2buildfreq,area))
                    elif (len(alts) == 1) & (restype not in ['SER', 'CYS']):
                        for count,k in enumerate(alts):
                            if len(k) > 1:
                                res2build,res2buildfreq = match_and_build(k,res,library_tmp,rotamer_geometry)
                                alt_confs.append((resn,restype,chain,res2build,res2buildfreq,1))
                    # for serine that only has chi1
                    elif (len(alts) == 1) & (restype in ['SER', 'CYS']):
                        for count,k in enumerate(alts[0]):
                            area = areas[0][count]
                            k = [k]
                            res2build,res2buildfreq = match_and_build(k,res,library_tmp,rotamer_geometry)
                            alt_confs.append((resn,restype,chain,res2build,res2buildfreq,area))
                if res in chi1_branched:
                    if len(alts[0]) > 2:
                        for count,k in enumerate(alts[0]):
                            area = areas[0][count]
                            k = [k]
                            res2build,res2buildfreq = match_and_build(k,res,library_tmp,rotamer_geometry)
                            alt_confs.append((resn,restype,chain,res2build,res2buildfreq,area))
                    elif len(alts[0]) == 2:
                        if res == 'VAL':
                            if alts[0][0] < 100:
                                k = [alts[0][1]]
                            else:
                                k = [alts[0][0]]
                            res2build,res2buildfreq = match_and_build(k,res,library_tmp,rotamer_geometry)
                            alt_confs.append((resn,restype,chain,res2build,res2buildfreq,1))
                        if res == 'THR':
                            if alts[0][0] < 100:
                                k = [alts[0][0]]
                            else:
                                k = [alts[0][1]]
                            res2build,res2buildfreq = match_and_build(k,res,library_tmp,rotamer_geometry)
                            alt_confs.append((resn,restype,chain,res2build,res2buildfreq,1))
                if res in chi2_branched:
                    if (len(angles[0]) > 1) or (len(angles[1]) > 2):
                        for count,k in enumerate(alts):
                            area = sum(areas[count])
                            if len(k) > 1:
                                res2build,res2buildfreq = match_and_build(k,res,library_tmp,rotamer_geometry)
                                alt_confs.append((resn,restype,chain,res2build,res2buildfreq,area))
                    elif (len(angles[0]) == 1) and (len(angles[1]) == 2):
                        k = alts[0]
                        res2build,res2buildfreq = match_and_build(k,res,library_tmp,rotamer_geometry)
                        alt_confs.append((resn,restype,chain,res2build,res2buildfreq,1))
                if res in chi2_allowed:
                    if (len(angles[0]) > 2) or (len(angles[1]) > 1):
                        for count,k in enumerate(alts):
                            area = sum(areas[count])
                            if len(k) > 1:
                                res2build,res2buildfreq = match_and_build(k,res,library_tmp,rotamer_geometry)
                                alt_confs.append((resn,restype,chain,res2build,res2buildfreq,area))
                    elif (len(angles[0]) == 2) and (len(angles[1]) == 1):
                        k = alts[1]
                        res2build,res2buildfreq = match_and_build(k,res,library_tmp,rotamer_geometry)
                        alt_confs.append((resn,restype,chain,res2build,res2buildfreq,1))
                if res in chi3_branched:
                    if len(angles) > 1:
                        if (len(angles[0]) > 1) or (len(angles[1]) > 2) or (len(angles[1]) > 3):
                            for count,k in enumerate(alts):
                                area = sum(areas[count])
                                if len(k) > 1:
                                    res2build,res2buildfreq = match_and_build(k,res,library_tmp,rotamer_geometry)
                                    alt_confs.append((resn,restype,chain,res2build,res2buildfreq,area))
                        elif (len(angles[0]) == 1) and (len(angles[1]) == 1) and (len(angles[2]) == 2):
                            k = alts[1]
                            res2build,res2buildfreq = match_and_build(k,res,library_tmp,rotamer_geometry)
                            alt_confs.append((resn,restype,chain,res2build,res2buildfreq,1))
                    else:
                        for count,k in enumerate(alts):
                            if len(k) > 1:
                                res2build,res2buildfreq = match_and_build(k,res,library_tmp,rotamer_geometry)
                                alt_confs.append((resn,restype,chain,res2build,res2buildfreq,1))
                if res in rings:
                    if len(angles[0]) > 1:
                        for count,k in enumerate(alts):
                            area = areas[count][0]
                            if len(k) > 1:
                                res2build,res2buildfreq = match_and_build(k,res,library_tmp,rotamer_geometry)
                                alt_confs.append((resn,restype,chain,res2build,res2buildfreq,area))
                    elif len(angles[0]) == 1:
                        k = [alts[0][0]]
                        res2build,res2buildfreq = match_and_build(k,res,library_tmp,rotamer_geometry)
                        alt_confs.append((resn,restype,chain,res2build,res2buildfreq,1))
        alt_confs = pd.DataFrame(alt_confs)
        return alt_confs
    else:
      print('No alts found')
      print('Done')

def output_build_list(filename,sigmathreshold,alt_confs,buildlimit,ligand,pdbin,distance,singleconfopt):
    """ organize output, apply some filters"""
    try:
        print(' ')
        alt_confs.columns = ['res_n','res_type','chain','rotamer','frequency%','area_sum']
        alt_confs = alt_confs.replace('',np.nan)
        alt_confs = alt_confs.dropna()
        alt_confs = alt_confs.drop_duplicates(subset=['res_n','res_type','chain','rotamer'])
        ## sort list by frequency, so when we trim we keep the most likely ones
        #alt_confs = alt_confs.sort_values(by=['res_n','chain','frequency%'],ascending=False)
        alt_confs = alt_confs.sort_values(by=['res_n','chain','area_sum'],ascending=False)
        ## continue with residues that n>buildlimit number of confs or include single confs
        alt_confs_limit = pd.DataFrame([])
        if (singleconfopt):
            print('Limiting build list to n<='+str(buildlimit)+' rotamers per residue.')
            print('Also including single conformers.')
            print(' ')
            alt_confs = alt_confs.groupby(['res_n','chain']).head(buildlimit).reset_index(drop=True)
            alt_confs_limit = alt_confs[['res_n','chain']].value_counts().loc[lambda x: x>0].index
            alt_confs = alt_confs[alt_confs[['res_n','chain']].apply(tuple,1).isin(list(alt_confs_limit))]
        else:
            print('Limiting build list to n<='+str(buildlimit)+' rotamers per residue.')
            print(' ')
            alt_confs = alt_confs.groupby(['res_n','chain']).head(buildlimit).reset_index(drop=True)
            alt_confs_limit = alt_confs[['res_n','chain']].value_counts().loc[lambda x: x>1].index
            alt_confs = alt_confs[alt_confs[['res_n','chain']].apply(tuple,1).isin(list(alt_confs_limit))]
            #print(alt_confs.values)
        ## buildlimit building to if residues are near ligand
        if ligand is not None:
            close_residues = find_residues_around_ligand(pdbin,ligand,distance)
            alt_confs = pd.merge(close_residues,alt_confs,on=['res_n','chain'],how='left')
        alt_confs = alt_confs.sort_values(by=['chain','res_n'])
        alt_confs.columns = ['res_n','res_type','chain','rotamer','frequency%','area_sum']
        del alt_confs['frequency%']
        del alt_confs['area_sum']
        alts_path = filename[:-4]+'_'+str(sigmathreshold)+'_alts.csv'
        alt_confs.to_csv(alts_path,header=True,index=False)
        print('Matching rotamers saved to: ringer_alts.csv...')
        print('Ready to build with Coot!!')
        print(' ')
        return alts_path
    except IndexError:
        print('Sorry, no matching rotamers found')
        print('Exiting...')
        print(' ')

def get_atom_info(pdbin):
    """ a function to convert a PDB to a Pandas dataframe. """
    if pdbin.endswith(".pdb"):
        ## radius
        pdb_info = []
        anisou_info = []
        with open(pdbin) as data:
            read_data = data.read()
            lines = read_data.splitlines()
            for line in lines:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atom_type = line[12:16].strip()
                    #remove hydrogens
                    if 'H' not in atom_type[0]:
                    #if True:
                        atom_het = line[0:7].strip()
                        atom_num = int(line[7:12].strip())
                        residue_name = line[17:20].strip()
                        chain = line[21].strip()
                        res_seq = int(line[22:26].strip())
                        alt_loc = line[16:17]
                        # occupancy
                        occ = float(line[56:60].strip())
                        b_fac = float(line[61:66].strip())
                        x,y,z = line[30:38].strip(),line[38:46].strip(),line[46:54].strip()
                        x,y,z = float(x),float(y),float(z)
                        atom_type2 = line[75:80].strip()
                        pdb_info.append((pdbin,atom_het,atom_num,atom_type,alt_loc,occ,residue_name,chain,res_seq,x,y,z,b_fac,atom_type2))
    pdb_info = pd.DataFrame(pdb_info,columns = ['model','atom_het','atom_num','atom_type','alt_loc','occupancy','res_type','chain','res_n','x','y','z','b_factor','atom_type2'])
    return pdb_info

def calculate_lig_protein_distance(pdb_info,distance,lig_res_n,lig_chain):
    """calculate distances of residues from ligand of interest. return array of residues within a certain distance"""
    pdb_info['dist_from_lig']=np.nan
    ## ignore hydrogens
    close_residues = []
    pdb_info = pdb_info[pdb_info['atom_type'].str[0]!='H']
    for chain in pdb_info['chain'].drop_duplicates():
        for resn in pdb_info[(pdb_info['atom_het']!='HETATM')&(pdb_info['chain']==chain)]['res_n'].drop_duplicates():
            res = pdb_info[pdb_info['res_n']==resn]
            lig = pdb_info[(pdb_info['res_n']==lig_res_n)&(pdb_info['chain']==lig_chain)]
            for index,patom in res.iterrows():
                for index,latom in lig.iterrows():
                    lx,ly,lz = latom[['x','y','z']]
                    lcoords = np.array((lx,ly,lz))
                    px,py,pz = patom[['x','y','z']]
                    pcoords = np.array((px,py,pz))
                    dist = np.linalg.norm(lcoords-pcoords)
                    dist = np.absolute(dist)
                    if (dist <= distance):
                        close_residues.append((resn,chain))
    close_residues = pd.DataFrame(close_residues)
    close_residues.columns = ['res_n','chain']
    close_residues = close_residues.drop_duplicates()
    return close_residues

def find_residues_around_ligand(PDB,resid,DI):
    """return list of residues near ligand to run FLEXR on."""
    try:
        print('Finding rotamers to build around ligand: '+resid)
        lig_res_n = int(resid.split(' ')[0])
        lig_chain = resid.split(' ')[1]

    except:
        print('Sorry, input ligand of interest not interpretable.')
        print('Exiting...')
        print(' ')
        sys.exit(1)
    distance = DI
    try:
        pdb_info = get_atom_info(PDB)
    except:
        print('Sorry, PDB file not found.')
        print('Exiting...')
        print(' ')
        sys.exit(1)
    try:
        close_residues = calculate_lig_protein_distance(pdb_info,distance,lig_res_n,lig_chain)
    except:
        print('Sorry, something happened while calculating protein-ligand distances.')
        print('Exiting...')
        print(' ')
        sys.exit(1)
    return close_residues

def render_by_attribute(file_in,file_out,res_t,res_n,attributes):
    """ change values in B-factor column to render in PyMOL """
    print('Adding values to B-factor column in: '+file_out)
    print(' ')
    df = pd.concat([res_t,res_n,attributes],axis=1).reset_index()
    df.columns = ['index','res_t','res_n','attribute']
    df['attribute'] = [round(x,2) for x in df['attribute']]
    with open(file_in,'r') as f_in:
        f_lines = f_in.readlines()
        with open(file_out,'w+') as f_out:
            for line in f_lines:
                if line.startswith('ATOM'):
                    res_t = line[17:20].strip()
                    res_n = line[22:26].strip()
                    value2write = \
                        df[(df['res_t']==res_t)&(df['res_n']==res_n)]\
                            ['attribute'].values
                    if len(value2write) == 1:
                        new_line = line[:60]+' '+str(value2write[0])+' '+\
                            line[69:]
                        f_out.write(new_line)
                    else:
                        new_line = line[:60]+' '+str(0)+'  '+line[67:]
                        f_out.write(new_line)
                if line.startswith('HETATM'):
                    f_out.write(line)
                #else:
                #    f_out.write(line)














