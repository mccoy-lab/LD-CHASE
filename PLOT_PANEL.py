#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PLOT_PANEL

Plots log-likelihood ratio vs. chromosomal position from a LLR file.

May 20, 2020
"""

import pickle, bz2, gzip, collections, math
from statistics import mean, variance
from math import log
from operator import attrgetter, itemgetter
import argparse, sys
from itertools import product, starmap, chain, accumulate, islice
from collections import ChainMap

def chr_length(chr_id):
    """ Return the chromosome length for a given chromosome, based on the reference genome hg38.""" 
    #The data of chromosome length was taken from https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38
    length_dict = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
                  'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422,
                  'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
                  'chr16': 90338345, 'chr17':  83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
                  'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415}
    return length_dict[chr_id]

def mean_and_var(x):
    """ Calculates the mean and variance. """
    cache = tuple(x)
    m = mean(cache)
    var = variance(cache, xbar=m)
    return m, var

def mean_and_std_of_mean_of_rnd_var(A):
    """ Calculates the mean and sample standard deviation of the mean of random variables.
        Each row of A represents a random variable, with observations in the columns."""
    if type(A)==dict:
        A = tuple(tuple(i) for i in A.values()) 
    
    M, N = len(A), len(A[0])  ### N is the number of random variables, while M is the number of samples.
    mu = sum(sum(likelihoods_in_window) for likelihoods_in_window in A) / N
    arg = ((sum(sampled_likelihoods) - mu)**2 for sampled_likelihoods in zip(*A))
    std = (sum(arg) / (N - 1))**.5 / M
    mean = mu / M
    return mean, std

def LLR(y,x):
    """ Calculates the logarithm of y over x and deals with edge cases. """
    if x and y:
        result = log(y/x)
    elif x and not y:
        result = -1.23456789 
    elif not x and y:
        result = +1.23456789 
    elif not x and not y:
        result = 0 
    else:
        result = None    
    return result

def load_likelihoods(filename):
    """ Loads from a file a dictionary that lists genomic windows that contain
    at least two reads and gives the bootstrap distribution of the 
    log-likelihood ratios (LLRs). """
    
    Open = {'bz2': bz2.open, 'gzip': gzip.open}.get(filename.rpartition('.')[-1], open)
    
    with Open(filename, 'rb') as f:
        likelihoods = pickle.load(f)
        info = pickle.load(f)
    return likelihoods, info

def show_info(info):
    S = info['statistics']
    ancestral_makeup = ", ".join("{:.1f}% {}".format(100*v, k) for k, v in info['ancestral_makeup'].items()) if type(info['ancestral_makeup'])==dict else ', '.join(info['ancestral_makeup'])
    matched_alleles = ", ".join("{}: {:.1f}%".format(k,100*v) for k, v in info['statistics']['matched_alleles'].items())

    print('\nFilename of the disomy observation table: %s' % info['disomy_obs_filename'])
    print('\nFilename of the monosomy observation table: %s' % info['monosomy_obs_filename'])
    print('\nSummary statistics:')
    print('-------------------')
    print('Chromosome ID: %s' % info['chr_id'])
    print('Depth of coverage of the disomy sequence: %.2f' % info['depth']['disomy'])
    print('Depth of coverage of the monosomy sequence: %.2f' % info['depth']['monosomy'])
    print('Number of genomic windows: %d, Mean and standard error of genomic window size: %d, %d.' % (S.get('num_of_windows',0),S.get('window_size_mean',0),S.get('window_size_std',0)))
    print('Mean and standard error of meaningful reads per genomic window from the disomy sequence: %.1f, %.1f.' % (S.get('disomy_reads_mean',0), S.get('disomy_reads_std',0)))
    print('Mean and standard error of meaningful reads per genomic window from the monosomy sequence: %.1f, %.1f.' % (S.get('monosomy_reads_mean',0), S.get('monosomy_reads_std',0)))
    print('Ancestral makeup: %s, Fraction of alleles matched to the reference panel: %s.' % (ancestral_makeup, matched_alleles))

    if S.get('LLRs_per_chromosome',None):
        L = S['LLRs_per_chromosome']
        print("--- Chromosome-wide LLR between UNMATCHED and MATCHED ----")
        print(f"Mean LLR: {L['mean_of_mean']:.3f}, Standard error of the mean LLR: {L['std_of_mean']:.3f}")
        print(f"Fraction of genomic windows with a negative LLR: {L['fraction_of_negative_LLRs']:.3f}")
        
    

def bin_genomic_windows(windows,chr_id,num_of_bins):
    """ Lists the bins and gives the genomic windows that they contain. """
    bin_size = chr_length(chr_id) / num_of_bins
    result = {}
    j = 0
    
    for i in range(num_of_bins): ### All bins before the first the genomic window are filled with Nones.
        if sum(windows[0])/2 < (i+1)*bin_size:
            break
        result[i/num_of_bins,(i+1)/num_of_bins] = None
    
    for k,(a,b) in enumerate(windows):
        if not bin_size*i <= (a+b)/2 < bin_size*(i+1):
            result[i/num_of_bins,(i+1)/num_of_bins] = (j,k)
            j = k
            for i in range(i+1,num_of_bins): #Proceed to the next non-empty bin; Empty bins are filled with Nones.
                if (a+b)/2 < (i+1)*bin_size:
                    break
                result[i/num_of_bins,(i+1)/num_of_bins] = None
    
    for i in range(i,num_of_bins): ### All bins after the last the genomic window are filled with Nones.
        result[i/num_of_bins,(i+1)/num_of_bins] = (j,k) if j != k else None
        j = k 
    return result

def binning(LLRs_per_window,info,num_of_bins):
    """ Genomic windows are distributed into bins. The LLRs in a genomic windows
    are regarded as samples of a random variable. Within each bin, we calculate
    the mean and population standard deviation of the mean of random variables. 
    The boundaries of the bins as well as the mean LLR and the standard-error
    per bin are returned. """
             
    #K,M,V = tuple(LLR_stat.keys()), *zip(*LLR_stat.values())
    list_of_windows = [*LLRs_per_window.keys()]
    bins = bin_genomic_windows(list_of_windows, info['chr_id'], num_of_bins)
    X = [*bins]
    
    LLR_matrix = [*LLRs_per_window.values()]
    Y, E = [], []
    for C in bins.values():
        if C:
            mean, std = mean_and_std_of_mean_of_rnd_var(LLR_matrix[C[0]:C[1]])
        else:
            mean, std = None, None
        
        Y.append(mean)
        E.append(std)
    
    return X,Y,E

def detect_transition(X, Y, E):
    """  Traces meiotic crossovers based on inferred switches between
         tracts of BPH and SPH trisomy. """
    
    A = [(l,j,k) for i,j,k in zip(X,Y,E) for l in i if j!=None and (abs(j)-k)>0 ]
    if len(A)>=4:
        x,y,e = zip(*A)
        result = [.5*(x0+x1) for x0,x1,y0,y1 in zip(x[1::2],x[2::2],y[1::2],y[2::2]) if y1/y0<0] 
    else:
        result = []
        
    return result

def detect_crossovers(genomic_windows, mean_of_LLRs, variance_of_LLRs, z_score=1.96, lookahead=20):
    """ Detecting crossovers by indetifying transitions between BPH and SPH 
        regions. """
        
    crossovers = {}
    #Scan the chromosome in the 5'-to-3' direction to find crossovers.
    x_coord = tuple(0.5*(a+b) for a,b in genomic_windows)
    acc_means = tuple(accumulate(mean_of_LLRs))
    acc_vars = tuple(accumulate(variance_of_LLRs))
    triple = tuple(zip(x_coord, acc_means, acc_vars))
    
    #maxima and minima candidates are temporarily stored in mx and mn, respectively.
    mx, mn, last_ind = None, None, 0
        
    for index, (x, y, v) in enumerate(triple):
        
        if  mx==None or y > mx:
            mx_index,mx_pos,mx,mx_var = index,x,y,v
        
        if  mn==None or y < mn:
            mn_index,mn_pos,mn,mn_var = index,x,y,v

        if mx!=None and 0 < (mx-y)-z_score*(v-mx_var)**.5 and index-mx_index>=lookahead:
            for x2, y2, v2 in triple[max(mx_index-lookahead,last_ind):last_ind:-1]:
                if 0 < (mx-y2)-z_score*(mx_var-v2)**.5:
                    kappa = min((mx-y2)/(mx_var-v2)**.5,(mx-y)/(v-mx_var)**.5)
                    crossovers[mx_pos] = kappa
                    mx, mn, last_ind = None, None, mx_index # set algorithm to find the next minima
                    break
                
        if mn!=None and 0 < (y-mn)-z_score*(v-mn_var)**.5  and index-mn_index>=lookahead:    
            for x2, y2, v2 in triple[max(mn_index-lookahead,last_ind):last_ind:-1]:
                if  0 < (y2-mn)-z_score*(mn_var-v2)**.5:
                    kappa = min((y2-mn)/(mn_var-v2)**.5,(y-mn)/(v-mn_var)**.5)
                    crossovers[mn_pos] = kappa
                    mx, mn, last_ind = None, None, mn_index # set algorithm to find the next maxima
                    break
                
    return crossovers

def detect_crossovers_v2(genomic_windows, mean_of_LLRs, variance_of_LLRs, z_score=1.96, lookahead=20):
    """ Detecting crossovers by indetifying transitions between BPH and SPH 
        regions. """
        
    crossovers = {}
    #Scan the chromosome in the 5'-to-3' direction to find crossovers.
    x_coord = tuple(0.5*(a+b) for a,b in genomic_windows)
    acc_means = tuple(accumulate(mean_of_LLRs))
    acc_vars = tuple(accumulate(variance_of_LLRs))
    triple = tuple(zip(x_coord, acc_means, acc_vars))
    
    #Maxima and minima candidates are temporarily stored in mx and mn, respectively.
    mx, mn, last_ind, last_extremum = None, None, 0, 0
        
    def recover_skipped_extremum(extremum_type,last_ind,ind):
        """ Recovers skipped extremum """
        ### extremum_typ is the max (minimum) function if a maximum (minimum) was skipped.
        ### last_ind is the index of the last detected extremum.
        ### ind is the index of the most recent detected extremum.
        M0, M2, V0, V2 = acc_means[last_ind],acc_means[ind], acc_vars[last_ind], acc_vars[ind]
        Z1, X1, M1, V1 = extremum_type((((M1-M0)/(V1-V0)**.5-(M2-M1)/(V2-V1)**.5, X1, M1, V1) 
                              for X1,M1,V1 in triple[last_ind+5:ind-4] if V0!=V1!=V2), key=itemgetter(0))
        kappa = min(abs(M2-M1)/(V2-V1)**.5,abs(M0-M1)/(V1-V0)**.5)
        print(f'Recovering skipped {extremum_type.__name__:s} point:',(Z1, X1, M1, V1))
        return {X1: kappa}
        
    for index, (x, y, v) in enumerate(triple):
        
        if  mx==None or y > mx:
            mx_index,mx_pos,mx,mx_var = index,x,y,v
        
        if  mn==None or y < mn:
            mn_index,mn_pos,mn,mn_var = index,x,y,v

        if mx!=None and 0 < (mx-y)-z_score*(v-mx_var)**.5 and index-mx_index>=lookahead: #maximum point detected
            for x2, y2, v2 in triple[max(mx_index-lookahead,last_ind):last_ind:-1]:
                if 0 < (mx-y2)-z_score*(mx_var-v2)**.5:
                   
                    if last_extremum == +1: #Recovering skipped minimum point.
                        crossovers.update(recover_skipped_extremum(min,last_ind,mx_index)) 

                    kappa = min((mx-y2)/(mx_var-v2)**.5,(mx-y)/(v-mx_var)**.5)
                    crossovers[mx_pos] = kappa
                    mx, mn, last_ind, last_extremum = None, None, mx_index, +1 # set algorithm to find the next minimum
                    break
                
        if mn!=None and 0 < (y-mn)-z_score*(v-mn_var)**.5  and index-mn_index>=lookahead: #minimum point detected
            for x2, y2, v2 in triple[max(mn_index-lookahead,last_ind):last_ind:-1]:
                if  0 < (y2-mn)-z_score*(mn_var-v2)**.5:
                    
                    if last_extremum == -1: #Recovering skipped maximum point.
                        crossovers.update(recover_skipped_extremum(max,last_ind,mn_index))
                        
                    kappa = min((y2-mn)/(mn_var-v2)**.5,(y-mn)/(v-mn_var)**.5)
                    crossovers[mn_pos] = kappa
                    
                    mx, mn, last_ind, last_extremum = None, None, mn_index, -1 # set algorithm to find the next maxima
                    break
                
    return crossovers
    
def clusters_crossovers(crossovers, dx, min_frac=0.5):
    """ Clusters crossovers in a region of width dx, when at least min_frac of 
        the samples admit this crossover. Then, the crossovers in each cluster
        are averaged and associated with a single crossover in the monosomy. """
    
    monosomy_crossovers = {}
    
    C_filtered = [*filter(len,crossovers.values())] # Ignore empty cases.
    lenC = len(C_filtered)
    
    if lenC:
        C_flatten, K_flatten = zip(*sorted(ChainMap(*C_filtered).items()) ) # Combine all crossovers into one sorted list.
       
    for k in range(lenC, int(min_frac * lenC)+1, -1):
        i = 0
        while(i<=len(C_flatten)-k):
            if C_flatten[i+k-1]-C_flatten[i] < dx:
                monosomy_crossovers[sum(C_flatten[i:i+k])/k] = (C_flatten[i:i+k], k/lenC, min(K_flatten[i:i+k])  )
                C_flatten = C_flatten[:i] + C_flatten[i+k:] #Equivalent of del C_flatten[i:i+k+1], but works on tuples.
                K_flatten = K_flatten[:i] + K_flatten[i+k:] #Equivalent of del K_flatten[i:i+k+1], but works on tuples.
                i += k
            else:        
                i += 1    
            
    # The keys are the averaged crossovers of each cluster. In addition, the value is a tuple that contains:
    # (a) crossovers in the cluster, (b) proportion crossovers that supported the cluster formation and
    # (c) The minimal kappa score of the supporting crossovers.
    return monosomy_crossovers

def clusters_crossovers_v2(crossovers, dx, min_frac=0.5):
    """ Clusters crossovers in a region of width dx, when at least min_frac of 
        the samples admit this crossover. Then, the crossovers in each cluster
        are averaged and associated with a single crossover in the monosomy. """
    
    ### disomy_id : position : kappa
    monosomy_crossovers = {}
    
    lenC = sum(1 for disomy_id, events in crossovers.items() if len(events))
    
    if lenC:
    
        crossovers_flatten = ((disomy_id,position,kappa) 
                                  for disomy_id, events in crossovers.items() 
                                      for position,kappa in events.items())
        
        D_flatten, C_flatten, K_flatten = zip(*sorted(crossovers_flatten,key=itemgetter(1)))
           
    
    for k in range(2*lenC, int(min_frac * lenC)+1, -1):
        i = 0
        while(i<=len(C_flatten)-k): #i is the index of the sliding window
            if C_flatten[i+k-1]-C_flatten[i] < dx:
                SUPPORTING_INFO = (*zip(islice(D_flatten,i,i+k), islice(C_flatten,i,i+k), islice(K_flatten,i,i+k)),)
                monosomy_crossovers[sum(C_flatten[i:i+k])/k] = (SUPPORTING_INFO, k/lenC, min(K_flatten[i:i+k])  )

                
                D_flatten = D_flatten[:i] + D_flatten[i+k:] #Equivalent of del D_flatten[i:i+k+1], but works on tuples.
                C_flatten = C_flatten[:i] + C_flatten[i+k:] #Equivalent of del C_flatten[i:i+k+1], but works on tuples.
                K_flatten = K_flatten[:i] + K_flatten[i+k:] #Equivalent of del K_flatten[i:i+k+1], but works on tuples.
               
                i += k #overlaps between crossover clusters are not allowed.
            else:        
                i += 1    
            
    # The keys are the averaged crossovers of each cluster. In addition, the value is a tuple that contains:
    # (a) crossovers in the cluster, (b) proportion crossovers that supported the cluster formation and
    # (c) The minimal kappa score of the supporting crossovers.
        
    return monosomy_crossovers

def capitalize(x):
    return x[0].upper() + x[1:]
    
def panel_plot(DATA,**kwargs):
    """ Creates a multi-panel figure. For each numbered chromosome, a figure 
        depicts the log-likelihood ratio vs. chromosomal position for BPH over
        SPH. """
    
    import matplotlib as mpl

    export = dict()
    

    scale = kwargs.get('scale', 0.5)
    bin_size = kwargs.get('bin_size', 4000000)
    z_score = kwargs.get('z_score', 1.96)
    lookahead = kwargs.get('lookahead', 25)

    save = kwargs.get('save', '')
    extension = kwargs.get('extension', 'svg')

    
    fs=28 * scale
    columns = 6
    rows = math.ceil(len(DATA)/columns)
    
    
    if save!='':
            mpl.use('Agg')
    else:
        #['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 'WXCairo', 'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']
        mpl.use('Qt5Agg')
    mpl.rcParams.update({'figure.max_open_warning': 0})
    import matplotlib.pyplot as plt
    
    num_of_bins = {'chr'+str(i): chr_length('chr'+str(i))//bin_size for i in [*range(1,23)]+['X','Y']}

    colors = {'purple': (177/255,122/255,162/255),
              'orange': (242/255,142/255,44/255),
              'red': (239/255,106/255,92/255),
              'blue': (104/255,162/255,183/255),
              'green':(104/255,162/255,104/255)}

    if len(DATA)>columns:
        fig,axs = plt.subplots(rows ,columns, sharex='col', sharey='row', figsize=(6.666 * columns * scale, 5.625 * rows * scale))
        fig.subplots_adjust(left=0.045, bottom=0.13-rows*0.015, right=.99, top=(0.81+rows*0.025 if kwargs.get('title',None) else 0.86+rows*0.025), wspace=None, hspace=None)
    else:
        fig,axs = plt.subplots(rows ,columns, sharex='none', sharey='row', figsize=( 6.666 * columns * scale, 1.25 * 5.625 * rows * scale))
        fig.subplots_adjust(left=0.05, bottom=0.3, right=.99, top=(0.82 if kwargs.get('title',None) else 0.86), wspace=None, hspace=None)
    
    
    AX = [i for j in axs for i in j] if len(DATA)>columns else axs
        
    YMAX = [0]*len(DATA)
    crossovers = {}
    for g,(ax1,(identifier,(likelihoods,info))) in enumerate(zip(AX,DATA.items())):


        LLRs = {window: tuple(starmap(LLR,likelihoods_in_window))
                           for window,likelihoods_in_window in likelihoods.items()} 
                                                            
        
        X,Y,E = binning(LLRs,info,num_of_bins[info['chr_id']])
        Y = [(y if y else 0) for y in Y]
        E = [(z_score*e if e else 0) for e in E]
        
        export[identifier] = {'x': X, 'y': Y, 'z_score * SD': E}

        
        T = [(x[1]+x[0])/2 for x in X]                
        steps_x = [X[0][0]]+[i[1] for i in X[:-1] for j in (1,2)]+[X[-1][1]]
        steps_y = [i for i in Y for j in (1,2)]
        ax1.plot(steps_x, steps_y, label=identifier ,color=colors['green'], linewidth=2, zorder=10, scalex=True, scaley=True, alpha=0.8)
        
        P = [(x[1]-x[0])/2 for x in X]                
        ax1.errorbar(T, Y, xerr = P, ecolor=colors['green'],marker=None, ls='none',alpha=1, zorder=13, linewidth=5*scale) 
        ax1.errorbar(T, Y, yerr = E, ecolor='black',marker=None, ls='none',alpha=0.2, zorder=15, linewidth=4*scale) 
        
        yabsmax = max(map(abs,Y))
        
        ###transitions.append(detect_transition(X,Y,E))
        
        l = chr_length(info['chr_id'])
        genomic_windows = info['statistics']['LLRs_per_genomic_window']
        mean_of_LLRs = [a for a,b in info['statistics']['LLRs_per_genomic_window'].values()]
        variance_of_LLRs = [b for a,b in info['statistics']['LLRs_per_genomic_window'].values()]
        unnormalized_crossovers = detect_crossovers_v2(genomic_windows, mean_of_LLRs, variance_of_LLRs, z_score=z_score, lookahead=lookahead)
        crossovers[g] = {pos/l: kappa for pos,kappa in unnormalized_crossovers.items()}
        export[identifier].update({'crossovers': [*crossovers[g]]})

        YMAX[g] = yabsmax if YMAX[g]< yabsmax else YMAX[g]

    for g,(ax1,(identifier,(likelihoods,info))) in enumerate(zip(AX,DATA.items())):
        mean_genomic_window_size = info['statistics']['window_size_mean']/chr_length(info['chr_id']) 
        ymax = max(YMAX[columns*(g//columns):columns*(g//columns+1)])
        ax1.errorbar( 0.88-mean_genomic_window_size, -0.76*ymax,marker=None, ls='none', xerr=25*mean_genomic_window_size, linewidth=2*scale, color='k', capsize=4*scale, zorder=20)
        ax1.text(     0.88-mean_genomic_window_size, -0.82*ymax, '25 GW',  horizontalalignment='center', verticalalignment='top',fontsize=2*fs//3, zorder=20)
        ax1.plot([0,1],[0,0],color='black', ls='dotted',alpha=0.7,zorder=0, linewidth=2*scale, scalex=False, scaley=False)
        ax1.set_title(identifier,fontsize=fs)
        export[identifier]['size_of_25_GW'] = 25*mean_genomic_window_size

    
    if len({info['chr_id'] for likelihoods,info in DATA.values()})==1:
        plot_monosomy_crossovers = True
        monosomy_crossovers = clusters_crossovers_v2(crossovers, dx=1000000/l) #1/num_of_bins[info['chr_id']])
        export['monosomy'] = {'crossovers': [*monosomy_crossovers]}
    else:
        plot_monosomy_crossovers = False
        
    
    for g,ax1 in enumerate(AX[:len(DATA)]):
        ymax = max(YMAX[columns*(g//columns):columns*(g//columns+1)])
        ax1.set_ylim((-1.01*ymax,+1.01*ymax))
        ax1.set_xlim((0,1)) 
        
        #Replace ticks along the x-axis 
        X_ticks = [i/10 for i in range(0,11,2)]
        X_labels = [('%g' % j) for j in X_ticks] 
        ax1.set_xticks(X_ticks)
        ax1.set_xticklabels(X_labels)
        
        ax1.tick_params(axis='x', labelsize=fs) 
        ax1.tick_params(axis='y', labelsize=fs)
        ax1.xaxis.set_tick_params(width=2*scale)
        ax1.yaxis.set_tick_params(width=2*scale)
        ###ax1.grid(color='black', linestyle='-.', linewidth=1,alpha=0.5)
        for axis in ['top','bottom','left','right']:
            ax1.spines[axis].set_linewidth(2*scale)
            
        for i in crossovers[g]:
            ax1.plot([i,i],[-1.01*ymax,1.01*ymax],color='purple', ls='dotted',alpha=0.7,zorder=19, linewidth=2*scale, scalex=False, scaley=False)
        
        if plot_monosomy_crossovers:
            for i in monosomy_crossovers:
                ax1.plot([i,i],[-1.01*ymax,1.01*ymax],color='red', ls='solid',alpha=0.7,zorder=19, linewidth=2*scale, scalex=False, scaley=False)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    
    plt.xlabel('Chromosomal position', fontsize=fs,labelpad=23*scale)
    plt.ylabel('Log-likelihood ratio', fontsize=fs,labelpad=45*scale)        
    if kwargs.get('title',None): 
        fig.suptitle(kwargs['title'], fontsize=int(1.2 * fs), color='black', fontweight="bold")    

    for l in range(1,len(AX)-len(DATA)+1):
        AX[-l].tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False, width=0)
        for axis in ['top','bottom','left','right']:
            AX[-l].spines[axis].set_visible(False)
    
    for l in range(len(DATA)-columns,len(DATA)):
        AX[l].xaxis.set_tick_params(labelbottom=True)
          
        
    if save!='':
        print('Saving plot...')
        #ax1.set_title(save.rpartition('/')[-1].removesuffix('.png'))
        #plt.tight_layout()
        #extension = 'svg'
        plt.savefig('.'.join([save,extension]), format=extension) #bbox_inches='tight'
        plt.close(fig)
    else:
       #plt.tight_layout()
       plt.show() 
       
    return export

def single_plot(likelihoods,info,**kwargs):
    """ Creates a figure  depicts the log-likelihood ratio vs. chromosomal
        position for (a) BPH over disomy, (b) disomy over SPH and (c) SPH over 
        monosomy. """
        
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    #from scipy.interpolate import interp1d
    
    scale = kwargs.get('scale', 1)
    z_score = kwargs.get('z_score', 1.96)
    bin_size = kwargs.get('bin_size', 4000000)
    save = kwargs.get('save', '')
    extension = kwargs.get('extension', 'svg')
    
    if save!='':
        mpl.use('Agg')
    else:
        #['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 'WXCairo', 'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']
        mpl.use('Qt5Agg')
    
    
    num_of_bins = {'chr'+str(i): chr_length('chr'+str(i))//bin_size for i in [*range(1,23)]+['X','Y']}

    fs = 24 * scale
    
    LLRs = {window: tuple(starmap(LLR,likelihoods_in_window))
                       for window,likelihoods_in_window in likelihoods.items()} 
                                                        
        
    fig,(ax1)=plt.subplots(1,1, figsize=(16 * scale, 9 * scale))
    fig.subplots_adjust(left=0, bottom=0, right=1, top=0.95, wspace=None, hspace=None)
    H = {}
    X,Y,E = binning(LLRs, info, num_of_bins[info['chr_id']])
    Y = [(y if y else 0) for y in Y]
    E = [(z_score*e if e else 0) for e in E]        
    T = [(x[1]+x[0])/2 for x in X]            
    
    
    steps_x = [X[0][0]]+[i[1] for i in X[:-1] for j in (1,2)]+[X[-1][1]]
    steps_y = [i for i in Y for j in (1,2)]
    H[0] = ax1.plot(steps_x, steps_y, label='BPH vs. SPH',color=(104/255,162/255,104/255), linewidth=2*scale, zorder=10, scalex=True, scaley=True, alpha=0.8)
    P = [(x[1]-x[0])/2 for x in X]                
    ax1.errorbar(T, Y, xerr = P, color=(104/255,162/255,104/255),marker=None, ls='none',alpha=1, zorder=13, linewidth=3*scale) 
    
    ax1.errorbar(T, Y, yerr = E, ecolor='black',marker=None, ls='none',alpha=0.2, linewidth=scale, zorder=15) 


    ax1.tick_params(axis='x', labelsize=fs) 
    ax1.tick_params(axis='y', labelsize=fs)
    ax1.xaxis.set_tick_params(width=scale)
    ax1.yaxis.set_tick_params(width=scale)
    ###ax1.grid(color='black', linestyle='-.', linewidth=1,alpha=0.5)
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(scale)
        
    ax1.set_title(kwargs.get('title',None), fontsize=int(fs), color='black', fontweight="bold")    
    ax1.set_ylabel('Log-likelihood ratio (normalized)', fontsize=fs,labelpad=2*scale)   
    ax1.set_xlabel('Chromosomal position (normalized)', fontsize=fs,labelpad=2*scale)

    
    #Replace ticks along the x-axis 
    X_ticks = [i/10 for i in range(0,11,2)]
    X_labels = [('%g' % j) for j in X_ticks] 
    ax1.set_xticks(X_ticks)
    ax1.set_xticklabels(X_labels)
    
    mean_genomic_window_size = info['statistics']['window_size_mean']/chr_length(info['chr_id']) 
    ymin,ymax = ax1.get_ylim()
    ax1.errorbar( 0.9-mean_genomic_window_size, ymin + 0.08*(ymax-ymin),marker=None, ls='none', xerr=25*mean_genomic_window_size, linewidth=2*scale, color='k', capsize=4*scale)
    ax1.text(   0.9-mean_genomic_window_size, ymin + 0.05*(ymax-ymin), '25 GW',  horizontalalignment='center', verticalalignment='top',fontsize=2*fs//3)
    ax1.plot([0,1],[0,0],color='black', ls='dotted',alpha=0.5)
    ax1.set_ylim((ymin,ymax))
    ax1.set_xlim((0,1))    
    
    if save!='':
        print('Saving plot...')
        #ax1.set_title(save.rpartition('/')[-1].removesuffix('.png'))
        #extension = 'svg'
        #plt.tight_layout()
        plt.savefig('.'.join([save,extension]), format=extension, bbox_inches='tight')
        plt.close(fig)
        
    else:
       plt.tight_layout()
       plt.show()
       
def wrap_panel_plot_many_cases(filenames, **kwargs):
    """ Wraps the function panel_plot to show a panel with many cases. """
    
    DATA = {}
    for llr_filename in filenames:
        likelihoods,info = load_likelihoods(llr_filename)
        if llr_filename[-6:]=='.LLR.p':
            identifer = llr_filename[:-6].rsplit('/',1).pop()
        elif llr_filename[-10:]=='.LLR.p.bz2':
            identifer = llr_filename[:-10].rsplit('/',1).pop()
        elif llr_filename[-9:]=='.LLR.p.gz':
            identifer = llr_filename[:9].rsplit('/',1).pop()
        else:
            identifer = llr_filename.rsplit('/',1).pop()
        DATA[identifer.split('.')[1]]=(likelihoods,info)
        show_info(info)
    title = f"{info['chr_id'].replace('chr', 'Chromosome '):s}; Contrasted with {identifer.split('.')[0]:s}"
    export = panel_plot(DATA, title=title, **kwargs)
    return export


def wrap_single_plot(llr_filename, **kwargs):
    """ Wraps the function single_plot. """
    likelihoods,info =  load_likelihoods(llr_filename)
    show_info(info)
    monosomy_id = info['disomy_obs_filename'].rsplit('/',1)[-1].split('.',1)[0]
    disomy_id = info['monosomy_obs_filename'].rsplit('/',1)[-1].split('.',1)[0]
    title = f"{disomy_id:s} contrasted w/ {monosomy_id:s}, {info['chr_id'].replace('chr', 'Chromosome '):s}"
    single_plot(likelihoods, info, title=title, **kwargs)
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Plots log-likelihood ratios (LLR) vs. chromosomal position from a LLR file.')
    parser.add_argument('llr_filename', metavar='LLR_FILENAME', type=str, nargs='+',
                        help='One or more LLR files created by CONTRAST_CROSSOVERS, containing likelihoods to observe reads under various aneuploidy landscapes .')
    parser.add_argument('-b', '--bin-size', type=int, metavar='INT', default=2000000,
                        help='The bin size in which the chromosome is divided. The default value is 2,000,000 bp.')
    parser.add_argument('-z', '--z-score', type=int, metavar='INT', default=1.96,
                        help='The z-score value for the confidence intervals. The default value is 1.96, which corresponds to confidence level of 95\%.')

    kwargs = vars(parser.parse_args())
    kwargs['pairs'] = [j.split(',') for j in kwargs.get('pairs','')]
    
    if  len(kwargs['llr_filename'])==1:
        kwargs['llr_filename'] = kwargs['llr_filename'].pop()
        wrap_single_plot(**kwargs)
    else:
        kwargs['filenames'] = kwargs['llr_filename']
        del kwargs['llr_filename']
        wrap_panel_plot_many_cases(**kwargs)
    sys.exit(0)

else:
    print('The module PLOT_PANEL was imported.')


####################################################
# Produce panel plots for all cases in the folders #
####################################################

#import os
#work_dir = 'results2/'
#identifiers = {i.split('.')[0] for i in os.listdir(work_dir) if i[-3:]=='bz2'}
#for identifier in identifiers:
#    try:
#        if not os.path.isfile(work_dir+identifier+'.svg'):
#            wrap_panel_plot(identifier,pairs=(('BPH','SPH'),),save=identifier,work_dir=work_dir, num_of_bins_in_chr21=20)
#    except Exception as e:
#        print(identifier,e)
filenames = [
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_monosomies/CHA-S-1-30-May-2020.CHA-S-9-30-May-2020.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_monosomies/CHA-S-1-30-May-2020.CHA-S-19-30-May-2020.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_monosomies/CHA-S-1-30-May-2020.CHA-S-2-30-May-2020.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_monosomies/CHA-S-1-30-May-2020.CHA-S-3-30-May-2020.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_monosomies/CHA-S-1-30-May-2020.CHA-S-18-30-May-2020.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_monosomies/CHA-S-1-30-May-2020.CHA-S-17-30-May-2020.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_monosomies/CHA-S-1-30-May-2020.CHA-S-16-30-May-2020.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_monosomies/CHA-S-1-30-May-2020.CHA-S-13-30-May-2020.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_monosomies/CHA-S-1-30-May-2020.CHA-S-7-30-May-2020.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_monosomies/CHA-S-1-30-May-2020.CHA-S-14-30-May-2020.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_monosomies/CHA-S-1-30-May-2020.CHA-S-11-30-May-2020.chr10.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_monosomies/CHA-S-1-30-May-2020.CHA-S-8-30-May-2020.chr10.LLR.p.bz2"]

#filenames = [
#"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe/TSA-M-1-21-Jul-2020.TSA-M-1-21-Jul-2020.chr15.LLR.p.bz2",
#"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe/TSA-M-1-21-Jul-2020.TSA-M-2-21-Jul-2020.chr15.LLR.p.bz2"]

filenames = [
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-2-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-3-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-4-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-6-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-7-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-8-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-9-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-10-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-11-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-14-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-15-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-16-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-17-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-18-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-19-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-20-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-21-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-22-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-23-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-25-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-28-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-29-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-30-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-31-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-35-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-36-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-40-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-42-29-Jul-2021.chr15.LLR.p.bz2",
"/home/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CC_analysis/results/CReATe_haploids/BEN-A-A-29-29-Jul-2021.BEN-A-A-44-29-Jul-2021.chr15.LLR.p.bz2"
]
#wrap_single_plot(llr_filename=filenames[0])

wrap_panel_plot_many_cases(filenames)