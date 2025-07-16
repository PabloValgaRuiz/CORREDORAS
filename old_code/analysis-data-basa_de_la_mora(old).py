import pyEDM as edm
import concurrent.futures
import time
import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pygam as gm
import copy as cp

import pyEDM as edm
import time
import math
from matplotlib.ticker import (MultipleLocator)

from PyIF import te_compute as te


def concurrent_TE(columns, target, fit_type, k, embedding):
    #___________________________________________ Load the data ___________________________________________
    species = pd.read_pickle('GAM_species/species_%s.pkl' %(fit_type))

    pair = pd.DataFrame({columns : species[columns]['y'], target : species[target]['y']}, index=species[columns]['x'])

    TE1 = te.te_compute(pair[columns].values, pair[target].values,
                           k=k, embedding=embedding, safetyCheck=False, GPU=False)
    TE2 = te.te_compute(pair[target].values, pair[columns].values,
                            k=k, embedding=embedding, safetyCheck=False, GPU=False)
    
    results = {}
    results[f'{columns}:{target}'] = TE2
    results[f'{target}:{columns}'] = TE1
    
    print(f'finished {columns}:{target}. Results: TE1={TE1}, TE2={TE2}')

    return results

def concurrent_CCM(columns, target, fit_type):
    #___________________________________________ Load the data ___________________________________________
    # pair = pd.read_pickle('pair_files/pair_%s_%s.pkl' %(columns, target))
    #______________________________ Find the optimal embedding dimension ______________________________
    
    species = pd.read_pickle('GAM_species/species_%s.pkl' %(fit_type))

    pair = pd.DataFrame({columns : species[columns]['y'], target : species[target]['y']}, index=species[columns]['x'])

    dfE_columns = edm.EmbedDimension(dataFrame=pair, lib = [1, len(pair)], pred = [1, len(pair)],
                            columns=columns, target=columns, showPlot=False, maxE=10)
    dfE_target = edm.EmbedDimension(dataFrame=pair, lib = [1, len(pair)], pred = [1, len(pair)],
                            columns=target, target=target, showPlot=False, maxE=10)
    dfE = pd.concat([dfE_columns, dfE_target])
    E = int(dfE.iloc[dfE['rho'].idxmax()]['E'])


    #________________________________________ Carry out CCM ___________________________________________
    
    result = edm.CCM(dataFrame=pair,
                    columns = columns, target = target,
                    E = E, Tp=0, libSizes =[E+2, len(pair)],
                    sample = 20, showPlot=False)

    return result
        
if __name__ == '__main__':

    mpl.rcParams['figure.figsize'] = (5,5)
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['mathtext.fontset'] = 'cm'
    mpl.rcParams['mathtext.rm'] = 'serif'
    mpl.rcParams["figure.autolayout"] = True

    print('libraries loaded!')


    nr_times = 154

    dirdatain = '../data/'

    species = {} # dictionary with the IDs of the species (as they change across networks)

    print('Common settings done!')


    mydf = pd.read_csv(dirdatain+'bsm_raw_pollen.csv')
    mydf = mydf.fillna(0)
    mydf['age'] = -mydf['age']
    mydf = mydf.sort_values(by=['age'])

    # rename the columns: remove the spaces and replace them with underscores

    mydf.columns = [col.replace(' ', '_') for col in mydf.columns]

    all_my_species = list(mydf.columns)
    myspecies = ['Betula','Artemisia']
    myspecies = ['Betula', 'Pinus', 'Corylus', 'Artemisia', 'Olea']
    myspecies = all_my_species[1:] # ALL SPECIES

    for spec in myspecies[:]: # ITERATE OVER A COPY
        if mydf[spec].sum() == 0:
            print('Species %s has all zeros!' %(spec))
            myspecies.remove(spec)
            mydf = mydf.drop(columns=[spec])
    
#____________________________________________________ GAM ____________________________________________________
    k = 2
    embedding = 2

    fit_type = 'fit_140_0.01'
    if 0: # DO NOT DO GAM IF IT ALREADY EXISTS
        print('\nComputing GAM ....')
        # for transfer entropy

        species = {}
        vx = gm.utils.make_2d(mydf['age'], verbose=False).astype('float')
        for spec in myspecies:
            print('\tIterating on species %s ...' %(spec))
            vy = mydf[spec]

            mygam = gm.GAM(terms='auto', n_splines=140, lam=0.01).fit(vx, vy)
            # mygam = gm.GAM(terms='auto').gridsearch(vx, vy)

            species[spec] = {}
            species[spec]['prev_x'] = mydf['age'].to_numpy()
            species[spec]['prev_y'] = mydf[spec].to_numpy()

            XX = mygam.generate_X_grid(term=0, n=140)
            YY = mygam.predict(XX)
            species[spec]['x'] = cp.copy(XX)[:,0]
            species[spec]['y'] = cp.copy(YY)

            # print(species[spec]['x'].shape)
            # print(species[spec]['y'].shape)
            # print(species[spec]['prev_x'].shape)
            # print(species[spec]['prev_y'].shape)

        print('Ho finito!\n')

        species_df = pd.DataFrame(species)
        species_df.to_pickle('GAM_species/species_%s.pkl' %(fit_type))
    
#____________________________________________________ PLOT GAM ________________________________________________
    # fig = plt.figure(figsize=(8,8))
    # for spec in myspecies:
    #     # scatter with x markers instead of points
    #     scatter = plt.scatter(species[spec]['prev_x'], species[spec]['prev_y'], label=spec, marker="x", s=20)
    #     # line same color as the scatter before
    #     color = scatter.get_facecolors()[0]
    #     plt.plot(species[spec]['x'], species[spec]['y'], label=spec, color=color)
    # ax = plt.gca()
    # ax.spines[['right', 'top']].set_visible(False)
    # plt.legend()
    # plt.savefig('plots-tests/basa_gam_%s.pdf' %(fit_type))
#______________________________________________________ CCM ____________________________________________________

    results_CCM = {}
    futures_CCM = {}

    results_te = {}
    futures_te = {}

    # time it
    start = time.time()
    with concurrent.futures.ProcessPoolExecutor(max_workers=24) as executor:
        for i,columns in enumerate(myspecies):
            for j,target in enumerate(myspecies):
                if i <= j:  # The crucial condition: i must be less than j
                    continue
                print('Computing %s:%s' %(columns, target))
                # futures_CCM[columns+':'+target] = executor.submit(concurrent_CCM, columns, target, fit_type)
                futures_te[columns+':'+target] = executor.submit(concurrent_TE, columns, target, fit_type, k, embedding)

                # futures[columns+':'+target] = concurrent_CCM(columns, target, fit_type)
    
    for key in futures_te.keys():
        results_te[key] = futures_te[key].result()

    # for key in futures_CCM:
    #     results_CCM[key] = futures_CCM[key].result()
        
    finish = time.time()

    # store the results in a pickle file
    # pd.to_pickle(results_CCM, f'results/results_basa_CCM_{fit_type}.pkl')
    pd.to_pickle(results_te, f'results/results_basa_TE_{fit_type}_k_{k}_emb_{embedding}.pkl')
    # print in seconds
    print(f'Finished in {finish-start} seconds')
