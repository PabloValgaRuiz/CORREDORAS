import pandas as pd
import numpy as np

import matplotlib as mpl
import pandas as pd
import numpy as np

import concurrent.futures
import time

from pyspi.calculator import Calculator

def stationary_bootstrap(x):
    # repeat the data once for periodicity, and not choose data out of the range
    x_period = np.concatenate((x, x))

    x_final = np.array([])
    while True:
        L = np.random.randint(low=1, high=len(x))
        x_0 = np.random.randint(low=0, high=len(x))
        
        if len(x_final) + L < len(x):
            # choose L data points from x_period starting from x_0
            x_final = np.concatenate((x_final, x_period[x_0:x_0+L]))
        else:
            # fill the remaining data points
            x_final = np.concatenate((x_final, x_period[x_0:x_0+len(x)-len(x_final)]))
            return x_final

# if main

def main():

    dirdatain = '../data/'

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


    fit_type = 'fit_140_0.01'
    species_df = pd.read_pickle('GAM_species/species_%s.pkl' %(fit_type))

    # create np.array of all the species y values
    # example of one of them: species_df['Betula']['y'] is one row
    array = np.zeros((len(myspecies), 140))
    species_index = {}
    for i, spec in enumerate(myspecies):
        array[i] = species_df[spec]['y']
        species_index[spec] = i
    
    # save index dictionary to a file
    pd.to_pickle(species_index, 'results2/species_index.pkl')

    p_values = np.zeros((len(myspecies), len(myspecies)))


    #____________________________________________________________________
    # CALCULATE THE ORIGINAL TABLE

    # We're using this: te_kraskov_NN-4_DCE_k-1_kt-1_l-1_lt-1
    # kraskov with 4 neighbors, DCE, k=1, kt=1, l=1, lt=1
    calc = Calculator(dataset = array, configfile='TE_configs/kraskov_k1.yaml')
    calc.compute()
    original_table = calc.table

    #____________________________________________________________________
    # BOOTSTRAPPING

    # generate 10 tables of random data the same shape
    n_bootstrap = 1000
    bootstrapped_tables = pd.Series(dtype=object, index=range(n_bootstrap))
    
    for i in range(n_bootstrap):

        bootstrapped_array = np.zeros((len(myspecies), 140))
        for i, spec in enumerate(myspecies):
            bootstrapped_array[i] = stationary_bootstrap(species_df[spec]['y'])
        
        # calculate the table for the bootstrapped data
        calc = Calculator(dataset = bootstrapped_array, configfile='TE_configs/kraskov_k1.yaml')
        calc.compute()

        bootstrapped_tables.append(calc.table)

    #____________________________________________________________________
    # calculate the p-values

    p_values = np.zeros(original_table.shape)
    for i, bootstrapped_table in enumerate(bootstrapped_tables):
        p_values += (original_table < bootstrapped_table).astype(int)
    p_values /= n_bootstrap

    np.save(f'results2/p_values_TE_bootstrap_both_{n_bootstrap}iter_2.npy', p_values)