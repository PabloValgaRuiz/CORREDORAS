import pandas as pd
import numpy as np
from statsmodels.tsa.api import VAR
import statsmodels.api as sm

import time

import sys
import os

# database = 'moossee'
# results_dir = f'results/{database}'
# dirdatain = f'GAM_species/{database}'
# fit_type = 'fit_815_0.01'

# database = 'Burgäschisee'
# results_dir = f'results/{database}'
# dirdatain = f'GAM_species/{database}'
# fit_type = 'fit_723_0.01'

database = 'GG'
results_dir = f'results/{database}'
dirdatain = f'GAM_species/{database}'
fit_type = 'fit_300_0.01'

# BASA RECORDAR AÑADIR EL D13GAM
# database = 'basa'
# results_dir = f'results/{database}'
# dirdatain = f'GAM_species/{database}'
# # dirdatain = f'linear_species_basa/{database}'
# fit_type = 'fit_140_0.01'

NORMALIZE = False
PAR = False
LAG = 1

def load_pickle_with_pandas(filepath, retry_delay=1):
    while True:
        try:
            df = pd.read_pickle(filepath)
            print("Pickle file loaded successfully.")
            return df
        except Exception as e:
            print(f"Error reading pickle file: {e}. Retrying in {retry_delay} second(s)...")
            time.sleep(retry_delay)

# returns an array of indices to bootstrap the data



def stationary_bootstrap(SIZE):
    x = np.arange(SIZE, dtype=int)
    # repeat the data once for periodicity, and not choose data out of the range
    x_period = np.concatenate((x, x), dtype=int)
    
    x_final = np.array([], dtype=int)
    while True:
        L = int(np.random.exponential(scale=len(x), size=None))
        while L > len(x):
            L = int(np.random.exponential(scale=len(x), size=None))

        x_0 = np.random.randint(low=0, high=len(x))
        
        if len(x_final) + L < len(x):
            # choose L data points from x_period starting from x_0
            x_final = np.concatenate((x_final, x_period[x_0:x_0+L]))
        else:
            # fill the remaining data points
            x_final = np.concatenate((x_final, x_period[x_0:x_0+len(x)-len(x_final)]))
            return x_final

def grangers_causation_matrix_OLS(data : np.array, variables : list, conditional_variables : np.array):
    
    lag = LAG

    df_causality = pd.DataFrame(np.zeros((len(variables), len(variables))), columns=variables, index=variables)
    df_p_values = pd.DataFrame(np.zeros((len(variables), len(variables))), columns=variables, index=variables)

    for i,r in enumerate(df_causality.index):
        for j,c in enumerate(df_causality.columns): 
            
            # Manually shift the data to use a specific lag
            # Target: Y from 'lag' to the end
            # Predictors: X, Y, and Conditionals from beginning to end-lag
            Y_target = data[j][lag:]
            Y_past = data[j][:-lag]
            X_past = data[i][:-lag]
            Cond_past = conditional_variables[j][:, :-lag].T

            # Stack predictors: Intercept, Y_t-k, X_t-k, and Conditionals_t-k
            # column_stack creates a matrix where each variable is a column
            X_matrix = np.column_stack([np.ones(len(Y_target)), Y_past, X_past, Cond_past])

            try:
                # Use OLS to find the coefficient of X_past on Y_target
                model = sm.OLS(Y_target, X_matrix)
                results = model.fit()
                
                # index 0: intercept, 1: Y_past, 2: X_past
                df_causality.loc[r, c] = results.params[2]
                df_p_values.loc[r, c] = results.pvalues[2]
            except Exception as e:
                df_p_values.loc[r, c] = 1
                df_causality.loc[r, c] = 0

    return df_causality.values, df_p_values.values

def grangers_causation_matrix(data : np.array, variables : list, conditional_variables : np.array):

    df_causality = pd.DataFrame(np.zeros((len(variables), len(variables))), columns=variables, index=variables)
    df_p_values = pd.DataFrame(np.zeros((len(variables), len(variables))), columns=variables, index=variables)

    # abundances_sum = np.sum(data, axis=0) # sum of all species in every time step

    # check sum
    

    for i,r in enumerate(df_causality.index):
        for j,c in enumerate(df_causality.columns): # data[[x,y]] opposite to the other granger causality test file, which returns y causing x
                                                                    # the rest of the species, not X or Y
            data_tuple = pd.DataFrame({'X': data[i], 'Y': data[j]})

            # add conditional variables BOOTSTRAPPED LIKE THE CAUSED VARIABLE Y (j), so that they can still explain it
            data_tuple = pd.concat([data_tuple, pd.DataFrame(conditional_variables[j].T)], axis=1)
            # if (i == 0) & (j == 1):
            #     print(data_tuple)
            #     exit(0)
            try:
                model = VAR(data_tuple)
                results = model.fit(maxlags=1)
                # causality = results.test_causality(caused='Y', causing=['X'], kind='f')

                df_p_values.loc[r, c] = 1 # causality.pvalue
                df_causality.loc[r, c] = results.coefs[0][1][0] # get the coefficient of X in the regression of Y
            except Exception as e:
                # likely if the time series are all or almost all zeros
                df_p_values.loc[r, c] = 1
                df_causality.loc[r, c] = 0
            # print(c)
            # print(r)
            # print(min_p_value)
            # print(df.loc[r, c])

    return df_causality.values, df_p_values.values

def load_conditional_variables(start, end):
    conditional_variables = {}

    # # load delta13C
    if database == 'basa':
        try:
            conditional_variables['delta13C'] = load_pickle_with_pandas(f'{dirdatain}/d13gam_{fit_type}.pkl')['d13C (permil)'].to_numpy()
            conditional_variables['delta13C'] = conditional_variables['delta13C'][start:end] # select the time window
        except FileNotFoundError as e:
            print(f"File not found: {dirdatain}/d13gam_{fit_type}.pkl\n check we're not in basa de la mora") # probably using garba guracha (GG) data and not basa

    # load multipliers
    try:
        conditional_variables['multipliers'] = load_pickle_with_pandas(f'{dirdatain}/gam_multipliers_{fit_type}.pkl').values[start:end][:,1]
    except FileNotFoundError as e:
        print(f"File not found: {results_dir}/multipliers.pkl\n check we're not in basa de la mora")

    return conditional_variables

def subset_dict(d, keys):
    return {k: d[k] for k in keys if k in d}

import matplotlib.pyplot as plt
def bootstrap_worker(bs_idx, array: np.array, myspecies: list, conditional_variables: dict, boot_dir:str, start:int, end:int):
        # set a seed for numpy random
        np.random.seed(bs_idx)
        
        bootstrapped_array = np.zeros((len(myspecies), end-start))
        bootstrapped_cond_vars = np.zeros((len(myspecies), len(conditional_variables), end-start))

        for i, spec in enumerate(myspecies):
            # generate bootstrapped indices for any series
            bootstrapped_indices = stationary_bootstrap(end-start)
            # bootstrap both the species data and conditionals IN THE SAME EXACT WAY SO THAT THE CORRELATION REMAINS
            bootstrapped_array[i] = array[i][bootstrapped_indices]
            for j, key in enumerate(conditional_variables):
                bootstrapped_cond_vars[i][j] = conditional_variables[key][bootstrapped_indices]
        
        table_causality_bs, table_p_values_bs = grangers_causation_matrix_OLS(
                    bootstrapped_array,
                    variables=myspecies,
                    conditional_variables=bootstrapped_cond_vars
                    )
        pd.to_pickle(table_causality_bs, f'{boot_dir}/table_{start}-{end}_bs{bs_idx}.pkl')
        # pd.to_pickle(table_p_values_bs,  f'{boot_dir}/table_p_values_{start}-{end}_bs{bs_idx}.pkl')
        return bs_idx

def main():
    
    # process ID
    proces_ID = sys.argv[1]
    start = int(sys.argv[2])
    end = int(sys.argv[3])

    out_base = f'{results_dir}/bootstrap_conditional_GC_OLS/abundances_cond1_lag{LAG}'

    os.makedirs(out_base, exist_ok=True)

    species_df = load_pickle_with_pandas(f'{dirdatain}/species_%s.pkl' %(fit_type))
    species_index = load_pickle_with_pandas(f'{dirdatain}/species_index.pkl')
    myspecies = list(species_df.columns)

    # check the df
    print(species_df.head())

    # create np.array of all the species y values
    # example of one of them: species_df['Betula']['y'] is one row
    array = np.zeros((len(myspecies), end-start))

    for spec in myspecies: # time window from start to end
        i = species_index[spec]
        array[i] = species_df[spec]['y'][start:end] # this goes over the end at 140, but ignores it.

        # CAN DO THIS IN PLOT RESULTS FILE
        # array[i][array[i] < 0] = 0 # set to 0 if the value is less than 0.001 (in abundances) cause GAM is tricky
                                        #                                  (0.01 in PAR without lycopodium)        
    
    # load the conditional variables
    conditional_variables = load_conditional_variables(start, end)
    # calculate total biomass
    conditional_variables['biomass'] = conditional_variables['multipliers'] * array.sum(axis=0)

    if database in ['basa']: # BASA DE LA MORA
        conditional_variables = subset_dict(conditional_variables, ['delta13C'])

    elif database in ['GG', 'moossee', 'Burgäschisee']: # GARBA GURACHA, MOOSSEE y Burgäschisee
        conditional_variables = subset_dict(conditional_variables, ['biomass'])

    # NORMALIZE ALL TIME SERIES
    if(NORMALIZE == True):
        for i,spec in enumerate(myspecies):
            if np.std(array[i]) > 0:
                array[i] = (array[i] - np.mean(array[i])) / np.std(array[i])
        for key in conditional_variables:
            if np.std(conditional_variables[key]) > 0:
                conditional_variables[key] = (conditional_variables[key] - np.mean(conditional_variables[key])) / np.std(conditional_variables[key])

    # USE PAR INSTEAD OF COUNT
    if (PAR == True):
        for spec in myspecies:
            i = species_index[spec]
            array[i] = array[i] * conditional_variables['multipliers']

    #____________________________________________________________________
    # CALCULATE THE ORIGINAL TABLE
    start_time = time.perf_counter()
    
    # to make it compatible with bootstrapped grangers causation matrix, give every species a conditional variable of
    # its own (even if in the original series they are all the same, copies)
    conditional_variables_per_species = np.zeros((len(myspecies), len(conditional_variables), end-start))
    for i, spec in enumerate(myspecies):
        for j, key in enumerate(conditional_variables):
            conditional_variables_per_species[i][j] = conditional_variables[key]
    
    # DO CAUSALITY
    table_causality, table_p_values = grangers_causation_matrix_OLS(
        array,
        variables=myspecies,
        conditional_variables=conditional_variables_per_species
        )

    # save the table to a file
    pd.to_pickle(table_causality, f'{out_base}/original_table_{start}-{end}.pkl')
    # pd.to_pickle(table_p_values,  f'{out_base}/original_table_p_values_{start}-{end}.pkl')

    end_time = time.perf_counter()
    print(f"Time taken: {end_time - start_time} seconds")

    # bootstrapping

    from concurrent.futures import ProcessPoolExecutor
    import functools

    boot_dir = out_base + '/bootstrap_tables'
    os.makedirs(boot_dir, exist_ok=True)

    n_bootstrap = range(0, 4000)

    start_time = time.perf_counter()
    with ProcessPoolExecutor(max_workers=None) as executor:
        func = functools.partial(
            bootstrap_worker,
            array=array,
            conditional_variables=conditional_variables,
            myspecies=myspecies,
            boot_dir=boot_dir,
            start=start,
            end=end
        )
        for b in executor.map(func, n_bootstrap):
            print(f"Completed bootstrap {b}")

    end_time = time.perf_counter()
    print(f"Bootstrapping completed. Time taken: {end_time - start_time} seconds")

if __name__ == "__main__":
    main()