import pandas as pd
import numpy as np
from statsmodels.tsa.api import VAR

import time

import sys
import os

results_dir = 'results_GG'
dirdatain = 'GAM_species_GG'
fit_type = 'fit_301_0.01'

# results_dir = 'results_basa'
# dirdatain = 'GAM_species'
# fit_type = 'fit_140_0.01'

def load_pickle_with_pandas(filepath, retry_delay=1):
    while True:
        try:
            df = pd.read_pickle(filepath)
            print("Pickle file loaded successfully.")
            return df
        except Exception as e:
            print(f"Error reading pickle file: {e}. Retrying in {retry_delay} second(s)...")
            time.sleep(retry_delay)

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

def grangers_causation_matrix(data : np.array, variables : list, conditional_variables : dict):

    df_causality = pd.DataFrame(np.zeros((len(variables), len(variables))), columns=variables, index=variables)
    df_p_values = pd.DataFrame(np.zeros((len(variables), len(variables))), columns=variables, index=variables)

    abundances_sum = np.sum(data, axis=0) # sum of all species in every time step

    for i,r in enumerate(df_causality.index):
        for j,c in enumerate(df_causality.columns): # data[[x,y]] opposite to the other granger causality test file, which returns y causing x
                                                                    # the rest of the species, not X or Y
            data_tuple = pd.DataFrame({'X': data[i], 'Y': data[j]})
            data_tuple = data_tuple.assign(**conditional_variables)  # add conditional variables to the DataFrame

            try:
                model = VAR(data_tuple)
                results = model.fit(maxlags=1)
                causality = results.test_causality(caused='Y', causing=['X'], kind='f')

                df_p_values.loc[r, c] = causality.pvalue
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

    # load delta13C
    try:
        conditional_variables['delta13C'] = load_pickle_with_pandas(f'{dirdatain}/d13gam_{fit_type}.pkl')['d13C (permil)'].to_numpy()
        conditional_variables['delta13C'] = conditional_variables['delta13C'][start:end] # select the time window
    except FileNotFoundError as e:
        print(f"File not found: {dirdatain}/d13gam_{fit_type}.pkl\n check we're not in basa de la mora") # probably using garba guracha (GG) data and not basa

    # load multipliers
    try:
        conditional_variables['multipliers'] = load_pickle_with_pandas(f'{dirdatain}/gam_multipliers_{fit_type}.pkl').values[start:end]
    except FileNotFoundError as e:
        print(f"File not found: {results_dir}/multipliers.pkl\n check we're not in basa de la mora")

    return conditional_variables

def subset_dict(d, keys):
    return {k: d[k] for k in keys if k in d}


def bootstrap_worker(bs_idx, array: np.array, myspecies: list, conditional_variables: dict, boot_dir:str, start:int, end:int):
        # set a seed for numpy random
        np.random.seed(bs_idx)
        bootstrapped_array = np.zeros((len(myspecies), end-start))
        for i, spec in enumerate(myspecies):
            bootstrapped_array[i] = stationary_bootstrap(array[i])
        table_causality_bs, table_p_values_bs = grangers_causation_matrix(
                    bootstrapped_array,
                    variables=myspecies,
                    conditional_variables=conditional_variables
                    )
        pd.to_pickle(table_causality_bs, f'{boot_dir}/table_{start}-{end}_bs{bs_idx}.pkl')
        pd.to_pickle(table_p_values_bs,  f'{boot_dir}/table_p_values_{start}-{end}_bs{bs_idx}.pkl')
        return bs_idx

def main():

    # process ID
    proces_ID = sys.argv[1]
    start = int(sys.argv[2])
    end = int(sys.argv[3])

    out_base = f'{results_dir}/bootstrap_conditional_GC/abundances_cond2'
    
    species_df = load_pickle_with_pandas(f'{dirdatain}/species_%s.pkl' %(fit_type))
    myspecies = list(species_df.columns)

    # create np.array of all the species y values
    # example of one of them: species_df['Betula']['y'] is one row
    array = np.zeros((len(myspecies), end-start))
    species_index = {}
    for i, spec in enumerate(myspecies): # time window from start to end
        array[i] = species_df[spec]['y'][start:end]
        array[i][array[i] < 0.001] = 0 # set to 0 if the value is less than 0.001 (in abundances) cause GAM is tricky
        species_index[spec] = i        #                                  (0.01 in PAR without lycopodium)

    # save index dictionary to a file
    pd.to_pickle(species_index, f'{results_dir}/species_index.pkl')
    # save the species dataframe to a file
    pd.to_pickle(species_df, f'{results_dir}/species_df.pkl')

    # load the conditional variables
    conditional_variables = load_conditional_variables(start, end)
    conditional_variables = subset_dict(conditional_variables, ['delta13C', 'multipliers'])
    #____________________________________________________________________
    # CALCULATE THE ORIGINAL TABLE
    start_time = time.perf_counter()
    
    table_causality, table_p_values = grangers_causation_matrix(
                array,
                variables=myspecies,
                conditional_variables=conditional_variables
                )
    
    # save the table to a file
    pd.to_pickle(table_causality, f'{out_base}/original_table_{start}-{end}.pkl')
    pd.to_pickle(table_p_values,  f'{out_base}/original_table_p_values_{start}-{end}.pkl')
    
    end_time = time.perf_counter()
    print(f"Time taken: {end_time - start_time} seconds")

    # bootstrapping

    from concurrent.futures import ProcessPoolExecutor
    import functools

    boot_dir = out_base + '/bootstrap_tables'
    os.makedirs(boot_dir, exist_ok=True)

    
    
    n_bootstrap = 1000
    
    start_time = time.perf_counter()
    with ProcessPoolExecutor(max_workers=20) as executor:
        func = functools.partial(
            bootstrap_worker,
            array=array,
            conditional_variables=conditional_variables,
            myspecies=myspecies,
            boot_dir=boot_dir,
            start=start,
            end=end
        )
        for b in executor.map(func, range(n_bootstrap)):
            print(f"Completed bootstrap {b}")

    end_time = time.perf_counter()
    print(f"Bootstrapping completed. Time taken: {end_time - start_time} seconds")

if __name__ == "__main__":
    main()