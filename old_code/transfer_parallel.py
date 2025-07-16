import pandas as pd
import numpy as np

import pandas as pd
import numpy as np

import time

import sys
from pyspi.calculator import Calculator

results_dir = 'results_GG'
dirdatain = '../data_GG/GAM'
fit_type = 'fit_301_0.01'

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


def load_pickle_with_pandas(filepath, retry_delay=1):
    while True:
        try:
            df = pd.read_pickle(filepath)
            print("Pickle file loaded successfully.")
            return df
        except Exception as e:
            print(f"Error reading pickle file: {e}. Retrying in {retry_delay} second(s)...")
            time.sleep(retry_delay)
            
def worker_process(index, start=0, end=140):
    species_df = load_pickle_with_pandas(f'{results_dir}/species_df.pkl')
    species_index = load_pickle_with_pandas(f'{results_dir}/species_index.pkl')
    myspecies = list(species_index.keys())

    bootstrapped_array = np.zeros((len(myspecies), end-start))
    for i, spec in enumerate(myspecies):
        bootstrapped_array[i] = stationary_bootstrap(species_df[spec]['y'][start:end])
    
    # calculate the table for the bootstrapped data
    calc = Calculator(dataset = bootstrapped_array, configfile='TE_configs/kraskov_k1.yaml')
    calc.compute()
    results = calc.table
    pd.to_pickle(results, f'{results_dir}/transfer_entropies/TE_matrix{index}_{start}:{end}.pkl')

def main():

    # process ID
    proces_ID = sys.argv[1]
    start = int(sys.argv[2])
    end = int(sys.argv[3])

    
    species_df = load_pickle_with_pandas(f'{dirdatain}/species_%s.pkl' %(fit_type))
    myspecies = list(species_df.columns)

    # create np.array of all the species y values
    # example of one of them: species_df['Betula']['y'] is one row
    array = np.zeros((len(myspecies), end-start))
    species_index = {}
    for i, spec in enumerate(myspecies): # time window from start to end
        array[i] = species_df[spec]['y'][start:end]
        species_index[spec] = i
    
    # save index dictionary to a file
    pd.to_pickle(species_index, f'{results_dir}/species_index.pkl')
    # save the species dataframe to a file
    pd.to_pickle(species_df, f'{results_dir}/species_df.pkl')

    #____________________________________________________________________
    # CALCULATE THE ORIGINAL TABLE
    
    if False:
        # We're using this: te_kraskov_NN-4_DCE_k-1_kt-1_l-1_lt-1
        # kraskov with 4 neighbors, DCE, k=1, kt=1, l=1, lt=1
        calc = Calculator(dataset = array, configfile='TE_configs/kraskov_k1.yaml')
        calc.compute()
        original_table = calc.table
        # save the original table to a file
        pd.to_pickle(original_table, f'{results_dir}/transfer_entropies/original_table_{start}:{end}.pkl')
        exit(0)

    #____________________________________________________________________
    # BOOTSTRAPPING

    start_time = time.perf_counter()

    worker_process(proces_ID, start, end)

    end_time = time.perf_counter()
    print(f"Time taken for bootstrapping: {end_time - start_time} seconds")
    #____________________________________________________________________
    # calculate the p-values

    # p_values = np.zeros(original_table.shape)
    # for i, bootstrapped_table in enumerate(bootstrapped_tables):
    #     p_values += (original_table < bootstrapped_table).astype(int)
    # p_values /= n_bootstrap

    # np.save(f'{results_dir}/p_values_TE_bootstrap_both_{n_bootstrap}iter_2.npy', p_values)

if __name__ == "__main__":
    main()