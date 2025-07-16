import pandas as pd
import numpy as np

import pandas as pd
import numpy as np
from statsmodels.tsa.stattools import grangercausalitytests
import time

import sys

results_dir = 'results_GG'
dirdatain = '../data_GG/GAM'
fit_type = 'fit_301_0.01'

# results_dir = 'results_basa'
# dirdatain = '../code/GAM_species'
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

def grangers_causation_matrix(data, variables, test='ssr_chi2test'):    

    df_causality = pd.DataFrame(np.zeros((len(variables), len(variables))), columns=variables, index=variables)
    df_p_values = pd.DataFrame(np.zeros((len(variables), len(variables))), columns=variables, index=variables)
    for i,c in enumerate(df_causality.columns):
        for j,r in enumerate(df_causality.index): # data[[y,x]] BECAUSE THIS FUNCTION RETURNS THE OPPOSITE, Y CAUSING X --> NOW, IT RETURNS X CAUSING Y
            data_pair = pd.DataFrame({'x': data[i], 'y': data[j]})
            
            try:
                test_result = grangercausalitytests(data_pair, maxlag=1, verbose=False)
                maxlag = 1
                p_values = [round(test_result[i+1][0][test][1],4) for i in range(maxlag)]
                min_p_value = np.min(p_values)
                df_p_values.loc[r, c] = min_p_value

                causality = test_result[1][1][1].params[1]
                df_causality.loc[r, c] = causality
            except Exception as e:
                # probably the time series are almost all zeros
                df_p_values.loc[r, c] = 1
                df_causality.loc[r, c] = 0
            # print(c)
            # print(r)
            # print(min_p_value)
            # print(df.loc[r, c])

    return df_causality.values, df_p_values.values

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
        array[i][array[i] < 0.001] = 0 # set to 0 if the value is less than 0.001 cause GAM is tricky
        species_index[spec] = i

    # save index dictionary to a file
    pd.to_pickle(species_index, f'{results_dir}/species_index.pkl')
    # save the species dataframe to a file
    pd.to_pickle(species_df, f'{results_dir}/species_df.pkl')

    #____________________________________________________________________
    # CALCULATE THE ORIGINAL TABLE
    start_time = time.perf_counter()
    
    original_table_causality, original_table_p_values = grangers_causation_matrix(array, variables=myspecies)
    # save the original table to a file
    pd.to_pickle(original_table_causality, f'{results_dir}/granger_causalities/original_table_{start}-{end}.pkl')
    pd.to_pickle(original_table_p_values, f'{results_dir}/granger_causalities/original_table_p_values_{start}-{end}.pkl')
    
    end_time = time.perf_counter()
    print(f"Time taken for bootstrapping: {end_time - start_time} seconds")

    exit(0)

    #____________________________________________________________________
    # BOOTSTRAPPING



    #____________________________________________________________________
    # calculate the p-values

    # p_values = np.zeros(original_table.shape)
    # for i, bootstrapped_table in enumerate(bootstrapped_tables):
    #     p_values += (original_table < bootstrapped_table).astype(int)
    # p_values /= n_bootstrap

    # np.save(f'{results_dir}/p_values_TE_bootstrap_both_{n_bootstrap}iter_2.npy', p_values)

if __name__ == "__main__":
    main()