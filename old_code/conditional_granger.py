import pandas as pd
import numpy as np

import pandas as pd
import numpy as np
from statsmodels.tsa.api import VAR

import time

import sys

# results_dir = 'results_GG'
# dirdatain = '../data_GG/GAM'
# fit_type = 'fit_301_0.01'

results_dir = 'results_basa'
dirdatain = '../code/GAM_species'
fit_type = 'fit_140_0.01'

def load_pickle_with_pandas(filepath, retry_delay=1):
    while True:
        try:
            df = pd.read_pickle(filepath)
            print("Pickle file loaded successfully.")
            return df
        except Exception as e:
            print(f"Error reading pickle file: {e}. Retrying in {retry_delay} second(s)...")
            time.sleep(retry_delay)

def grangers_causation_matrix(data, variables):    

    df_causality = pd.DataFrame(np.zeros((len(variables), len(variables))), columns=variables, index=variables)
    df_p_values = pd.DataFrame(np.zeros((len(variables), len(variables))), columns=variables, index=variables)

    abundances_sum = np.sum(data, axis=0) # sum of all species in every time step
    
    for i,r in enumerate(df_causality.index):
        for j,c in enumerate(df_causality.columns): # data[[x,y]] opposite to the other granger causality test file, which returns y causing x
                                                                    # the rest of the species, not X or Y
            data_tuple = pd.DataFrame({'X': data[i], 'Y': data[j], 'C': abundances_sum - (data[i] + data[j])})

            try:
                model = VAR(data_tuple[['X', 'Y', 'C']])
                results = model.fit(maxlags=1)
                causality = results.test_causality(caused='Y', causing=['X'], kind='f')

                df_p_values.loc[r, c] = causality.pvalue
                df_causality.loc[r, c] = results.coefs[0][1][0] # get the coefficient of X in the regression of Y
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
    pd.to_pickle(original_table_causality, f'{results_dir}/conditional_GC/original_table_{start}:{end}.pkl')
    pd.to_pickle(original_table_p_values, f'{results_dir}/conditional_GC/original_table_p_values_{start}:{end}.pkl')
    
    end_time = time.perf_counter()
    print(f"Time taken for bootstrapping: {end_time - start_time} seconds")


if __name__ == "__main__":
    main()