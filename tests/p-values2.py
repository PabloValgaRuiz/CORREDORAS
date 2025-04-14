import concurrent.futures
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyEDM as edm
from PyIF import te_compute as te
from pyspi.calculator import Calculator
from pyunicorn.timeseries import Surrogates

# Paramecium and Didinium (empirical predator prey)
# data from https://doi.org/10.2307/4195

def Paramecium_Didinium_0_375(CC = '0.375', timesteps = 1000, **kwargs):
    kwargs.setdefault('data')
    # read from predator_prey_data/CC0.375.txt
    data = pd.read_csv(f'predator_prey_data/CC{CC}.txt', sep = ',')

    if timesteps < len(data):
        return data['t'].values[0:timesteps], data['x'].values[0:timesteps], data['y'].values[0:timesteps]
    else:
        return data['t'].values, data['x'].values, data['y'].values

def Paramecium_Didinium_0_5(CC = '0.5', timesteps = 1000, **kwargs):
    kwargs.setdefault('data')
    # read from predator_prey_data/CC0.375.txt
    data = pd.read_csv(f'predator_prey_data/CC{CC}.txt', sep = ',')
    if timesteps < len(data):
        return data['t'].values[0:timesteps], data['x'].values[0:timesteps], data['y'].values[0:timesteps]
    else:
        return data['t'].values, data['x'].values, data['y'].values

def chaotic_system(corr_yx, corr_xy=0, timesteps = 1000, dt = 1, noise = 0):
    x = np.array([.5])
    y = np.array([.5])
    t = np.arange(0, timesteps, dt)
    for i in range(1, timesteps):
        xtp1 = x[i-1]*(3.8-3.8*x[i-1] - corr_xy*y[i-1]) + np.random.normal(0, noise)
        ytp1 = y[i-1]*(3.5-3.5*y[i-1] - corr_yx*x[i-1]) + np.random.normal(0, noise)
        x = np.append(x, xtp1)
        y = np.append(y, ytp1)
    return t, x, y

def CAM(corr_yx, corr_xy = 0, timesteps = 1000, dt = 1, noise = 0.3):
    x = np.array([0.5])
    y = np.array([0.5])
    t = np.arange(0, timesteps, dt)
    for i in range(1, timesteps):
        xtp1 = 0.5*x[i-1] + corr_xy*y[i-1] + np.random.normal(0, noise)
        ytp1 = corr_yx*x[i-1] + 0.7*y[i-1] + np.random.normal(0, noise)
        x = np.append(x, xtp1)
        y = np.append(y, ytp1)
    return t, x, y

def CAM_external(corr_yx = 0, corr_xy = 0.2, timesteps = 1000, dt = 1, noise = 0.3):
    x = np.array([0.5])
    y = np.array([0.5])
    t = np.arange(0, timesteps, dt)
    force = np.sin(t/10) + np.sin(t/30)
    for i in range(1, timesteps):
        xtp1 = 0.5*x[i-1] + corr_xy*y[i-1] + np.random.normal(0, noise) + 0.1 * force[i]
        ytp1 = corr_yx*x[i-1] + 0.7*y[i-1] + np.random.normal(0, noise) + 0.3 * force[i]
        x = np.append(x, xtp1)
        y = np.append(y, ytp1)
    return t, x, y

def method_CCM(x, y, **kwargs):
    pair = pd.DataFrame({'x': x, 'y': y})

    dfE1 = edm.EmbedDimension(dataFrame=pair, lib = [1, int(len(pair)/2)], pred = [int(len(pair)/2), len(pair)],
                                columns='x', target='x', showPlot=False, maxE=6, tau = 1)
    dfE2 = edm.EmbedDimension(dataFrame=pair, lib = [1, int(len(pair)/2)], pred = [int(len(pair)/2), len(pair)],
                                columns='y', target='y', showPlot=False, maxE=6, tau = 1)
    dfE = pd.concat([dfE1, dfE2])
    E = int(dfE.iloc[dfE['rho'].idxmax()]['E'])
    ccm = edm.CCM(dataFrame=pair, E=E, tau=1, columns='x', target='y', libSizes=[len(pair)], sample=1, showPlot=False)

    return ccm['y:x'].iloc[-1]

def method_TE(x, y, **kwargs):
    k = kwargs.get('k', 1)
    embedding = kwargs.get('embedding', 1)
    return te.te_compute(y, x, safetyCheck=False, k=k, embedding=embedding)

def method_TE_kraskov_k1(x, y, **kwargs):
    kwargs.setdefault('data')
    dataset = np.array([x, y])

    print("dataset:   ", dataset)
    print("dataset shape:   ", dataset.shape)
    print("dataset dtype:   ", dataset.dtype)
    calc = Calculator(dataset = dataset, configfile='TE_configs/kraskov_k1.yaml')
    calc.compute()
    
    results = calc.table
    print(results)

    for column in results.columns.levels[0]:
        return results[column]['proc-1']['proc-0']
    
    raise ValueError("No results found in the table")


def method_TE_kraskov_DCE(x, y, **kwargs):
    kwargs.setdefault('data')
    dataset = np.array([x, y])

    print("dataset:   ", dataset)
    print("dataset shape:   ", dataset.shape)
    print("dataset dtype:   ", dataset.dtype)
    calc = Calculator(dataset = dataset, configfile='TE_configs/kraskov_DCE.yaml')
    calc.compute()
    
    results = calc.table
    print(results)

    for column in results.columns.levels[0]:
        return results[column]['proc-1']['proc-0']
    
    raise ValueError("No results found in the table")

def method_TE_kraskov(x, y, **kwargs):
    kwargs.setdefault('data')
    dataset = np.array([x, y])
    calc = Calculator(dataset = dataset, configfile='TE_configs/kraskov.yaml')
    calc.compute()
    
    results = calc.table
    print(results)

    for column in results.columns.levels[0]:
        return results[column]['proc-1']['proc-0']
    
    raise ValueError("No results found in the table")

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
            

def compute_pvalues_parallel(x1, x2, x3, y, causal_method, **kwargs):
    # corr1 = causal_method(x1, y, **kwargs)
    # corr2 = causal_method(x2, y, **kwargs)
    corr3 = causal_method(x3, y, **kwargs)
    # return corr1, corr2, corr3
    return corr3
def obtain_pvalues(t, x, y, causal_method, nshuff = 10,  **kwargs): # arguments for the causal method
    
    corr = causal_method(x, y, **kwargs)

    lfrom = np.random.randint(low=0, high=len(x), size=nshuff)

    futures = {}
    for i in range(nshuff):
        x1 = np.random.permutation(x)
        x2 = np.concatenate((x[lfrom[i]:], x[:lfrom[i]]))
        # x3 = Surrogates(np.array([x])).correlated_noise_surrogates()[0]
        x3 = stationary_bootstrap(x)
        # y3 = stationary_bootstrap(y)
        # x3 = x1                                   # Y IS ALSO BOOTSTRAPPED
        futures[i] = compute_pvalues_parallel(x1, x2, x3, y, causal_method, **kwargs)

    
    pval1 = 0
    pval2 = 0
    pval3 = 0
    for i in range(nshuff):
        # corr1, corr2 = futures[i].result()
        # if corr1 > corr:
        #     pval1 = pval1 + 1
        # if corr2 > corr:
        #     pval2 = pval2 + 1
        corr3 = futures[i]
        if corr3 > corr:
            pval3 = pval3 + 1
    
    # pval1 = pval1/nshuff
    # pval2 = pval2/nshuff
    pval3 = pval3/nshuff

    # return corr, pval1, pval2
    return corr, pval3

def signif(pval, alpha = 0.05):
    if pval < alpha:
        return f"**{pval:.4f}**"
    else:
        return f"{pval:.4f}"

if __name__ == '__main__':

    systems = {
               "chaotic_system":chaotic_system,
               "CAM":CAM,
               "CAM (external common influence)":CAM_external,
               "Paramecium_Didinium_0.375":Paramecium_Didinium_0_375,
               "Paramecium_Didinium_0.5":Paramecium_Didinium_0_5
               }
    
    # Table for different entropies

    methods = {
            #    "TE":method_TE,
               "TE_kraskov_k1":method_TE_kraskov_k1
            #    "TE_kraskov_DCE":method_TE_kraskov_DCE,
            #    "TE_kraskov":method_TE_kraskov
               }
    
    # methods = {"TE":method_TE,
    #            "CCM":method_CCM}

    nshuff = 10000
    # system = "chaotic_system"

# methods in columns

    # open file
    with open('pvalues_bootstrap_source.txt', 'w') as f:
        for system in systems:
            t,x,y = systems[system](corr_yx = 0.4, timesteps = 200, dt = 1, noise = 0.01)
            f.write("|   "+system+"   |   "+r"$x\ \rarr \ y$"+"   |   IDK   |   ")
            for method in methods:
                corr, pval3 = obtain_pvalues(t, x, y, methods[method], nshuff = nshuff, k = 1, embedding = 1)
                f.write(signif(pval3)+"   |   ")
            f.write("\n")

            f.write("|   "+system+"   |   "+r"$y\ \rarr \ x$"+"   |   IDK   |   ")
            for method in methods:
                corr, pval3 = obtain_pvalues(t, y, x, methods[method], nshuff = nshuff, k = 1, embedding = 1)
                f.write(signif(pval3)+"   |   ")
            f.write("\n")
    exit()

# methods in rows
    # open file
    with open('NAHBROSTOP.txt', 'w') as f:
        for system in systems:
            for method in methods:
                t,x,y = systems[system](corr_yx = 0.4, timesteps = 200, dt = 1, noise = 0.01)
                
                corr, pval3 = obtain_pvalues(t, x, y, methods[method], nshuff = nshuff, k = 1, embedding = 1)
                print("|   ",system, "   |   ", method, "   |   ", r"$x\ \rarr \ y$", "   |   IDK   |   ", f"{corr:.3f}", "   |   ", signif(pval3), "   |")
                # write to file
                f.write("|   "+system+"   |   "+method+"   |   "+r"$x\ \rarr \ y$"+"   |   IDK   |   "+f"{corr:.3f}"+"   |   "+signif(pval3)+"   |\n")
                corr, pval3 = obtain_pvalues(t, y, x, methods[method], nshuff = nshuff, k = 1, embedding = 1)
                print("|   ",system, "   |   ", method, "   |   ", r"$y\ \rarr \ x$", "   |   IDK   |   ", f"{corr:.3f}", "   |   ", signif(pval3), "   |")
                # write to file
                f.write("|   "+system+"   |   "+method+"   |   "+r"$y\ \rarr \ x$"+"   |   IDK   |   "+f"{corr:.3f}"+"   |   "+signif(pval3)+"   |\n")
        # close file
        f.close()