
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pygam as gm
import copy as cp


def clean_basa_pollen_concentration(mydf_pollen, myspecies):
    mydf_pollen['cal BP'] = - mydf_pollen['cal BP']
    mydf_pollen = mydf_pollen.fillna(0)
    mydf_pollen = mydf_pollen.sort_values(by=['cal BP'])

    all_species = list(mydf_pollen.columns)

    # create a dictionary to map conflicting species names
    species_mapping = {
        "Dec_Querc": "Quercus caducifolio",
        "Ever_Querc": "Quercus perennifolio",
        "Ephedra dist": "Ephedra distachya",
        "Ephedra frag": "Ephedra fragilis",
        "Lygeum": "Lygeum spartum",
        "Cicho": "Cichorioideae",
        "Astroi": "Asteroideae",
        "Carduaceae": "Cardueae",
        "Rubiac": "Rubiaceae",
        "Chenopo": "Chenopodiaceae",
        "Caryphy": "Caryophyllaceae",
        "Brassicac": "Brassicaceae",
        "Saxifrag": "Saxifragaceae",
        "Boraginac": "Boraginaceae",
        "Helianthem": "Helianthemum",
        "Euphorbiac": "Euphorbiaceae",
        "Primulac": "Primulaceae",
        "Scrophulari": "Scrophulariaceae",
        "Campanulac": "Campanulaceae",
        "Valerian": "Valerianaceae",
        "Cerealia": "Cerealia type",
        "Polygon": "Polygonaceae"
    }
    # Sanguisorba and Ribes do not appear. Sanguisorba has 0 pollen, Ribes is not zero but not in the list
    myspecies.remove('Sanguisorba')
    myspecies.remove('Ribes')

    # replace the column names in mydf_pollen_concentration
    mydf_pollen = mydf_pollen.rename(columns=species_mapping)

    metrics_list = ['accrate', 'volume', 'lycadd', 'lyc']

    mydf_pollen = mydf_pollen[['cal BP'] + metrics_list + myspecies] # not taking into account the lycadd / lyc factor
    mydf_pollen = mydf_pollen.set_index('cal BP')

    multipliers = mydf_pollen['accrate'] / mydf_pollen['volume']
    # print(mydf_pollen[myspecies])

    mydf_pollen_PAR = mydf_pollen[myspecies].copy()
    for spec in myspecies:
        mydf_pollen_PAR[spec] = mydf_pollen[spec] * multipliers
    
    # pickle the multipliers
    pd.to_pickle(multipliers, 'results_basa/multipliers.pkl')

    return mydf_pollen[myspecies] # mydf_pollen or mydf_pollen_PAR

def plot_abundances(mydf = pd.DataFrame, myspecies = list, dir = 'plots-tests/basa_total_abundances_TESTS.pdf'):
        # sort myspecies by highest aboundance
    myspecies = sorted(myspecies, key=lambda x: mydf[x].sum(), reverse=True)

    sum = mydf[myspecies].sum(axis=1)
    # plot the sum of the abundances over time
    plt.figure(figsize=(10, 5))
    # choose a different mpl color sequence
    print(plt.style.available)
    plt.style.use('classic')
    plt.plot(mydf.index, mydf[myspecies])
    plt.xlim(mydf.index.min(), mydf.index.max())
    plt.legend(myspecies[:15], loc='upper right', fontsize=7)
    plt.grid(True)
    plt.xlabel('cal BP')
    plt.ylabel('Total abundances')
    plt.savefig(dir, dpi=300)


if __name__ == '__main__':

    dirdatain = '../data/'    

    myspecies = 'Abies,Pinus,Juniperus,Taxus,Betula,Corylus,Alnus,Carpinus,Salix,Ulmus,Populus,Acer,Fraxinus,Fagus,Tilia,Juglans,Castanea,Quercus caducifolio,Quercus perennifolio,Pistacia,Rhamnus,Phillyrea,Buxus,Sambucus,Viburnum,Sanguisorba,Tamarix,Thymelaeaceae,Ephedra distachya,Ephedra fragilis,Ericaceae,Hereda helix,Ilex aquifolium,Viscum album,Lonicera,Vitis,Oleaceae,Myrtus,Olea,Poaceae,Lygeum spartum,Artemisia,Cichorioideae,Asteroideae,Cardueae,Rubiaceae,Centaurea,Chenopodiaceae,Caryophyllaceae,Plantago,Brassicaceae,Saxifragaceae,Fabaceae,Genista,Lotus type,Trifolium type,Rosaceae,Ribes,Boraginaceae,Sedum,Helianthemum,Lamiaceae,Urticaceae,Rumex,Berberidaceae,Euphorbiaceae,Primulaceae,Scrophulariaceae,Papaver,Campanulaceae,Convolvulaceae,Liliaceae,Iridaceae,Crassulaceae,Ranunculaceae,Cistaceae,Galium,Apiaceae,Valerianaceae,Cerealia type,Polygonaceae,Ranunculus'.split(',')
    
    mydf = pd.read_csv(dirdatain+'basa_original.csv') #'basa_char_par.csv'
    mydf = clean_basa_pollen_concentration(mydf, myspecies=myspecies)

    for spec in myspecies[:]: # ITERATE OVER A COPY
        if mydf[spec].sum() == 0:
            print('Species %s has all zeros!' %(spec))
            myspecies.remove(spec)
            mydf = mydf.drop(columns=[spec])

    # plot_abundances(mydf, myspecies)
    
#____________________________________________________ GAM ____________________________________________________
    n_splines = 140
    fit_type = f'fit_{n_splines}_0.01'
    if 1: # DO NOT DO GAM IF IT ALREADY EXISTS
        print('\nComputing GAM ....')
        # for transfer entropy

        species = {}
        vx = gm.utils.make_2d(mydf.index, verbose=False).astype('float')
        for spec in myspecies:
            print('\tIterating on species %s ...' %(spec))
            vy = mydf[spec]

            mygam = gm.GAM(terms='auto', n_splines=n_splines, lam=0.01).fit(vx, vy)
            # mygam = gm.GAM(terms='auto').gridsearch(vx, vy)

            species[spec] = {}
            species[spec]['prev_x'] = mydf.index.to_numpy()
            species[spec]['prev_y'] = mydf[spec].to_numpy()

            XX = mygam.generate_X_grid(term=0, n=140)
            YY = mygam.predict(XX)
            species[spec]['x'] = cp.copy(XX)[:,0]
            species[spec]['y'] = cp.copy(YY)

        species_df = pd.DataFrame(species)
        species_df.to_pickle('GAM_species/species_%s.pkl' %(fit_type))

#___________________________________ GAM over Mendukilo (delta13 -> temperatures) ___________________________________

dirdatain = '../data/'
delta13 = pd.read_csv(f'{dirdatain}/Mendukilo.csv')
delta13['age (kyr BP)'] = - delta13['age (kyr BP)'] * 1000  # Convert from kyr to years
delta13 = delta13.fillna(0)
delta13 = delta13.sort_values(by=['age (kyr BP)'])

vx = gm.utils.make_2d(delta13['age (kyr BP)'], verbose=False).astype('float')
vy = delta13['d13C (permil)']

d13gam = gm.GAM(terms='auto', n_splines=n_splines, lam=0.01).fit(vx, vy)

YY = d13gam.predict(XX)

d13gam_df = pd.DataFrame({
    'age (kyr BP)': XX[:, 0],
    'd13C (permil)': YY
})

d13gam_df.to_pickle('GAM_species/d13gam_%s.pkl' %(fit_type))

#__________________________________________ GAM over multipliers ______________________________________________________

multipliers = pd.read_pickle('results_basa/multipliers.pkl')
vx = gm.utils.make_2d(multipliers.index, verbose=False).astype('float')
vy = multipliers.values

gam_multipliers = gm.GAM(terms='auto', n_splines=n_splines, lam=0.01).fit(vx, vy)

YY = gam_multipliers.predict(XX)

gam_multipliers_df = pd.DataFrame({
    'cal BP': XX[:, 0],
    'multipliers': YY
})

gam_multipliers_df.to_pickle('GAM_species/gam_multipliers_%s.pkl' %(fit_type))