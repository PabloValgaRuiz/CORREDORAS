
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pygam as gm
import copy as cp


def clean_basa_pollen_concentration(mydf_pollen_concentration : pd.DataFrame, myspecies : list) -> pd.DataFrame:
    """
    Cleans and processes a DataFrame containing pollen concentration data.
    This function performs the following operations:
    1. Converts the 'cal BP' column values to negative and sorts them (to get ascending time).
    2. Maps conflicting species names to standardized names using a predefined dictionary.
    3. Removes specific species ('Sanguisorba' and 'Ribes') from the provided species list.
    4. Renames the columns of the DataFrame based on the species mapping.
    5. Sets the 'cal BP' column as the index of the DataFrame.
    Args:
        mydf_pollen_concentration (pd.DataFrame): A DataFrame containing pollen concentration data.
            It must include a 'cal BP' column and other columns representing species.
        myspecies (list): A list of species names to be retained in the DataFrame.
    Returns:
        pd.DataFrame: A cleaned and processed DataFrame with standardized species names and 'cal BP' as the index.
    """
    mydf_pollen_concentration['cal BP'] = - mydf_pollen_concentration['cal BP']
    mydf_pollen_concentration = mydf_pollen_concentration.fillna(0)
    mydf_pollen_concentration = mydf_pollen_concentration.sort_values(by=['cal BP'])

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
    mydf_pollen_concentration = mydf_pollen_concentration.rename(columns=species_mapping)

    mydf_pollen_concentration = mydf_pollen_concentration[['cal BP'] + myspecies]
    mydf_pollen_concentration = mydf_pollen_concentration.set_index('cal BP')

    return mydf_pollen_concentration

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
    
    mydf = pd.read_csv(dirdatain+'basa_char_par.csv') #'bsm_raw_pollen.csv'
    mydf = clean_basa_pollen_concentration(mydf, myspecies=myspecies)

    for spec in myspecies[:]: # ITERATE OVER A COPY
        if mydf[spec].sum() == 0:
            print('Species %s has all zeros!' %(spec))
            myspecies.remove(spec)
            mydf = mydf.drop(columns=[spec])

    # plot_abundances(mydf, myspecies)
    
#____________________________________________________ GAM ____________________________________________________

    fit_type = 'fit_140_0.01'
    if 1: # DO NOT DO GAM IF IT ALREADY EXISTS
        print('\nComputing GAM ....')
        # for transfer entropy

        species = {}
        vx = gm.utils.make_2d(mydf.index, verbose=False).astype('float')
        for spec in myspecies:
            print('\tIterating on species %s ...' %(spec))
            vy = mydf[spec]

            mygam = gm.GAM(terms='auto', n_splines=140, lam=0.01).fit(vx, vy)
            # mygam = gm.GAM(terms='auto').gridsearch(vx, vy)

            species[spec] = {}
            species[spec]['prev_x'] = mydf.index.to_numpy()
            species[spec]['prev_y'] = mydf[spec].to_numpy()

            XX = mygam.generate_X_grid(term=0, n=140)
            YY = mygam.predict(XX)
            species[spec]['x'] = cp.copy(XX)[:,0]
            species[spec]['y'] = cp.copy(YY)

        print('Ho finito!\n')

        species_df = pd.DataFrame(species)
        species_df.to_pickle('GAM_species/species_%s.pkl' %(fit_type))
    