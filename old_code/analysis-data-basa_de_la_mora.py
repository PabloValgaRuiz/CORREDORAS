
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pygam as gm
import copy as cp

        
if __name__ == '__main__':

    dirdatain = '../data/'

    species = {} # dictionary with the IDs of the species (as they change across networks)

    print('Common settings done!')


    mydf = pd.read_csv(dirdatain+'bsm_raw_pollen.csv')
    mydf = mydf.fillna(0)
    mydf['age'] = -mydf['age']
    mydf = mydf.sort_values(by=['age'])

    # total counts of pollen for total abundances in conditional GC
    mydf_total_pollen = pd.read_csv(dirdatain+'basa_raw_counts.csv')
    mydf_total_pollen = mydf_total_pollen.fillna(0)
    mydf_total_pollen['age'] = -mydf_total_pollen['age']
    mydf_total_pollen = mydf_total_pollen.sort_values(by=['age'])

    # rename the columns: remove the spaces and replace them with underscores

    mydf.columns = [col.replace(' ', '_') for col in mydf.columns]
    mydf_total_pollen.columns = [col.replace(' ', '_') for col in mydf_total_pollen.columns]

    all_my_species = list(mydf.columns)
    myspecies = all_my_species[1:] # ALL SPECIES
    myspecies = 'Abies,Pinus,Juniperus,Taxus,Betula,Corylus,Alnus,Carpinus,Salix,Ulmus,Populus,Acer,Fraxinus,Fagus,Tilia,Juglans,Castanea,Quercus caducifolio,Quercus perennifolio,Pistacia,Rhamnus,Phillyrea,Buxus,Sambucus,Viburnum,Sanguisorba,Tamarix,Thymelaeaceae,Ephedra distachya,Ephedra fragilis,Ericaceae,Hereda helix,Ilex aquifolium,Viscum album,Lonicera,Vitis,Oleaceae,Myrtus,Olea,Poaceae,Lygeum spartum,Artemisia,Cichorioideae,Asteroideae,Cardueae,Rubiaceae,Centaurea,Chenopodiaceae,Caryophyllaceae,Plantago,Brassicaceae,Saxifragaceae,Fabaceae,Genista,Lotus type,Trifolium type,Rosaceae,Ribes,Boraginaceae,Sedum,Helianthemum,Lamiaceae,Urticaceae,Rumex,Berberidaceae,Euphorbiaceae,Primulaceae,Scrophulariaceae,Papaver,Campanulaceae,Convolvulaceae,Liliaceae,Iridaceae,Crassulaceae,Ranunculaceae,Cistaceae,Galium,Apiaceae,Valerianaceae,Cerealia type,Polygonaceae,Ranunculus'.replace(' ','_').split(',')

    # print which species in myspecies are not in mydf_total_pollen
    for spec in myspecies[:]:
        if spec not in mydf_total_pollen.columns:
            print('Species %s is not in mydf_total_pollen!' %(spec))
            myspecies.remove(spec)
    # exit(0)


    mydf_total_pollen = mydf_total_pollen.rename(columns={'Quercus_(decidious)': 'Quercus_caducifolio'})
    mydf_total_pollen = mydf_total_pollen.rename(columns={'Quercus_(evergreen)': 'Quercus_perennifolio'})
  # mydf_total_pollen = mydf_total_pollen.rename(columns={'': 'Sanguisorba'})
    mydf_total_pollen = mydf_total_pollen.rename(columns={'Hedera_helix': 'Hereda_helix'}) # it's hedera helix, but in the pollen abundances dataset it's written as hereda helix
  # mydf_total_pollen = mydf_total_pollen.rename(columns={'':'Lygeum_spartum'})
    mydf_total_pollen = mydf_total_pollen.rename(columns={'Chenopodiaceae':'Amaranthaceae'})
    mydf_total_pollen = mydf_total_pollen.rename(columns={'Lotus-type':'Lotus_type'})
    mydf_total_pollen = mydf_total_pollen.rename(columns={'Trifolium-type':'Trifolium_type'})
    mydf_total_pollen = mydf_total_pollen.rename(columns={'Valerianoideae':'Valerianaceae'})
    mydf_total_pollen = mydf_total_pollen.rename(columns={'Cerealia-type':'Cerealia_type'})


    for spec in myspecies[:]: # ITERATE OVER A COPY
        if mydf[spec].sum() == 0:
            print('Species %s has all zeros!' %(spec))
            myspecies.remove(spec)
            mydf = mydf.drop(columns=[spec])
            mydf_total_pollen = mydf_total_pollen.drop(columns=[spec])

    # sum all of the abundances over time
    sum = mydf_total_pollen[myspecies].sum(axis=1)
    # plot the sum of the abundances over time
    plt.figure(figsize=(10, 5))
    plt.plot(mydf['age'], sum, label='Total abundances', color='black')
    plt.savefig('plots-tests/basa_total_abundances_TESTS.pdf', dpi=300)

    exit(0)
    
#____________________________________________________ GAM ____________________________________________________

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

        print('Ho finito!\n')

        species_df = pd.DataFrame(species)
        species_df.to_pickle('GAM_species/species_%s.pkl' %(fit_type))
    