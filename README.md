# HOW TO RUN

It is recommended to run using conda (installed through miniconda). Open a terminal on the repository folder.

```
conda create -p ./.conda python=3.11
conda activate ./.conda
pip install -r requirements.txt
```

```
python src/1_preprocess_basa.py
python src/1_preprocess_gg.py
```
```
python src/2_stationary_bootstrap_GC.py 0 51
python src/2_stationary_bootstrap_GC.py 51 85
python src/2_stationary_bootstrap_GC.py 85 140
```

For the plots, additional packages must be installed. If using a jupyter notebook (in vscode, for example), install ipykernel.
```
conda install -p ./.conda ipykernel --update-deps --force-reinstall
```
Then additional packages.
```
pip install networkx leidenalg iplotx
```

Then run the 3_plot_results_BS_GC.ipynb cells.
