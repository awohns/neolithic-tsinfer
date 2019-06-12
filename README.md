# neolithic-tsinfer
Run modern descendants analyses for "Kinship and burial in a late Neolithic Cambridgeshire site."

## tsinfer with Trumpington Individuals as Potential Ancestors to SGDP and TGP
First, you will need to download the .trees files for all autosomes of the:

Thousand Genomes Project (1kg/TGP): https://zenodo.org/record/3051855#.XQEYm9NKi7M

Simons Genome Diversity Project (SGDP): https://zenodo.org/record/3052359#.XQEYQtNKi7M

You will need to change the paths in infer_trees/augment_ancestors_pipeline.py to point to directory containing the .trees files downloaded from Zenodo. 

Then run (for example): 
```
python augment_ancestors_pipeline.py --dataset 1kg --chromosome 20
```

To create .trees files generated using the Trumpington individuals as potential ancestors. 

## Generate plots
The code implementing this is in a Jupyter Notebook in figures/ directory. 
