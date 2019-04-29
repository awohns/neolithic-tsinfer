# neolithic-tsinfer
Run modern descendants analyses for "Kinship and burial in a late Neolithic Cambridgeshire site."

## tsinfer with Trumpington Individuals as Potential Ancestors to SGDP and TGP
First, you will need to clone https://github.com/mcveanlab/treeseq-inference. In treeseq-inference/human-data, use the makefile to create .trees files for all autosomes of the Thousand Genomes Project (1kg/TGP) and Simons Genome Diversity Project (SGDP). 

You will need to change the paths in infer_trees/augment_ancestors_pipeline.py to point to the treeseq-inference/human-data directory containing the .trees files created using the makefile. 

Then, just run, for example: 
```
python augment_ancestors_pipeline.py --dataset 1kg --chromosome 20
```

To create .trees files generated using the Trumpington individuals as potential ancestors. 

## Generate plots
The code implementing this is in a Jupyter Notebook in figures/ directory. 
