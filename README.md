# neolithic-tsinfer
Run modern descendants analyses for "East Anglian Early Neolithic monument burial linked to contemporary Megaliths."

## tsinfer with Trumpington Individuals as Potential Ancestors to SGDP and TGP
First, you will need to clone https://github.com/mcveanlab/treeseq-inference. In treeseq-inference/human-data, use the makefile to create .trees files for all autosomes of the Thousand Genomes Project (1kg/TGP) and Simons Genome Diversity Project (SGDP). 

Then, just run, for example: 
```
python augment_ancestors_pipeline.py --dataset 1kg --chromosome 20 --path your/path/
```

To create .trees files generated using the Trumpington individuals as potential ancestors. 

Be sure that the path you provide is to a directory containing both the neolithic-tsinfer/ and treeseq-inference/ directories.

## Generate plots
The code implementing this is in a Jupyter Notebook in figures/ directory. 
