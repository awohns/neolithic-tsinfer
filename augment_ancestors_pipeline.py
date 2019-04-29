import pandas as pd
import numpy as np
import json
import collections 
import time
import subprocess
import argparse
import os

import msprime
import tskit
import tsinfer



def get_sites_alleles(samples, dataset,chromosome, path):
    """
    Create list of sites and alleles from the sgdp samples file
    """
    sites_alleles = pd.DataFrame(index = range(0,samples.num_sites), columns= ["Chromosome","Position","Ref","Alt"])

    sites_alleles.Position = [int(site) for site in samples.sites_position[:]]
    sites_alleles.Chromosome = chromosome
    sites_alleles.to_csv(path + dataset + "_chr" + chromosome + "_sites.txt",sep=" ",index=False, header=False)
    sites_alleles.Ref = [allele[0] for allele in samples.sites_alleles[:]] 
    sites_alleles.Alt = [allele[1] for allele in samples.sites_alleles[:]] 

    sites_alleles.to_csv(path + dataset + "_chr" + chromosome + "_sites_alleles.txt",sep=" ")


def assign_genotype(tm611, tm613, aligned): 
    genos_aligned = [[0,1],[1,0],[0,0],[1,1]]
    genos_reversed = [[1,0],[0,1],[1,1],[0,0]]

    if aligned:
        genos = genos_aligned
    elif not aligned:
        genos = genos_reversed

    if tm611.TM611.values[0][0:3] == "0|1":
        tm611_geno = genos[0]
    elif tm611.TM611.values[0][0:3] == "1|0":
        tm611_geno = genos[1]
    elif tm611.TM611.values[0][0:3] == "0|0":
        tm611_geno = genos[2]
    elif tm611.TM611.values[0][0:3] == "1|1":
        tm611_geno = genos[3]
        
    if tm613.TM613.values[0][0:3] == "0|1":
        tm613_geno = genos[0]
    elif tm613.TM613.values[0][0:3] == "1|0":
        tm613_geno = genos[1]
    elif tm613.TM613.values[0][0:3] == "0|0":
        tm613_geno = genos[2]
    elif tm613.TM613.values[0][0:3] == "1|1":
        tm613_geno = genos[3]

    return(tm611_geno + tm613_geno)

def create_neolithic_samples_file(dataset, chromosome, path):
    #Match the neolithic samples to the UKBB file
    tm611_df = pd.read_csv(path + 'neolithic_vcfs/TM611_' + dataset + '_chr'+ chromosome + '_sites.recode.vcf', sep="\t",skiprows=34)
    tm613_df = pd.read_csv(path + 'neolithic_vcfs/TM613_' + dataset + '_chr'+ chromosome + '_sites.recode.vcf', sep="\t",skiprows=34)

    sites = pd.read_csv(path + "augment_ancestors/neolithic_sites_alleles_overlap/" + dataset + "_chr" + chromosome + "_sites_alleles.txt", sep=" ")
       

    #Create sample_file
    with tsinfer.SampleData(path=path + "augment_ancestors/neolithic_sample_data/" + dataset + "_chr" + chromosome + "_neolithic.samples") as sample_data:
        # Define populations
        sample_data.add_population(metadata={"name": "Neolithic"})
        # Define individuals
        sample_data.add_individual(
            ploidy=2, population=0, metadata={"name": "TM611"})
        sample_data.add_individual(
            ploidy=2, population=0, metadata={"name": "TM613"})

        # Define sites and genotypes
        for index, site in sites.iterrows():
            if index % 10000 ==0:
                print("Adding sample data file: ",index," of ", len(sites.index))
            if (site.Position in tm611_df.POS.values) and (site.Position in tm613_df.POS.values):
                tm611_value = tm611_df[tm611_df.POS == site.Position]
                tm613_value = tm613_df[tm613_df.POS == site.Position]

                ref = tm611_value.REF.values[0]
                alt = tm613_value.ALT.values[0]

                #If the refs and alts line up between aDNA and UKBB
                if (ref == site.Ref and alt == site.Alt):
                    sample_data.add_site(site.Position, assign_genotype(tm611_value, tm613_value, True), [ref, alt])
                    
                elif (ref == site.Alt and alt == site.Ref):
                    sample_data.add_site(site.Position, assign_genotype(tm611_value, tm613_value, False), [alt, ref])
    return(sample_data)

def intersect_samples_ancestors_ts(ancient_samples, ancestors_ts, output_name, path, inference_sites_only = True):
    tables = ancestors_ts.dump_tables()
    sample_data_position = ancient_samples.sites_position[:]
    if inference_sites_only:
        neolithic_sites = set(sample_data_position[
           ancient_samples.sites_inference[:] == 1])
    else:
        neolithic_sites = set(sample_data_position)
    ancestors_sites = set(tables.sites.position[:])
    intersecting_sites = ancestors_sites & neolithic_sites
    
    print("Intersecting sites = ", len(intersecting_sites))
    tables.sites.clear()
    tables.mutations.clear()
    for site in ancestors_ts.sites():
        if site.position in intersecting_sites:
            # Sites must be 0/1 for the ancestors ts.
            site_id = tables.sites.add_row(
                position=site.position, ancestral_state=site.ancestral_state)
            try:
                assert len(site.mutations) == 1
            except AssertionError:
                print(site)
            mutation = site.mutations[0]
            tables.mutations.add_row(
                site=site_id, node=mutation.node, derived_state="1")
    # Reduce this to the site topology now to make things as quick as possible.
    tables.simplify(reduce_to_site_topology=True, filter_sites=False)
    reduced_ts = tables.tree_sequence()
    reduced_ts.dump(path + "/" + output_name + ".trees")

    return(reduced_ts)


def downsample_dataset_samples(dataset_samples, augmented_ancestors_ts, dataset, chromosome, path):
    dataset_sample_positions = set(dataset_samples.sites_position[:])
    tables = augmented_ancestors_ts.dump_tables()
    ancestors_sites = set(tables.sites.position[:])
    intersecting_sites = ancestors_sites & dataset_sample_positions

    with tsinfer.SampleData(
            path=path+"augment_ancestors/downsampled_dataset_sample_data/" + dataset + "_chr" + chromosome + "_downsampled.samples", sequence_length=dataset_samples.sequence_length) as sample_data:
        for ind in dataset_samples.individuals():
            sample_data.add_individual(ploidy=2, location=ind.location, metadata=ind.metadata)

        for var in dataset_samples.variants():
            if var.site.position in intersecting_sites:
                sample_data.add_site(var.site.position, var.genotypes, var.alleles)
    return(sample_data)


if __name__ == "__main__":
    path = "/gpfs0/well/mcvean/ukbb12788/wohns/"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset', '-d', type=str, default="1kg", help="1kg, ukbb, or sgdp")
    parser.add_argument('--chromosome', '-c', type=str, default=20, help="chromosome number")

    args = parser.parse_args()
    dataset_samples = tsinfer.load("/well/mcvean/ukbb12788/wohns/treeseq-inference/human-data/" + args.dataset + "_chr" + args.chromosome + ".samples")
    get_sites_alleles(dataset_samples, args.dataset, args.chromosome, path + "augment_ancestors/neolithic_sites_alleles_overlap/")

    if not os.path.isfile(path + "neolithic_vcfs/TM611_" + args.dataset + "_chr" + args.chromosome + "_sites.recode.vcf"):
        subprocess.call("vcftools --gzvcf " + path + "neolithic_vcfs/TM611prefiltered.vcf.gz --positions " + path + "augment_ancestors/neolithic_sites_alleles_overlap/" + args.dataset + "_chr" + args.chromosome + "_sites.txt --recode --recode-INFO-all --out " + path + "neolithic_vcfs/TM611_" + args.dataset + "_chr" + args.chromosome + "_sites", shell=True)
    if not os.path.isfile(path + "neolithic_vcfs/TM613_" + args.dataset + "_chr" + args.chromosome + "_sites.recode.vcf"):
        subprocess.call("vcftools --gzvcf " + path + "neolithic_vcfs/TM613prefiltered.vcf.gz --positions " + path + "augment_ancestors/neolithic_sites_alleles_overlap/" + args.dataset + "_chr" + args.chromosome + "_sites.txt --recode --recode-INFO-all --out " + path + "neolithic_vcfs/TM613_" + args.dataset + "_chr" + args.chromosome + "_sites", shell=True)

    print("creating samples file")

    if not os.path.isfile(path + "augment_ancestors/neolithic_sample_data/" + args.dataset + "_chr" + args.chromosome + "_neolithic.samples"):
        neolithic_samples = create_neolithic_samples_file(args.dataset, args.chromosome, path)
    else:
        neolithic_samples = tsinfer.load(path + "augment_ancestors/neolithic_sample_data/" + args.dataset + "_chr" + args.chromosome + "_neolithic.samples")

    print("Loading Ancestors Tree Sequence")

    ancestors_ts = msprime.load("/well/mcvean/ukbb12788/wohns/treeseq-inference/human-data/" + args.dataset + "_chr" + args.chromosome + ".ancestors.trees")

    neolithic_samples_copy = neolithic_samples.copy()
    neolithic_samples_copy.sites_inference[:] = 1

    if not os.path.isfile(path + "augment_ancestors/intersected_ancestors_ts/" + args.dataset + "_chr" + args.chromosome + "_reduced.trees"):
        reduced_ancestors_ts = intersect_samples_ancestors_ts(neolithic_samples_copy, ancestors_ts, args.dataset + "_chr" + args.chromosome + "_reduced", path + "augment_ancestors/intersected_ancestors_ts/", False)
    else:
        reduced_ancestors_ts = tskit.load(path + "augment_ancestors/intersected_ancestors_ts/" + args.dataset + "_chr" + args.chromosome + "_reduced.trees")

    if not os.path.isfile(path + "augment_ancestors/augmented_ancestors_ts/" + args.dataset + "_chr" + args.chromosome + "_augmented_ancestors.trees"):
        # subprocess.call("tsinfer augment-ancestors -A " + path + "augment_ancestors/intersected_ancestors_ts/" + args.dataset + "_chr" + args.chromosome + "_reduced.trees " + path + "augment_ancestors/neolithic_sample_data/" + args.dataset + "_chr" + args.chromosome + "_neolithic.samples " + path + "augment_ancestors/augmented_ancestors_ts/" + args.dataset + "_chr" + args.chromosome + "_augmented_ancestors.trees", shell=True)
        augmented_ancestors_ts = tsinfer.augment_ancestors(neolithic_samples_copy, reduced_ancestors_ts, [0,1,2,3])
        augmented_ancestors_ts.dump(path + "augment_ancestors/augmented_ancestors_ts/" + args.dataset + "_chr" + args.chromosome + "_augmented_ancestors.trees")
    else:
        augmented_ancestors_ts = tskit.load(path + "augment_ancestors/augmented_ancestors_ts/" + args.dataset + "_chr" + args.chromosome + "_augmented_ancestors.trees")

    if not os.path.isfile(path+"augment_ancestors/downsampled_dataset_sample_data/" + args.dataset + "_chr" + args.chromosome + "_downsampled.samples"):
        downsampled_dataset_samples = downsample_dataset_samples(dataset_samples, augmented_ancestors_ts, args.dataset, args.chromosome, path)

    if not os.path.isfile(path + "augment_ancestors/matched_dataset_samples_augmented_ts/" + args.dataset + "_chr" + args.chromosome + "_matched_samples_augmented_ts.trees"):

        subprocess.call("nice tsinfer match-samples -A " + path + "augment_ancestors/augmented_ancestors_ts/" + args.dataset + "_chr" + args.chromosome + "_augmented_ancestors.trees --no-simplify " + path + "augment_ancestors/downsampled_dataset_sample_data/" + args.dataset + "_chr" + args.chromosome + "_downsampled.samples -O " + path + "augment_ancestors/matched_dataset_samples_augmented_ts/" + args.dataset + "_chr" + args.chromosome + "_matched_samples_augmented_ts.trees --progress --num-threads 64", shell=True)
 
  