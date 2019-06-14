"""
Implements full pipeline from Neolitic VCFs and 1kg/SGDP data to inferred
tree sequences with ancient samples Requires three arguments: a dataset (1kg or
sgdp), a chromosome (autosomes only), and the args.path to the
directory containing both neolithic-tsinfer and treeseq-inference.
"""
import pandas as pd
import subprocess
import argparse
import os

import tskit
import tsinfer


def get_sites_alleles(samples, dataset, chromosome, path):
    """
    Create list of sites and alleles from the specified samples file
    """
    sites_alleles = pd.DataFrame(
        index=range(0, samples.num_sites),
        columns=["Chromosome", "Position", "Ref", "Alt"])

    sites_alleles.Position = [int(site) for site in samples.sites_position[:]]
    sites_alleles.Chromosome = chromosome
    sites_alleles.to_csv(args.path + dataset + "_chr" + chromosome +
                         "_sites.txt", sep=" ", index=False, header=False)
    sites_alleles.Ref = [allele[0] for allele in samples.sites_alleles[:]]
    sites_alleles.Alt = [allele[1] for allele in samples.sites_alleles[:]]

    sites_alleles.to_csv(args.path + dataset + "_chr" + chromosome +
                         "_sites_alleles.txt", sep=" ")


def assign_genotype(tm611, tm613, aligned):
    genos_aligned = [[0, 1], [1, 0], [0, 0], [1, 1]]
    genos_reversed = [[1, 0], [0, 1], [1, 1], [0, 0]]

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


def intersect_samples_ancestors_ts(ancient_samples, ancestors_ts,
                                   output_name, path,
                                   inference_sites_only=True):
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
    reduced_ts.dump(args.path + "/" + output_name + ".trees")

    return(reduced_ts)


def downsample_dataset_samples(dataset_samples, augmented_ancestors_ts,
                               dataset, chromosome, path):
    dataset_sample_positions = set(dataset_samples.sites_position[:])
    tables = augmented_ancestors_ts.dump_tables()
    ancestors_sites = set(tables.sites.position[:])
    intersecting_sites = ancestors_sites & dataset_sample_positions
    with tsinfer.SampleData(
        path=args.path + "neolithic-tsinfer/data/"
        "downsampled_dataset_sample_data/" +
        dataset + "_chr" + chromosome +
        "_downsampled.samples",
        sequence_length=dataset_samples.sequence_length) as sample_data:
        for ind in dataset_samples.individuals():
            sample_data.add_individual(ploidy=2,
                                       location=ind.location,
                                       metadata=ind.metadata)

        for var in dataset_samples.variants():
            if var.site.position in intersecting_sites:
                sample_data.add_site(var.site.position,
                                     var.genotypes, var.alleles)
    return(sample_data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset', '-d', type=str, default="1kg",
                        help="1kg, ukbb, or sgdp")
    parser.add_argument('--chromosome', '-c', type=str, default=20,
                        help="chromosome number")
    parser.add_argument('--path', '-p', type=str,
                        default="/gpfs0/well/mcvean/ukbb12788/wohns/")
    args = parser.parse_args()

    # Load the sample_data file for the specified dataset
    dataset_samples = tsinfer.load(args.path +
                                   "treeseq-inference/human-data/" +
                                   args.dataset + "_chr" + args.chromosome +
                                   ".samples")

    print("Loading Neolithic sample_data")
    neolithic_samples = tsinfer.load(
        args.path + "neolithic-tsinfer/data/neolithic_sample_data/" +
        args.dataset + "_chr" + args.chromosome + "_neolithic.samples")

    print("Loading Ancestors Tree Sequence")
    ancestors_ts = tskit.load(
        args.path + "treeseq-inference/human-data/" + args.dataset +
        "_chr" + args.chromosome + ".ancestors.trees")

    neolithic_samples_copy = neolithic_samples.copy()
    neolithic_samples_copy.sites_inference[:] = 1

    if not os.path.isfile(args.path +
                          "neolithic-tsinfer/data/"
                          "intersected_ancestors_ts/" +
                          args.dataset + "_chr" + args.chromosome +
                          "_reduced.trees"):
        reduced_ancestors_ts = intersect_samples_ancestors_ts(
            neolithic_samples_copy, ancestors_ts, args.dataset +
            "_chr" + args.chromosome + "_reduced", args.path +
            "neolithic-tsinfer/data/intersected_ancestors_ts/", False)
    else:
        reduced_ancestors_ts = tskit.load(
            args.path + "neolithic-tsinfer/data/intersected_ancestors_ts/" +
            args.dataset + "_chr" + args.chromosome + "_reduced.trees")

    if not os.path.isfile(args.path +
                          "neolithic-tsinfer/data/"
                          "augmented_ancestors_ts/" +
                          args.dataset + "_chr" + args.chromosome +
                          "_augmented_ancestors.trees"):
        augmented_ancestors_ts = tsinfer.augment_ancestors(
            neolithic_samples_copy, reduced_ancestors_ts,
            [0, 1, 2, 3])
        augmented_ancestors_ts.dump(
            args.path + "neolithic-tsinfer/data/augmented_ancestors_ts/" +
            args.dataset + "_chr" + args.chromosome +
            "_augmented_ancestors.trees")
    else:
        augmented_ancestors_ts = tskit.load(
            args.path + "neolithic-tsinfer/data/augmented_ancestors_ts/" +
            args.dataset + "_chr" + args.chromosome +
            "_augmented_ancestors.trees")

    if not os.path.isfile(args.path +
                          "neolithic-tsinfer/data/"
                          "downsampled_dataset_sample_data/" +
                          args.dataset + "_chr" + args.chromosome +
                          "_downsampled.samples"):
        downsampled_dataset_samples = downsample_dataset_samples(
            dataset_samples, augmented_ancestors_ts, args.dataset,
            args.chromosome, args.path)

    if not os.path.isfile(args.path +
                          "neolithic-tsinfer/data/"
                           "matched_dataset_samples_augmented_ts/" +
                           args.dataset + "_chr" + args.chromosome +
                           "_matched_samples_augmented_ts.trees"):

        subprocess.call("nice tsinfer match-samples -A " + args.path +
                        "neolithic-tsinfer/data/augmented_ancestors_ts/" +
                        args.dataset + "_chr" + args.chromosome +
                        "_augmented_ancestors.trees --no-simplify " +
                        args.path +
                        "neolithic-tsinfer/data/"
                        "downsampled_dataset_sample_data/" +
                        args.dataset + "_chr" + args.chromosome +
                        "_downsampled.samples -O " + args.path +
                        "neolithic-tsinfer/data/"
                        "matched_dataset_samples_augmented_ts/" +
                        args.dataset + "_chr" + args.chromosome +
                        "_matched_samples_augmented_ts.trees"
                        " --progress --num-threads 64", shell=True)
    print('finished')
