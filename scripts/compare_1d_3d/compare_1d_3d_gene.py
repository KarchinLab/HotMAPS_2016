import argparse
import csv
import re
import time

# import modules needed for logging
import logging
import os
import random

logger = logging.getLogger(__name__)  # module logger
DROPPED_TTYPES = ["COAD", "READ", "THYM", "CHOL", "UVM", "TGCT", "SARC"]

def parse_arguments():
    """
    Function to parse command line arguements
    from the user

    Returns
    -------
    opts : dict
        command line arguements from the user
    """
    info = 'Compares the 3D clustering results to 1D'
    parser = argparse.ArgumentParser(description=info)
    # program arguments
    parser.add_argument('-d', '--regions-dir',
                        type=str, required=True,
                        help='Directory containing the regions .txt files for 1D and 3D')
    parser.add_argument('-r', '--radius-1d',
                        type=str, required=True,
                        help='The number of codons both upstream and downstream used for 1D algorithm')
    parser.add_argument('-o', '--output-dir',
                        type=str, required=True,
                        help='Directory to write the output to')
    parser.add_argument('-q', '--q-value',
                        type=float, default=.01,
                        help='Significance level of q-value in txt files to analyse')
    args = parser.parse_args()
    opts = vars(args)
    return opts


def read_regions_file(filepath):
    '''
    Convenience function to read a regions file
    and return its contents as a dictionary
    in the form tumour types as keys, pointing to
    dictionaries with genes as keys and a list of hotspot
    regions in these genes as values
    '''

    f = open(filepath)

    gene_or_struct_ind = 0
    ttype_ind = 1

    file_dict = {}

    for line in f:
        line = (line.strip()).split('\t')
        gene_or_struct = line[gene_or_struct_ind]
        ttype = line[ttype_ind]

        if ttype in DROPPED_TTYPES:
            continue

        if not (ttype, gene_or_struct) in file_dict:
            file_dict[(ttype, gene_or_struct)] = []
        #if not gene_or_struct in file_dict[ttype]:
        #    file_dict[ttype][gene_or_struct] = []

        for region in line[2:]:
            file_dict[(ttype, gene_or_struct)].append(region)

    return file_dict


def main(opts):
    """
    Main function
    """
    # number of codons that were examined up/downstream
    # for the 1D hotspot algorithm
    num_codons = opts['radius_1d']

    # get specified parameters
    output_dir = opts['output_dir']
    regions_dir = opts['regions_dir']
    q_val = opts['q_value']
    q_val = (str(q_val)).split('.')[1]

    # get all the required filenames
    regions_gene_1d = regions_dir + "/hotspot_regions_gene_." + q_val + "_1d_{0}.txt".format(num_codons)
    regions_gene_3d = regions_dir + "/hotspot_regions_gene_." + q_val + ".txt"
    #regions_struct_1d = regions_dir + "/hotspot_regions_structure_." + q_val + "_1d.txt"
    #regions_struct_3d = regions_dir + "/hotspot_regions_structure_." + q_val + ".txt"

    gene_1d_dict = read_regions_file(regions_gene_1d)
    gene_3d_dict = read_regions_file(regions_gene_3d)
    #struct_1d_dict = read_regions_file(regions_struct_1d)
    #struct_3d_dict = read_regions_file(regions_struct_3d)

    gene_tumour_pairs_1d = set(gene_1d_dict.keys())
    gene_tumour_pairs_3d = set(gene_3d_dict.keys())

    gene_tumour_unique_1d = gene_tumour_pairs_1d - gene_tumour_pairs_3d
    gene_tumour_unique_3d = gene_tumour_pairs_3d - gene_tumour_pairs_1d
    gene_tumour_both = gene_tumour_pairs_3d & gene_tumour_pairs_1d


    # if the output directory doesn't exist, create it
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # create files that contain examples of genes
    # present in 1d and not 3d and vice versa
    ttype_ind = 0
    gene_ind = 1
    output = open(output_dir + "/genes_1d_not_3d_{0}.txt".format(num_codons), 'w')

    gene_unique_1d = [pair[gene_ind] for pair in gene_tumour_unique_1d]
    gene_unique_3d = [pair[gene_ind] for pair in gene_tumour_unique_3d]
    gene_unique_1d = set(sorted(gene_unique_1d))
    gene_unique_3d = set(sorted(gene_unique_3d))

    for gene in gene_unique_1d:
        output.write(gene +'\n')

    output.close()

    output = open(output_dir + "/genes_3d_not_1d_{0}.txt".format(num_codons), 'w')
    for gene in gene_unique_3d:
        output.write(gene +'\n')

    output.close()

    # create files that contain examples of gene, tumour pairs
    # present in 1d and not 3d and vice versa
    ttype_ind = 0
    gene_ind = 1
    output = open(output_dir + "/tumour_genes_1d_not_3d_{0}.txt".format(num_codons), 'w')

    for pair in gene_tumour_unique_1d:
        output.write(pair[ttype_ind] + '\t' + pair[gene_ind] +'\n')

    output.close()

    output = open(output_dir + "/tumour_genes_3d_not_1d_{0}.txt".format(num_codons), 'w')
    for pair in gene_tumour_unique_3d:
        output.write(pair[ttype_ind] + '\t' + pair[gene_ind] +'\n')

    output.close()

    output = open(output_dir + "/tumour_genes_all_{0}.txt".format(num_codons), 'w')
    for pair in gene_tumour_unique_3d:
        output.write(pair[ttype_ind] + '\t' + pair[gene_ind] + '\t3D' +'\n')
    for pair in gene_tumour_unique_1d:
        output.write(pair[ttype_ind] + '\t' + pair[gene_ind] + '\t1D' +'\n')
    for pair in gene_tumour_both:
        output.write(pair[ttype_ind] + '\t' + pair[gene_ind] + '\tboth' +'\n')
    output.close()


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
