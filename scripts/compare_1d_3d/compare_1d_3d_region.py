import argparse
import csv
import re
import time
import IPython

# import modules needed for logging
import logging
import os
import random

logger = logging.getLogger(__name__)  # module logger
DROPPED_TTYPES = ["COAD", "READ", "THYM", "CHOL", "UVM", "TGCT", "SARC", "PANCAN12"]

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
    parser.add_argument('-r', '--reference-region',
                        type=str, required=True,
                        help='Reference set of gene-level hotspot regions')
    parser.add_argument('-c', '--compare-region',
                        type=str, required=True,
                        help='Set of gene-level hotspot regions to compare to '
                        'the reference set')
    #parser.add_argument('-d', '--regions-dir',
                        #type=str, required=True,
                        #help='Directory containing the regions .txt files for 1D and 3D')
    #parser.add_argument('-r', '--radius-1d',
                        #type=str, required=True,
                        #help='The number of codons both upstream and downstream used for 1D algorithm')
    parser.add_argument('-q', '--q-value',
                        type=float, default=.01,
                        help='Significance level of q-value in txt files to analyse')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='File to write the output to')
    args = parser.parse_args()
    opts = vars(args)
    return opts


def main(opts):
    """
    Main function
    """
    # get specified parameters
    #regions_dir = opts['regions_dir']

    # get all the required filenames
    regions_gene_ref = opts['reference_region']
    regions_gene_compare = opts['compare_region']
    #regions_gene_1d = regions_dir + "/hotspot_regions_gene_." + q_val + "_1d_{0}.txt".format(num_codons)
    #regions_gene_3d = regions_dir + "/hotspot_regions_gene_." + q_val + ".txt"

    # read in hotspot region file for 3D
    with open(regions_gene_ref) as handle:
        myreader = csv.reader(handle, delimiter='\t')
        hotspots_ref = []
        for line in myreader:
            if line[1] in DROPPED_TTYPES:
                continue
            else:
                hotspots_ref.append(line)

    # read in hotspot region file for 1D
    with open(regions_gene_compare) as handle:
        myreader = csv.reader(handle, delimiter='\t')
        hotspots_compare = list(myreader)

    # identify all unique genes
    uniq_genes = set([h[0] for h in hotspots_ref])

    # check which genes have difference in their location for hotspots
    output = [['Gene', 'Tumor Type (reference)', 'Hotspot ID',
               'Hotspot Region Residues', 'Tumor Type (compare)', 'Overlap']]
    for gene in uniq_genes:
        # get hotspots for a single gene
        hotspots_ref_gene = [h for h in hotspots_ref if h[0]==gene]
        hotspots_compare_gene = [h for h in hotspots_compare if h[0]==gene]

        # get clusters for reference
        clust_list = []
        ttype_list = []
        for myhotspot in hotspots_ref_gene:
            clusts = [set(c.split(';')) for c in myhotspot[2:]]
            clust_list.append(clusts)
            ttype_list.append(myhotspot[1])

        # get clusters for compared regions
        clust_list_compare = []
        ttype_list_compare = []
        for myhotspot in hotspots_compare_gene:
            clusts = [set(c.split(';')) for c in myhotspot[2:]]
            clust_list_compare.append(clusts)
            ttype_list_compare.append(myhotspot[1])

        # compare clusters in different tumor types to see if there
        # are differences
        for i in range(len(clust_list)):
            found_match_ttype = False
            for j in range(len(clust_list_compare)):
                # skip if not matching tumor type
                if ttype_list[i] != ttype_list_compare[j]:
                    continue

                # set flag indicating found a matching tumor type
                # for this gene
                found_match_ttype = True

                # get num of hotspot regions
                num_i = len(clust_list[i])
                num_j = len(clust_list_compare[j])

                # iterate over each region
                for l in range(num_i):
                    intersect_all = 0
                    for m in range(num_j):
                        intersect = len(clust_list[i][l] & clust_list_compare[j][m])
                        intersect_all += intersect
                    res_str = ';'.join(clust_list[i][l])
                    out_line = [gene, ttype_list[i], l, res_str,
                                ttype_list_compare[j], intersect_all]
                    output.append(out_line)

            # catch case where compared list does not have any hotspots
            # in the same gene
            if not found_match_ttype:
                for l in range(len(clust_list[i])):
                    intersect_all = 0
                    res_str = ';'.join(clust_list[i][l])
                    out_line = [gene, ttype_list[i], l, res_str,
                                ttype_list[i], intersect_all]
                    output.append(out_line)

    # write output
    with open(opts['output'], 'w') as handle:
        mywriter = csv.writer(handle, delimiter='\t', lineterminator='\n')
        mywriter.writerows(output)


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
