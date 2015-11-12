import argparse
import csv
import re
import IPython


# import modules needed for logging
import logging
import os
import random

logger = logging.getLogger(__name__)  # module logger

def parse_arguments():
    """
    Function to parse command line arguements
    from the user

    Returns
    -------
    opts : dict
        command line arguements from the user
    """

    info = 'Divides pdb info files for parallelization'
    parser = argparse.ArgumentParser(description=info)

    # program arguments
    parser.add_argument('-f', '--in-file',
                        type=str,
                        required=True,
                        help='PDB info file to divide')
    parser.add_argument('-n', '--num-splits',
                        default=1000,
                        type=int,
                        help='Number of splits to perform (Default: 1000)')
    parser.add_argument('-m', '--mut-file',
                        type=str,
                        required=True,
                        help='File containing mutation information')
    parser.add_argument('--split-dir',
                        default = "../data/split_pdbs/",
                        type=str,
                        help='Output directory for split PDB info files')

    args = parser.parse_args()
    opts = vars(args)
    return opts


def read_file(in_file):
    '''
    Read the file and output into a dictionary
    '''

    dictionary = {}
    # iterate through file
    for line in in_file:
        # tab delimited file
        fields = line.split('\t')
        # leave out the ones with no info
        if fields[1] == '':
            continue
        # check if id is already in the dictionary
        if not fields[0] in dictionary:
            dictionary[fields[0]] = []
        dictionary[fields[0]].append(line)
    return dictionary


def split_file(in_file, num_splits, split_dir, mut_file):
    """
    Splits the file and creates files in the desired directory
    """

    # create the output directory if it does
    # not exist
    if not os.path.exists(split_dir):
        os.mkdir(split_dir)

    # open the info file
    f = open(in_file)
    pdb_header = f.readline()

    # open the mutation file
    m = open(mut_file)
    mut_header = m.readline()

    # read into a dictionary containing
    # structure ids as keys and lines pertaining
    # to it as values
    pdb_dict = read_file(f)
    mut_dict = read_file(m)

    # determine total num of ids in file
    total_ids = len(pdb_dict.keys())
    print(total_ids)
    # determine num of ids to put in each split
    num_ids = int(total_ids/num_splits)

    # counters
    count_file = 0
    count_id = num_ids

    # randomize order of insertions
    keys = pdb_dict.keys()
    random.shuffle(keys)

    # iterate through dict and write to files
    #for key in sorted(pdb_dict):
    for key in keys:

        # check if we need a new file
        if (count_id == num_ids and count_file < num_splits):
            count_id = 0
            pdb_out = open(split_dir + "/pdb_info_split_" + str(count_file) + ".txt", 'w')
            pdb_out.write(pdb_header)
            mut_out = open(split_dir + "/mut_info_split_" + str(count_file) + ".txt", 'w')
            mut_out.write(mut_header)
            count_file += 1

        # write all lines pertaining to the structure id
        for line in pdb_dict[key]:
            pdb_out.write(line)
        if key in mut_dict:
            for line in mut_dict[key]:
                mut_out.write(line)

        count_id += 1


def main(opts):
    """
    Splits a PDB file into the desired number of splits

    Parameters
    ----------
    opts: dict
        command line arguements from the user
    """

    # split the file
    split_file(opts['in_file'], opts['num_splits'], opts['split_dir'], opts['mut_file'])


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)

