#!/usr/bin/env python
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))

# package import
import src.utils as utils

import csv
import argparse
import logging

logger = logging.getLogger(__name__)  # module logger

# directories for various pdb files
import ConfigParser
config = ConfigParser.SafeConfigParser()
config.read(os.path.join(file_dir, '../config.txt'))
refseq_dir = config.get('PDB', 'refseq_homology')
ensembl_dir = config.get('PDB', 'ensembl_homology')
biounit_dir = config.get('PDB', 'biological_assembly')
pdb_dir = config.get('PDB', 'non_biological_assembly')

def parse_arguments():
    info = 'Adds column for path to the correct PDB file'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-p', '--pdb-info',
                        type=str, required=True,
                        help='PDB Info file from mupit_modbase')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Updated PDB info file with path to PDBs')

    # logging arguments
    parser.add_argument('-ll', '--log-level',
                        type=str,
                        action='store',
                        default='',
                        help='Write a log file (--log-level=DEBUG for debug mode, '
                        '--log-level=INFO for info mode)')
    parser.add_argument('-l', '--log',
                        type=str,
                        action='store',
                        default='',
                        help='Path to log file. (accepts "stdout")')
    args = parser.parse_args()

    # handle logging
    if args.log_level or args.log:
        if args.log:
            log_file = args.log
        else:
            log_file = ''  # auto-name the log file
    else:
        log_file = os.devnull
    log_level = args.log_level
    utils.start_logging(log_file=log_file,
                        log_level=log_level)  # start logging
    return vars(args)


def main(opts):
    output = []
    missing_files = []
    missing_non_bio_files = []
    count = 0
    with open(opts['pdb_info']) as handle:
        # handle header line
        header = handle.readline().strip().split('\t')
        header.append('NonBioPDBPath')
        header.append('PDBPath')
        output.append(header)

        # iterate over each line
        for line in csv.reader(handle, delimiter='\t'):
            struct_id, chain, gene = line
            count += 1

            # figure out where the pdb file is
            non_biounit_path = ''
            if struct_id.startswith('ENSP'):
                # pdb_path = os.path.abspath(os.path.join(ensemble_dir, '{0}.pdb'.format(struct_id[4:].lstrip('0'))))
                pdb_path = os.path.abspath(os.path.join(ensembl_dir, '{0}.pdb'.format(struct_id)))
            elif struct_id.startswith('NP'):
                pdb_path = os.path.abspath(os.path.join(refseq_dir, '{0}.pdb'.format(struct_id)))
            else:
                # try to get non-biological assembly pdb file path
                putative_path = os.path.join(pdb_dir, 'pdb{0}.ent.gz'.format(struct_id))
                # the second path may occur if user downloads PDB from RCSB
                putative_pat2 = os.path.join(pdb_dir, '{0}.pdb.gz'.format(struct_id))
                if os.path.exists(putative_path):
                    non_biounit_path = putative_path
                elif os.path.exists(putative_path2):
                    non_biounit_path = putative_path2
                else:
                    missing_non_bio_files.append(putative_path)

                # try version numbers from 1 to 32
                for i in range(33):
                    pdb_path = os.path.abspath(os.path.join(biounit_dir, '{0}.pdb{1}.gz'.format(struct_id, i)))
                    if os.path.exists(pdb_path):
                        break
                #pdb_path = os.path.abspath(os.path.join(pdb_dir, 'pdb{0}.ent'.format(struct_id)))

            line.append(non_biounit_path)
            if not os.path.exists(pdb_path):
                logger.debug('Could not find PDB: {0}'.format(pdb_path))
                # print(pdb_path)
                line.append('')
                missing_files.append(pdb_path)
            else:
                line.append(pdb_path)

            # add line to output
            output.append(line)

    # write output to file
    with open(opts['output'], 'w') as handle:
        csv_writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        csv_writer.writerows(output)

    print('There are {0} missing pdb files'.format(len(set(missing_files))))
    print('There are {0} missing non-bioassembly pdb file'.format(len(set(missing_non_bio_files))))
    print('NOTE: Several dozen missing PDBs are expected because the PDB database is constantly changing')



if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
