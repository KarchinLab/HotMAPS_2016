import os
from Bio.PDB.parse_pdb_header import parse_pdb_header
import argparse
import csv
import re
import gzip

def parse_arguments():
    info = 'Add chain description to PDB information file'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Input PDB info file without chain descriptions')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output PDB info file WITH chain description')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in PDB info file
    with open(opts['input']) as handle:
        myreader = csv.reader(handle, delimiter='\t')
        header = next(myreader)
        header.append('ChainDescription')
        header2ix = {h: i for i, h in enumerate(header)}

        # found vs missing lines
        found_lines = []
        missing_lines = []
        missing_lines_no_pdb = []

        # iterate over each line
        pdb_id = ''
        for line in myreader:
            # figure out if the PDB file exists
            non_bio_pdb_path = line[header2ix['NonBioPDBPath']]
            is_empty =  non_bio_pdb_path == ''
            is_valid_path = os.path.exists(non_bio_pdb_path)

            # skip over lines without a good path
            if is_empty or not is_valid_path:
                if line[header2ix['PDBPath']]:
                    # Likely case of homology model
                    line.append('')
                    missing_lines.append(line)
                else:
                    # no PDB file found
                    line.append('')
                    missing_lines_no_pdb.append(line)
                continue

            # catch parsing errors in Bio.PDB
            try:
                # parse the header if you have seen the PDB before
                if line[header2ix['PDBId']] != pdb_id:
                    # read in pdb
                    pdb_id = line[header2ix['PDBId']]
                    if non_bio_pdb_path.endswith('.gz'):
                        # handle compressed case
                        with gzip.open(non_bio_pdb_path, 'rb') as gzhandle:
                            pdb_header = parse_pdb_header(gzhandle)
                    else:
                        # normal .ent file
                        pdb_header = parse_pdb_header(non_bio_pdb_path)

                    # process pdb header from COMPND lines
                    chain2desc = {}
                    compound = pdb_header['compound']
                    for key in compound:
                        # get the list of chain ids that have the same description
                        mychains = map(lambda s: s.strip(),
                                       re.split('\s*,\s*', compound[key]['chain']))
                        # PDB files with missing chain IDs are assigned as 'A'
                        mychains = [m.upper() if m else 'A' for m in mychains]
                        # get the chain description from the molecule field
                        chain_desc = compound[key]['molecule'].strip().replace(' ', '_')

                        # iterate over each chain so they are assigned a chaind desc
                        for c in mychains:
                            chain2desc[c] = chain_desc

                # add chain description to line
                chain_letter = line[header2ix['chain']]
                tmp_chain_desc = chain2desc[chain_letter]
                line.append(tmp_chain_desc)
                found_lines.append(line)
            except KeyboardInterrupt:
                raise
            except:
                print('Skipping structure: {0}'.format(pdb_id))

    # sort result
    pdb_ix = header2ix['PDBId']
    desc_ix = header2ix['ChainDescription']
    chain_ix = header2ix['chain']
    found_lines = sorted(found_lines, key=lambda x: (x[pdb_ix], x[desc_ix], x[chain_ix]))
    missing_lines = sorted(missing_lines, key=lambda x: (x[pdb_ix], x[desc_ix], x[chain_ix]))
    missing_lines_no_pdb = sorted(missing_lines_no_pdb, key=lambda x: (x[pdb_ix], x[desc_ix], x[chain_ix]))

    # write output
    with open(opts['output'], 'w') as write_handle:
        mywriter = csv.writer(write_handle, delimiter='\t', lineterminator='\n')
        out = [header] + found_lines + missing_lines + missing_lines_no_pdb
        mywriter.writerows(out)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

