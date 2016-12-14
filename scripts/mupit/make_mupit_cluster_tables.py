"""Parses my hotspot regions for PDB structures and converts it
to a format that can be uploaded into the MUPIT database."""
import argparse
import csv


def parse_arguments():
    info = ('Convert hotspot regions for structures in to a format ready '
            'for input into the mupit mysql database')
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-r', '--region',
                        type=str, required=True,
                        help='Hotspot regions for structures')
    parser.add_argument('-b', '--blacklist',
                        type=str, default=None,
                        help='Blacklisted residues due to multi-mapping (optional).')
    parser.add_argument('-reg', '--region-table',
                        type=str, required=True,
                        help='Formated output text file describing hotspot regions')
    parser.add_argument('-res', '--residue-table',
                        type=str, required=True,
                        help='Formated output text file describing the '
                        ' residues within hotspot regions')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    """This function flattens out the text file output from
    find_hotspot_regions_struct.py into a format useable by mupit.
    """
    # read in blacklist, if provided
    if opts['blacklist']:
        with open(opts['blacklist']) as handle:
            myreader = csv.reader(handle, delimiter='\t')
            header = next(myreader)
            struct_ix = header.index('pdb_id')
            res_ix = header.index('Structure Residue Info')

            # add every pdb residue in file into blacklist
            blacklist = set()
            for row in myreader:
                res_info = (row[struct_ix], row[res_ix])
                blacklist.add(res_info)
    else:
        # no residues to blacklist
        blacklist = set()

    # read in hotspot regions found in PDB structures
    with open(opts['region']) as handle:
        myreader = csv.reader(handle, delimiter='\t')

        # iterate over each region
        source = 'tcga'
        region_output = []
        residue_output = []
        reg_num = 0
        for line in myreader:
            # get info about regions
            struct_id, tumor_type = line[:2]

            # add information about residues
            already_found_res = set()
            for myregion in line[2:]:
                # set object created to de-duplicate residues from
                # different bioassembly models
                res_info_set = set([tuple(res.split(':')[1:])
                                   for res in myregion.split(';')
                                   if (struct_id, res[2:]) not in blacklist])

                # eliminate duplicate reporting
                res_info_set_new =  res_info_set - already_found_res

                # only add region if its new
                if res_info_set_new:
                    # add residue info to output
                    res_info_list = [[reg_num]+list(r) for r in res_info_set]
                    residue_output.extend(res_info_list)

                    # add region info to output
                    reg_info = [reg_num, source, struct_id, tumor_type]
                    region_output.append(reg_info)

                    # update hotspot region unique ID
                    reg_num += 1

                    # add residues to already found set
                    already_found_res |= res_info_set

    # write output for region table
    with open(opts['region_table'], 'w') as region_handle:
        mywriter = csv.writer(region_handle, delimiter='\t', lineterminator='\n')
        mywriter.writerows(region_output)

    # write output for residue table
    with open(opts['residue_table'], 'w') as residue_handle:
        mywriter = csv.writer(residue_handle, delimiter='\t', lineterminator='\n')
        mywriter.writerows(residue_output)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
