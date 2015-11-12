import csv
import argparse

def parse_arguments():
    info = (
        'Finds the necessary mutation count threshold for significance in '
        'each structure.'
    )
    parser = argparse.ArgumentParser(description=info)

    # program arguments
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Hotspot output from hotspot.py script')
    parser.add_argument('-alpha', '--alpha',
                        type=str, required=True,
                        help='file containing tumor type specific alpha levels.')
    parser.add_argument('-o', '--output',
                        default='output.txt',
                        type=str,
                        help='Output result file of hotspots')
    args = parser.parse_args()
    opts = vars(args)
    return opts


def main(opts):
    # read in significance level thresholds
    with open(opts['alpha']) as handle:
        data = [l.strip().split('\t') for l in handle]
        thresholds = {d[0]: float(d[1]) for d in data}

    with open(opts['input']) as handle:
        myreader = csv.reader(handle, delimiter='\t')

        # define column positions
        header = next(myreader)
        dens_ix = header.index('Mutation Density')
        pval_ix = header.index('Hotspot P-value')
        ttype_ix = header.index('Tumor Type')
        struct_ix = header.index('Structure')

        # iterate over every structure/tumor type pair
        output = [['Structure', 'Tumor Type', 'Count Threshold']]
        for line in myreader:
            # additional info
            ttype = line[ttype_ix]
            struct = line[struct_ix]

            # no p-value case, this can happen because of our filtering
            # out of heteroatoms from PDB structures
            if not line[pval_ix]:
                output.append([struct, ttype, 'No Mutations'])
                continue

            # extract mutation counts and associated p-values
            pvals = map(float, line[pval_ix].split(','))
            dens = map(int, line[dens_ix].split(','))

            # make sure nothing is significant if no residue
            # is found with a p-value equal/below cutoff
            default_cutoff = max(dens) + 1000

            if ttype in thresholds:
                # find if there is any significant counts
                signif_ix = [i for i in range(len(pvals))
                             if pvals[i]<=thresholds[ttype]]

                # set cutoff
                if signif_ix:
                    cutoff = min(dens[i] for i in signif_ix)
                else:
                    cutoff = default_cutoff  # no res signif
            else:
                # this is the case where no residues are significant
                # for this particular tumor type
                cutoff = default_cutoff

            output.append([struct, ttype, cutoff])

    # write output
    with open(opts['output'], 'wb') as out_file:
        mywriter = csv.writer(out_file, delimiter='\t', lineterminator='\n')
        mywriter.writerows(output)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
