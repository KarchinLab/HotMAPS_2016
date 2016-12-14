"""Makes a blacklist of residues that map to multiple postions within a single
PDB chain.
"""
import pandas as pd
import os
import argparse


def parse_arguments():
    info = ('This script makes a residue blacklist for mutations that maps to more than one '
            'position within a PDB protein chain')
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-a', '--annotation-dir',
                        type=str, required=True,
                        help='Annotation directory containing both the CRAVAT and '
                        'mupit mapping to the PDB structure')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--gene", action="store_true",
                       default = False,
                       help="Flag to identify black list residues based on the gene "
                       "reference transcript and codon position")
    group.add_argument("--structure", action="store_true",
                       default=False,
                       help="Flag to identify black list residues based on the "
                       "pdb structure positions")
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file containing multi-mapping mutations')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in all the mutation info
    useful_columns = ['pdb_id', 'chain', 'residue', 'ID', 'Sample ID',
                      'HUGO symbol', 'Chromosome', 'Position',
                      'Reference Transcript', 'Reference Codon Position']
    df_list = []
    for file_name in os.listdir(opts['annotation_dir']):
        file_path = os.path.join(opts['annotation_dir'], file_name)
        ttype_df = pd.read_csv(file_path, sep='\t', usecols=useful_columns)
        df_list.append(ttype_df)
    df = pd.concat(df_list)

    # remove bad residue numbers
    bad_res = df['residue'].str.contains('[A-Z]', na=False)
    df = df[~bad_res]

    # count the duplicates
    grp_cols = ['pdb_id', 'chain', 'ID', 'Sample ID', 'HUGO symbol',
                'Chromosome', 'Position',
                'Reference Transcript', 'Reference Codon Position'
                ]
    grp = df[grp_cols+['residue']].copy().groupby(grp_cols)
    cts = grp['residue'].apply(lambda x: len(x.unique()))
    repeats = cts[cts>1].reset_index()
    repeats = repeats.rename(columns={0: 'repeats'})

    # figure out which residues match those that can have multiple mappings
    tx_pos_cols = ['Reference Transcript', 'Reference Codon Position']
    dup_pos = set(map(tuple,
                      repeats[tx_pos_cols].values))
    cts_noindex = cts.reset_index().rename(columns={0: 'repeats'})
    is_dup_codon = cts_noindex[tx_pos_cols].apply(lambda x, dup_pos: tuple(x) in dup_pos,
                                                  args=(dup_pos,), axis=1)
    dups = cts_noindex[is_dup_codon]

    # figure out gene-level information
    if opts['gene']:
        # get the residue information at the gene-level
        res_info = dups['Reference Transcript'] + ":" + dups['Reference Codon Position'].astype(int).astype(str)
        dups['Reference Residue Info'] = res_info

        # save the blacklist to a file
        out_cols = ['HUGO symbol', 'Reference Residue Info']
        dups = dups.drop_duplicates(out_cols)
        dups[out_cols].to_csv(opts['output'], sep='\t', index=False)
    # figure out structure-level information
    else:
        # parse out which structure residues are duplicates
        # from the original dataframe object df
        def is_dup(row, dup_info):
            """Function to find if row is in duplicates."""
            row_tuple = tuple(row.tolist())
            return row_tuple in dup_info
        struct_cols = ['pdb_id', 'chain', 'ID']
        duplicate_info = set([tuple(l) for l in dups[struct_cols].values.tolist()])
        dup_row_flag = df[struct_cols].apply(is_dup, args=(duplicate_info,), axis=1)
        df = df[dup_row_flag]

        # drop duplicate res info
        df = df.drop_duplicates(['pdb_id', 'chain', 'residue'])

        # get the residue information at the gene-level
        res_info = df['chain'] + ":" + df['residue'].astype(int).astype(str)
        df['Structure Residue Info'] = res_info
        not_null = ~df['Reference Codon Position'].isnull()
        res_info = df.loc[not_null, 'Reference Transcript'] + ":" + df.loc[not_null, 'Reference Codon Position'].astype(int).astype(str)
        df['Reference Residue Info'] = ''
        df.loc[not_null, 'Reference Residue Info'] = res_info

        # save the blacklist to a file
        out_cols = ['HUGO symbol', 'pdb_id', 'Structure Residue Info',
                    'Reference Residue Info']
        df[out_cols].to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
