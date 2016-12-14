"""This script loads significant clusters into the MuPIT cluster tables."""
import MySQLdb
import os
import sys
import argparse

def parse_arguments():
    info = "This script loads significant clusters into the MuPIT cluster tables."
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-c', '--cluster',
                        type=str, required=True,
                        help='tab-delimited text file in the format of MuPIT\'s '
                        'cluster table')
    parser.add_argument('-r', '--residues',
                        type=str, required=True,
                        help='tab-delimited text file in the format of MuPIT\'s '
                        'cluster_residues table')
    parser.add_argument('-u', '--update-table',
                        action='store_true',
                        help='Prevent droping entire table before loading new data. '
                        'Only deletes data has the same tumor type as the input file.')
    parser.add_argument('--host',
                        type=str, default='karchin-db01',
                        help='MySQL host (default: karchin-db01)')
    parser.add_argument('--db',
                        type=str, default='mupit_modbase',
                        help='MySQL MuPIT database name (default: mupit_modbase)')
    parser.add_argument('--mysql-user',
                        type=str, required=True,
                        help='MySQL user name')
    parser.add_argument('--mysql-passwd',
                        type=str, required=True,
                        help='MySQL password')
    args = parser.parse_args()
    return vars(args)


def make_cluster_tables(cursor, db, cluster_file_name, residue_file_name):
    """Completely drop mutations table and re-load from file."""
    raise NotImplementedError  # not implemented yet
    cursor.execute('drop table if exists mutations')
    cursor.execute('create table mutations (source varchar(10), tissue varchar(10), structure_id varchar(20), residues varchar(10), occurrence int) engine=innodb')
    cursor.execute('load data local infile \'' + file_name + '\' into table mutations')
    db.commit()


def update_cluster_tables(cursor, db, cluster_file_name, residue_file_name):
    """Update the MuPIT mutations table with new data."""
    # find which tumor types are being added
    with open(cluster_file_name) as handle:
        uniq_ttypes = set(l.strip().split('\t')[3] for l in handle)

    # delete clusters that occur in the same tumor types
    # as those being loaded
    cursor.execute("SET @delTypes = '{0}'".format(','.join(uniq_ttypes)))
    delete_res = (
        "DELETE FROM cluster_residues "
        "WHERE cluster_residues.cluster_id IN ( "
            "SELECT * FROM ( "
                "SELECT cluster_id "
                "FROM cluster "
                "WHERE FIND_IN_SET(cluster.tissue, @delTypes)>0 "
            ") AS TMP "
        ")"
    )
    cursor.execute(delete_res)
    delete_clust = (
        "DELETE FROM cluster  "
        "WHERE cluster_id IN ( "
            "SELECT * FROM ( "
                "SELECT cluster_id "
                "FROM cluster c "
                "WHERE FIND_IN_SET(c.tissue, @delTypes)>0 "
            ") AS TMP "
        ")"
    )
    cursor.execute(delete_clust)

    # load input file now that we know there won't be duplicate mutations
    cursor.execute('load data local infile \'' + cluster_file_name + '\' into table cluster')
    cursor.execute('load data local infile \'' + residue_file_name + '\' into table cluster_residues')
    db.commit()


def main(opts):
    # input file to load
    cluster_path = opts['cluster']
    residue_path = opts['residues']

    # make mysql connection
    db = MySQLdb.connect(host=opts['host'],
                         user=opts['mysql_user'],
                         passwd=opts['mysql_passwd'],
                         db=opts['db'])
    cursor = db.cursor()

    # update mutations table
    if opts['update_table']:
        # just delete mutations that have the same tumor types as the
        # input file (out dated mutation counts).
        update_cluster_tables(cursor, db, cluster_path, residue_path)
    else:
        # drop table and completely re-load
        make_cluster_tables(cursor, db, cluster_path, residue_path)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
