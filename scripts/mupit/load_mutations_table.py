"""This script loads mutation counts into the MuPIT mutations table."""
import MySQLdb
import os
import sys
import argparse

def parse_arguments():
    info = "This script loads mutation counts into the MuPIT table called mutations."
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-m', '--mutations',
                        type=str, required=True,
                        help='tab-delimited text file in the format of MuPIT\'s '
                        'mutations table')
    parser.add_argument('-u', '--update-table',
                        action='store_true',
                        help='Prevent droping entire table before loading new data. '
                        'Only deletes data has the same tumor type as the input file.')
    parser.add_argument('--host',
                        type=str, required=True,
                        help='MySQL host')
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


def make_mutations_table(cursor, db, file_name):
    """Completely drop mutations table and re-load from file."""
    cursor.execute('drop table if exists mutations')
    cursor.execute('create table mutations (source varchar(10), tissue varchar(10), structure_id varchar(20), residues varchar(10), occurrence int) engine=innodb')
    cursor.execute('load data local infile \'' + file_name + '\' into table mutations')
    db.commit()


def update_mutations_table(cursor, db, file_name):
    """Update the MuPIT mutations table with new data."""
    # find which tumor types are being added
    with open(file_name) as handle:
        uniq_ttypes = set(l.split('\t')[1] for l in handle)

    # delete mutations that occur in the same tumor types
    # as those being loaded
    cursor.execute("SET @delTypes = '{0}'".format(','.join(uniq_ttypes)))
    cursor.execute("DELETE FROM mutations "
                   "WHERE FIND_IN_SET(mutations.tissue, @delTypes)>0")

    # load input file now that we know there won't be duplicate mutations
    cursor.execute('load data local infile \'' + file_name + '\' into table mutations')
    db.commit()


def main(opts):
    # input file to load
    data_path = opts['mutations']

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
        update_mutations_table(cursor, db, data_path)
    else:
        # drop table and completely re-load
        make_mutations_table(cursor, db, data_path)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
