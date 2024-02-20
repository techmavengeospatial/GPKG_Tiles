import sqlite3
from sqlite3 import Error
import os
import json


def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)

    return conn


def updateGPKGContent(conn, table_name):
    """
    :param conn:
    :param table_name:
    """
    res = ''
    try:
        sql = "UPDATE gpkg_contents set data_type='vector-tiles',identifier='" + table_name + "',description='TechMoven' where table_name = '" + table_name + "'"
        cur = conn.cursor()
        cur.execute(sql)
        conn.commit()
        # conn.close()
        res = 'updated'
    except Error as e:
        conn.close()
        res = e
        # print(e)

    return res


def modifieGPKGFORVT(arg_list):
    database = arg_list.output_file
    source_folder = arg_list.source_folder
    # create a database connection
    conn = create_connection(database)
    # with conn:
    res = updateGPKGContent(conn, arg_list.table_name)
    create_necessary_tabls(conn, source_folder, arg_list.table_name)
    return res


def create_necessary_tabls(conn, source_folder, table_name):
    cur = conn.cursor()
    create_geopackage_vt_fields(cur)
    create_geopackage_vt_layers(cur)
    # conn.close()
    if os.path.isfile(source_folder):
        print(f"{source_folder} exists.")
    else:
        source_folder = os.path.join(source_folder, 'metadata.json')
        # print(f"{source_folder} does not exist.")

    if os.path.exists(source_folder):
        f = open(source_folder)
        # returns JSON object as
        data = json.load(f)
        if data is not None:
            r = json.loads(data['json'])
            for r in r['vector_layers']:
                f = r['fields']
                cols = list(f.keys())
                vals = list(f.values())
                cur.execute(
                    'INSERT INTO gpkgext_vt_layers (table_name,name,description,minzoom,maxzoom,attributes_table_name) VALUES (?,?,?,?,?,?)',
                    (table_name, r['id'], r['description'], r['minzoom'], r['maxzoom'], table_name))
                seq_id = cur.lastrowid
                i = 0
                for fld in cols:
                    cur.execute(
                        'INSERT INTO gpkgext_vt_fields (layer_id,name,type) VALUES (?,?,?)',
                        (seq_id, fld, vals[i]))
                    i = i+1

            conn.commit()
            conn.close()
    else:
        print(f"{source_folder} does not exist.")


def create_geopackage_vt_fields(cursor):
    cursor.execute("""
                     CREATE TABLE IF NOT EXISTS {table_name}
                     (
                     id INTEGER PRIMARY KEY AUTOINCREMENT,
                      layer_id  int NOT NULL,
                      name  TEXT,
                      type TEXT
                      )
                   """.format(table_name='gpkgext_vt_fields'))


def create_geopackage_vt_layers(cursor):
    cursor.execute("""
                     CREATE TABLE IF NOT EXISTS {table_name}
                     (
                     id INTEGER PRIMARY KEY AUTOINCREMENT,
                      table_name  TEXT NOT NULL,
                      name  TEXT,
                      description TEXT,
                      minzoom int,
                      maxzoom int,
                      attributes_table_name TEXT,
                      geometry_type_name TEXT
                      )
                   """.format(table_name='gpkgext_vt_layers'))
