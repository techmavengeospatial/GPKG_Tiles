#!/usr/bin/python2.7
"""
Copyright (C) 2014 Reinventing Geospatial, Inc.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>,
or write to the Free Software Foundation, Inc., 59 Temple Place -
Suite 330, Boston, MA 02111-1307, USA.

Author: Steven D. Lander, Reinventing Geospatial Inc (RGi)
Date: 2013-07-12
   Requires: sqlite3, argparse
   Optional: Python Imaging Library (PIL or Pillow)
Description: Converts a TMS folder into a geopackage with
 PNGs for images with transparency and JPEGs for those
 without.
Credits:
  MapProxy imaging functions: http://mapproxy.org
  gdal2mb on github: https://github.com/developmentseed/gdal2mb

Version:
"""

from glob import glob
import gzip

# import buffer as buffer

from rgi.geopackage.common.zoom_metadata import ZoomMetadata
from rgi.geopackage.geopackage import Geopackage, PRAGMA_MINIMUM_SQLITE_VERSION
from rgi.geopackage.nsg_geopackage import NsgGeopackage
from rgi.geopackage.srs.ellipsoidal_mercator import EllipsoidalMercator
from rgi.geopackage.srs.geodetic import Geodetic
from rgi.geopackage.srs.geodetic_nsg import GeodeticNSG
from rgi.geopackage.srs.mercator import Mercator
from rgi.geopackage.srs.scaled_world_mercator import ScaledWorldMercator
from rgi.packaging.temp_db import TempDB
import tile2gpkg_vt

import json
try:
    from cStringIO import StringIO as ioBuffer
except ImportError:
    from io import BytesIO as ioBuffer
from time import sleep
from sys import stdout
from sys import version_info

if version_info[0] == 3:
    xrange = range

from sqlite3 import sqlite_version
from argparse import ArgumentParser
from sqlite3 import Binary as sbinary
from os import walk
from os.path import split, join, exists
from multiprocessing import cpu_count, Pool
# from distutils.version import LooseVersion
from looseversion import LooseVersion

try:
    from PIL.Image import open as IOPEN
except ImportError:
    IOPEN = None

# JPEGs @ 75% provide good quality images with low footprint, use as a default
# PNGs should be used sparingly (mixed mode) due to their high disk usage RGBA
# Options are mixed, jpeg, and png
IMAGE_TYPES = '.png', '.jpeg', '.jpg', '.pbf'


def write_geopackage_header(file_path):
    """
    writes geopackage header bytes to the sqlite database at file_path
    Args:
        file_path:

    Returns:
        nothing
    """
    header = 'GP10'
    with open(file_path, 'r+b') as file:
        file.seek(68, 0)
        file.write(header.encode())


def img_to_buf(img, img_type, jpeg_quality=75):
    """
    Returns a buffer array with image binary data for the input image.
    This code is based on logic implemented in MapProxy to convert PNG
    images to JPEG then return the buffer.

    Inputs:
    img -- an image on the filesystem to be converted to binary
    img_type -- the MIME type of the image (JPG, PNG)
    """
    # print("img_to_buf")
    defaults = {}
    buf = ioBuffer()
    if img_type == 'jpeg':
        img.convert('RGB')
        # Hardcoding a default compression of 75% for JPEGs
        defaults['quality'] = jpeg_quality
    elif img_type == 'source':
        img_type = img.format
    img.save(buf, img_type, **defaults)
    buf.seek(0)
    return buf


def img_has_transparency(img):
    """
    Returns a 0 if the input image has no transparency, 1 if it has some,
    and -1 if the image is fully transparent. Tiles *should be a perfect
    square (e.g, 256x256), so it can be safe to assume the first dimension
    will match the second.  This will ensure compatibility with different
    tile sizes other than 256x256.  This code is based on logic implemented
    in MapProxy to check for images that have transparency.

    Inputs:
    img -- an Image object from the PIL library
    """
    size = img.size[0]
    if img.mode == 'P':
        # For paletted images
        if img.info.get('transparency', False):
            return True
        # Convert to RGBA to check alpha
        img = img.convert('RGBA')
    if img.mode == 'RGBA':
        # Returns the number of pixels in this image that are transparent
        # Assuming a tile size of 256, 65536 would be fully transparent
        transparent_pixels = img.histogram()[-size]
        if transparent_pixels == 0:
            # No transparency
            return 0
        elif 0 < transparent_pixels < (size * size):
            # Image has some transparency
            return 1
        else:
            # Image is fully transparent, and can be discarded
            return -1
            # return img.histogram()[-size]
    return False


def file_count(base_dir):
    """
    A function that finds all image tiles in a base directory.  The base
    directory should be arranged in TMS format, i.e. z/x/y.

    Inputs:
    base_dir -- the name of the TMS folder containing tiles.

    Returns:
    A list of dictionary objects containing the full file path and TMS
    coordinates of the image tile.
    """
    print("Calculating number of tiles, this could take a while...")
    file_list = []
    # Avoiding dots (functional references) will increase performance of
    #  the loop because they will not be reevaluated each iteration.
    for root, sub_folders, files in walk(base_dir):
        temp_list = [join(root, f) for f in files if f.endswith(IMAGE_TYPES)]
        file_list += temp_list
    # print(file_list)
    print("Found {} total tiles.".format(len(file_list)))
    return [split_all(item) for item in file_list]


def split_all(path):
    """
    Function that parses TMS coordinates from a full images file path.

    Inputs:
    path -- a full file path to an image tile.

    Returns:
    A dictionary containing the TMS coordinates of the tile and its full
    file path.
    """
    parts = []
    full_path = path
    # Parse out the tms coordinates
    for i in xrange(3):
        head, tail = split(path)
        parts.append(tail)
        path = head
    file_dict = dict(y=int(parts[0].split('.')[0]),
                     x=int(parts[1]),
                     z=int(parts[2]),
                     path=full_path)
    return file_dict


def worker_map(temp_db, tile_dict, extra_args, invert_y):
    """
    Function responsible for sending the correct oriented tile data to a
    temporary sqlite3 database.

    Inputs:
    temp_db -- a temporary sqlite3 database that will hold this worker's tiles
    tile_dict -- a dictionary with TMS coordinates and file path for a tile
    tile_info -- a list of ZoomMetadata objects pre-generated for this tile set
    imagery -- the type of image format to send to the sqlite3 database
    invert_y -- a function that will flip the Y axis of the tile if present
    """
    tile_info = extra_args['tile_info']
    imagery = extra_args['imagery']
    jpeg_quality = extra_args['jpeg_quality']
    zoom = tile_dict['z']
    if extra_args['renumber']:
        zoom -= 1

    # print(tile_info)
    level = next((item for item in tile_info if item.zoom == zoom), None)
    # fiddle with offsets based on absolute (NSG profile) vs relative row/column numbering
    x_row = tile_dict['x'] if extra_args['nsg_profile'] else tile_dict['x'] - level.min_tile_row
    # print("y : " + str(tile_dict['y']))
    if invert_y is not None:
        y_column = invert_y(zoom, tile_dict['y'])
        # print("inverted y : " + str(y_column))
        if not extra_args['nsg_profile']:
            y_offset = invert_y(zoom, level.max_tile_col)
            y_column -= y_offset
    else:
        y_column = tile_dict['y'] if extra_args['nsg_profile'] else tile_dict['y'] - level.min_tile_col

    # print(imagery)
    if imagery == 'vt':
        # print(tile_dict['path'])
        file_handle = open(tile_dict['path'], 'rb')
        data = file_handle.read()
        # print(data)
        # Compress gzip
        compressed_data = gzip.compress(data)
        temp_db.insert_image_blob(zoom, x_row, y_column, compressed_data)
        file_handle.close()
    else:
        if IOPEN is not None:
            # print(tile_dict['path'])
            data = ioBuffer()
            img = IOPEN(tile_dict['path'], 'r')

            # TODO add options for "mvt" and "GeoJson"

            # print(imagery)
            if imagery == 'mixed':
                if img_has_transparency(img):
                    data = img_to_buf(img, 'png', jpeg_quality).read()
                else:
                    data = img_to_buf(img, 'jpeg', jpeg_quality).read()
            else:
                # print("other than mixed")
                data = img_to_buf(img, imagery, jpeg_quality).read()
            temp_db.insert_image_blob(zoom, x_row, y_column, sbinary(data))
        else:
            file_handle = open(tile_dict['path'], 'rb')
            data = buffer(file_handle.read())
            temp_db.insert_image_blob(zoom, x_row, y_column, data)
            file_handle.close()


def sqlite_worker(file_list, extra_args):
    """
    Worker function called by asynchronous processes.  This function
    iterates through a set of tiles to process them into a TempDB object.

    Inputs:
    file_list -- an array containing a subset of tiles that will be processed
                 by this function into a TempDB object
    base_dir -- the directory in which the geopackage will be created,
                .gpkg.part files will be generated here
    metadata -- a ZoomLevelMetadata object containing information about
                the tiles in the TMS directory
    """
    # TODO create the tempDB by adding the table name and telling which type (tiles/vectortiles)
    temp_db = TempDB(extra_args['root_dir'], extra_args['table_name'])
    with TempDB(extra_args['root_dir'],  extra_args['table_name']) as temp_db:
        invert_y = None
        if extra_args['lower_left']:
            if extra_args['srs'] == 3857:
                invert_y = Mercator.invert_y
            elif extra_args['srs'] == 4326:
                if extra_args['nsg_profile']:
                    invert_y = GeodeticNSG.invert_y
                else:
                    invert_y = Geodetic.invert_y
            elif extra_args['srs'] == 3395:
                invert_y = EllipsoidalMercator.invert_y
            elif extra_args['srs'] == 9804:
                invert_y = ScaledWorldMercator.invert_y
                #TODO update for retile
        # print(invert_y)
        [worker_map(temp_db, item, extra_args, invert_y) for item in file_list]


def allocate(cores, pool, file_list, extra_args):
    """
    Recursive function that fairly distributes tiles to asynchronous worker
    processes.  For N processes and C cores, N=C if C is divisible by 2.  If
    not, then N is the largest factor of 8 that is still less than C.
    """
    # if cores is 1:
    if cores == 1:
        print("Spawning worker with {} files".format(len(file_list)))
        return [pool.apply_async(sqlite_worker, [file_list, extra_args])]
    else:
        files = len(file_list)
        head = allocate(
            int(cores / 2), pool, file_list[:int(files / 2)], extra_args)
        tail = allocate(
            int(cores / 2), pool, file_list[int(files / 2):], extra_args)
        return head + tail


def build_lut(file_list, lower_left, srs):
    """
    Build a lookup table that aids in metadata generation.

    Inputs:
    file_list -- the file_list dict made with file_count()
    lower_left -- bool indicating tile grid numbering scheme (tms or wmts)
    srs -- the spatial reference system of the tile grid

    Returns:
    An array of ZoomLevelMetadata objects that describe each zoom level of the
    tile grid.
    """
    # Initialize a projection class
    if srs == 3857:
        projection = Mercator()
    elif srs == 4326:
        projection = Geodetic()
    elif srs == 9804:
        projection = ScaledWorldMercator()
    else:
        projection = EllipsoidalMercator()
    # Create a list of zoom levels from the base directory
    zoom_levels = list(set([int(item['z']) for item in file_list]))
    zoom_levels.sort()
    matrix = []
    # For every zoom in the list...
    for zoom in zoom_levels:
        # create a new ZoomMetadata object...
        level = ZoomMetadata()
        level.zoom = zoom
        # Sometimes, tiling programs do not generate the folders responsible
        # for the X axis if no tiles are being made within them.  This results
        # in tiles "shifting" because of the way they are renumbered when
        # placed into a geopackage.
        # To fix, is there a zoom level preceding this one...
        if zoom - 1 in [item for item in zoom_levels if item == (zoom - 1)]:
            # there is, now retrieve it....
            (prev,) = ([item for item in matrix if item.zoom == (zoom - 1)])
            # and fix the grid alignment values
            level.min_tile_row = 2 * prev.min_tile_row
            level.min_tile_col = 2 * prev.min_tile_col
            level.max_tile_row = 2 * prev.max_tile_row + 1
            level.max_tile_col = 2 * prev.max_tile_col + 1
            # Calculate the width and height
            level.matrix_width = prev.matrix_width * 2
            level.matrix_height = prev.matrix_height * 2
        else:
            # Get all possible x and y values...
            x_vals = [int(item['x'])
                      for item in file_list if int(item['z']) == zoom]
            y_vals = [int(item['y'])
                      for item in file_list if int(item['z']) == zoom]
            # then get the min/max values for each.
            level.min_tile_row, level.max_tile_row = min(x_vals), max(x_vals)
            level.min_tile_col, level.max_tile_col = min(y_vals), max(y_vals)
            # Fill in the matrix width and height for this top level
            x_width_max = max([item[
                                   'x'] for item in file_list if item['z'] == level.zoom])
            x_width_min = min([item[
                                   'x'] for item in file_list if item['z'] == level.zoom])
            level.matrix_width = (x_width_max - x_width_min) + 1
            y_height_max = max([item[
                                    'y'] for item in file_list if item['z'] == level.zoom])
            y_height_min = min([item[
                                    'y'] for item in file_list if item['z'] == level.zoom])
            level.matrix_height = (y_height_max - y_height_min) + 1
        level.min_x, level.min_y, level.max_x, level.max_y = calculate_top_left(level, projection, lower_left)
        # Finally, add this ZoomMetadata object to the list
        matrix.append(level)
    return matrix


def calculate_top_left(level, projection, lower_left):
    if lower_left:
        # TMS-style tile grid, so to calc the top left corner of the grid,
        # you must get the min x (row) value and the max y (col) value + 1.
        # You are adding 1 to the y value because the math to calc the
        # coord assumes you want the bottom left corner, not the top left.
        # Similarly, to get the bottom right corner, add 1 to x value.
        min_x, max_y = projection.get_coord(
            level.zoom, level.min_tile_row, level.max_tile_col + 1)
        max_x, min_y = projection.get_coord(
            level.zoom, level.max_tile_row + 1, level.min_tile_col)
    else:
        # WMTS-style tile grid, so to calc the top left corner of the grid,
        # you must get the min x (row value and the min y (col) value + 1.
        # You are adding 1 to the y value because the math to calc the
        # coord assumes you want the bottom left corner, not the top left.
        # Similarly, to get the bottom right corner, add 1 to x value.
        # -- Since this is WMTS, we must invert the Y axis before we calc
        inv_min_y = projection.invert_y(level.zoom, level.min_tile_col)
        inv_max_y = projection.invert_y(level.zoom, level.max_tile_col)
        min_x, max_y = projection.get_coord(
            level.zoom, level.min_tile_row, inv_min_y + 1)
        max_x, min_y = projection.get_coord(
            level.zoom, level.max_tile_row + 1, inv_max_y)
    return min_x, min_y, max_x, max_y


def build_lut_nsg(file_list, lower_left, srs, renumber):
    """
    Build a lookup table that aids in metadata generation, alternate method for NSG profile
    packages, which do not use relative coordinates.

    Inputs:
    file_list -- the file_list dict made with file_count()
    lower_left -- bool indicating tile grid numbering scheme (tms or wmts)
    srs -- the spatial reference system of the tile grid

    Returns:
    An array of ZoomLevelMetadata objects that describe each zoom level of the
    tile grid.
    """
    # Initialize a projection class
    if srs != 4326:
        print("NSG Profile requires that -srs be set to 4326")
        exit(1)
    # Currently, NSG profile support is only provided for epsg:4326.
    projection = GeodeticNSG()
    # Create a list of zoom levels from the base directory
    zoom_levels = list(set([int(item['z']) for item in file_list]))
    zoom_levels.sort()
    # If renumbering tiles we cannot have the old zoom level 0, so we remove it completely.
    if renumber:
        zoom_levels = [z for z in zoom_levels if z != 0]
    matrix = []
    # For every zoom in the list...
    for zoom in zoom_levels:
        # create a new ZoomMetadata object...
        level = ZoomMetadata()
        level.zoom = zoom
        if renumber:
            level.zoom -= 1

        # NSG profile geopackages do not use relative coordinates, so much of the
        # math involved in calculating conversions is no longer needed. However we do need
        # to update tile matrix calculations and still keep the min and max tile lists for use in the
        # contents table bounding box. we use the actual zoom rather than the (possibly) renumbered zoom
        # to  make sure we grab the correct tile locations
        # Get all possible x and y values...
        x_vals = [int(item['x'])
                  for item in file_list if int(item['z']) == zoom]
        y_vals = [int(item['y'])
                  for item in file_list if int(item['z']) == zoom]
        # Fill in the matrix width and height for this top level
        x_width_max = max(x_vals)
        x_width_min = min(x_vals)
        level.matrix_width = (x_width_max - x_width_min) + 1
        y_height_max = max(y_vals)
        y_height_min = min(y_vals)
        level.matrix_height = (y_height_max - y_height_min) + 1
        # then get the min/max available tiles for each. - for use in metadata
        level.min_tile_row, level.max_tile_row = min(x_vals), max(x_vals)
        level.min_tile_col, level.max_tile_col = min(y_vals), max(y_vals)
        # Fill in the matrix width and height for this top level
        # Because of tiling differences, we need to set the min and max based on the tiling format (TMS vs WMTS)
        level.min_x, level.min_y, level.max_x, level.max_y = calculate_top_left(level, projection, lower_left)
        #  use the renumbered zoom now to assign appropriate values.
        level.matrix_width, level.matrix_height = projection.get_matrix_size(level.zoom)
        level.min_tile_row, level.max_tile_row = 0, 2 ** (level.zoom + 1)
        level.min_tile_col, level.max_tile_col = 0, 2 ** level.zoom
        # Finally, add this ZoomMetadata object to the list
        matrix.append(level)
    return matrix


def combine_worker_dbs(out_geopackage):
    """
    Searches for .gpkg.part files in the base directory and merges them
    into one Geopackage file

    Inputs:
    out_geopackage -- the final output geopackage file
    """
    base_dir = split(out_geopackage.file_path)[0]
    if base_dir == "":
        base_dir = "."
    glob_path = join(base_dir + '/*.gpkg.part')
    file_list = glob(glob_path)
    print("Merging temporary databases...")
    # [out_geopackage.assimilate(f) for f in file_list]
    itr = len(file_list)
    status = ["|", "/", "-", "\\"]
    counter = 0
    for tdb in file_list:
        comp = len(file_list) - itr
        itr -= 1
        out_geopackage.assimilate(tdb)
        if tdb == file_list[-1]:
            stdout.write("\r[X] Progress: [" + "==" * comp + "  " * itr + "]")
        else:
            stdout.write("\r[" + status[counter] + "] Progress: [" + "==" *
                         comp + "  " * itr + "]")
        stdout.flush()
        if counter != len(status) - 1:
            counter += 1
        else:
            counter = 0
    print(" All geopackages merged!")


def main(arg_list):
    """
    Create a geopackage from a directory of tiles arranged in TMS or WMTS
    format.

    Inputs:
    arg_list -- an ArgumentParser object containing command-line options and
    flags
    """
    # TODO add argument for vector-tile format under imagery options
    # TODO add optional argument for "tiles" table name

    # Build the file dictionary
    files = file_count(arg_list.source_folder)
    if len(files) == 0:
        # If there are no files, exit the script
        print(" Ensure the correct source tile directory was specified.")
        exit(1)
    # Is the input tile grid aligned to lower-left or not?
    lower_left = arg_list.tileorigin == 'll' or arg_list.tileorigin == 'sw'
    # Get the output file destination directory
    root_dir, _ = split(arg_list.output_file)
    # Build the tile matrix info object
    if arg_list.nsg_profile:
        tile_info = build_lut_nsg(files, lower_left, arg_list.srs, arg_list.renumber)
    else:
        # print(files)
        tile_info = build_lut(files, lower_left, arg_list.srs)
        # print(tile_info[0].min_x)
        # print(tile_info[0].min_y)

    # Initialize the output file
    if arg_list.threading:
        # Enable tiling on multiple CPU cores
        cores = cpu_count()
        pool = Pool(cores)
        # Build allocate dictionary
        extra_args = dict(root_dir=root_dir,
                          tile_info=tile_info,
                          lower_left=lower_left,
                          srs=arg_list.srs,
                          imagery=arg_list.imagery,
                          jpeg_quality=arg_list.q,
                          nsg_profile=arg_list.nsg_profile,
                          renumber=arg_list.renumber,
                          table_name=arg_list.table_name)
        results = allocate(cores, pool, files, extra_args)
        status = ["|", "/", "-", "\\"]
        counter = 0
        try:
            while True:
                rem = sum([1 for item in results if not item.ready()])
                if rem == 0:
                    stdout.write("\r[X] Progress: [" + "==" * (cores - rem) +
                                 "  " * rem + "]")
                    stdout.flush()
                    print(" All Done!")
                    break
                else:
                    stdout.write("\r[" + status[counter] + "] Progress: [" +
                                 "==" * (cores - rem) + "  " * rem + "]")
                    stdout.flush()
                    if counter != len(status) - 1:
                        counter += 1
                    else:
                        counter = 0
                sleep(.25)
            pool.close()
            pool.join()
        except KeyboardInterrupt:
            print(" Interrupted!")
            pool.terminate()
            exit(1)
    else:
        # Debugging call to bypass multiprocessing (-T)
        extra_args = dict(root_dir=root_dir,
                          tile_info=tile_info,
                          lower_left=lower_left,
                          srs=arg_list.srs,
                          imagery=arg_list.imagery,
                          jpeg_quality=arg_list.q,
                          nsg_profile=arg_list.nsg_profile,
                          renumber=arg_list.renumber,
                          table_name=arg_list.table_name)
        sqlite_worker(files, extra_args)
    # Combine the individual temp databases into the output file
    # TODO GeoPackage and NSGGeoPackage need to be re-written to add Tiles or Vector tiles specifically
    if not arg_list.nsg_profile:
        with Geopackage(arg_list.output_file, arg_list.srs, arg_list.table_name) as gpkg:
            gpkg.initialize()
            combine_worker_dbs(gpkg)
            # Using the data in the output file, create the metadata for it
            gpkg.update_metadata(tile_info)

    else:
        with NsgGeopackage(arg_list.output_file, arg_list.srs, arg_list.table_name) as gpkg:
            gpkg.initialize()
            combine_worker_dbs(gpkg)
            # Using the data in the output file, create the metadata for it
            # print("else")
            gpkg.update_metadata(tile_info)

    # we do a late write of the applicaiton id if its needed, to allow time for the database  connections to clear out
    if LooseVersion(sqlite_version) < LooseVersion(PRAGMA_MINIMUM_SQLITE_VERSION):
        write_geopackage_header(arg_list.output_file)

    if arg_list.imagery == 'vt':
        out = tile2gpkg_vt.modifieGPKGFORVT(arg_list)
        if out == 'updated':
            print("Complete")
        else:
            print(out)
    else:
        print("Complete")


if __name__ == '__main__':
    print("""
        tiles2gpkg_parallel.py  Copyright (C) 2014  Reinventing Geospatial, Inc
        This program comes with ABSOLUTELY NO WARRANTY.
        This is free software, and you are welcome to redistribute it
        under certain conditions.
    """)
    PARSER = ArgumentParser(description="Convert TMS folder into geopackage")
    PARSER.add_argument("source_folder",
                        metavar="source",
                        help="Source folder of TMS files.")
    PARSER.add_argument("output_file",
                        metavar="dest",
                        help="Destination file path.")
    PARSER.add_argument("-tileorigin",
                        metavar="tile_origin",
                        help="Tile point of origin location. Valid options " +
                             "are ll, ul, nw, or sw.",
                        choices=["ll", "ul", "sw", "nw"],
                        default="ll")
    PARSER.add_argument("-srs",
                        metavar="srs",
                        help="Spatial reference " + "system. Valid options are"
                             + "3857, 4326, 3395, and 9804.",
                        type=int,
                        choices=[3857, 4326, 3395, 9804],
                        default=3857)

    # TODO: to support vector tiles, expand the choices to include "MVT" and "GeoJSON"
    PARSER.add_argument("-imagery",
                        metavar="imagery",
                        help="Imagery type. Valid options are mixed, " +
                             "jpeg, png, , source or vt.",
                        choices=["mixed", "jpeg", "png", "source", "vt"],
                        default="source")

    PARSER.add_argument("-table_name",
                        metavar="table_name",
                        help="The name of the tiles table.",
                        default="tiles")

    PARSER.add_argument("-q",
                        metavar="quality",
                        type=int,
                        default=75,
                        help="Quality for jpeg images, 0-100. Default is 75",
                        choices=list(range(100)))
    PARSER.add_argument("-a",
                        metavar="append",
                        # dest="append",
                        # action="store_true",
                        default=False,
                        help="Append tile set to existing geopackage")
    PARSER.add_argument("-T",
                        dest="threading",
                        action="store_false",
                        default=True,
                        help="Disable multiprocessing.")
    PARSER.add_argument("-renumber",
                        dest="renumber",
                        action="store_true",
                        default=False,
                        help="Enable re-numbering tiles/zoom levels from standard Geodetic to NSG geodetic. Only valid"
                             "if the NSG Profile is enabled.")
    group = PARSER.add_mutually_exclusive_group(required=False)
    group.add_argument("-nsg",
                       dest="nsg_profile",
                       help="Enforce NSG Profile Requirements on output GeoPackage. Currently Requires data to"
                            "use the EPSG:4326 Global Geodetic projection. Note: it will not convert tiles to the"
                            "proper Tile Matrix, they MUST tiled correctly to be packaged correctly.",
                       action='store_true',
                       )
    group.add_argument("-ogc",
                       dest="nsg_profile",
                       help="Follow OGC GeoPackage specification without NSG Profile additions",
                       action='store_false')
    PARSER.set_defaults(nsg_profile=False)

    ARG_LIST = PARSER.parse_args()
    if not ARG_LIST.a:
        if not exists(ARG_LIST.source_folder) or exists(ARG_LIST.output_file):
            PARSER.print_usage()
            print("Ensure that TMS directory exists and out file does not.")
            exit(1)
    if ARG_LIST.q is not None and ARG_LIST.imagery == 'png':
        PARSER.print_usage()
        print("-q cannot be used with png")
        exit(1)
    if ARG_LIST.nsg_profile and ARG_LIST.srs != 4326:
        PARSER.print_usage()
        print("-nsg requires that -srs be set to 4326")
        exit(1)
    if not ARG_LIST.nsg_profile and ARG_LIST.renumber:
        PARSER.print_usage()
        print("-renumber requires that the -nsg flag also be active")

    main(ARG_LIST)
