# TMG_GeoPackage-VectorTiles_python_automation

This solution builds OGC GPKG GeoPackage Map Tiles

## Supported Input Tiles

> **Raster Tiles PNG/JPEG**

> **Heightmap**

> **PNG Elevation-Terrain Tiles**

> **Vector Tiles**

It accepts MBTILES as your input and stores that into gpkg as seperate table/layer.

## Please see our Commercial Offerings (Windows Apps)

**Tile Utilities**

> https://portfolio.techmaven.net/apps/tile-utilities/

**Map Tiling**

> https://maptiling.techmaven.net/

## Geospatial Data Serving

> https://geodataserver.techmaven.net/

> https://tileserver.techmaven.net/

## Mobile Apps

**We also have mobile apps that work with GeoPackage**

> https://portfolio.techmaven.net/apps/

> https://mapexplorer.techmaven.net/

> https://earthexplorer.techmaven.net/

> https://mapdiscovery.techmaven.net/

> https://geonamesmapexplorer.techmaven.net/

> https://geodatacollector.techmaven.net

## Getting Started

1. Options to use the command line tool

   - -h, --help show this help message and exit

   - -i INPUT, --input=INPUT Input path of mbtiles

   - -o OUTPUT, --output=OUTPUT Output path of gpkg

   - -p PROJECTION, --proj=PROJECTION ie. 3395/3857/4326

   - -r RESOURCE, --resource=RESOURCE Directory path of resources folder

   - -t TABLE, --table=TABLE Table name in geopackage

## Example

`--proj possible values can be 3395 3957 4326`

**3395 Projection**

> python tiles_to_rbt_gpkg.py --input "input file path of mbtiles" --output "output gpkg path" --proj 3395 --resource "Resource folder directory path" --table "tbl_hillshades"

**3857 Projection**

> python tiles_to_rbt_gpkg.py --input "input file path of mbtiles" --output "output gpkg path" --proj 3857 --resource "Resource folder directory path" --table "tbl_cultural"

**4326 Projection**

> python tiles_to_rbt_gpkg.py --input "input file path of mbtiles" --output "output gpkg path" --proj 4326 --resource "Resource folder directory path" --table "tbl_physical"

## Dependencies

`optparse`

> pip install optparse

`pyproj`

> pip install pyproj

`sqlite3`

> pip instal sqlite3
