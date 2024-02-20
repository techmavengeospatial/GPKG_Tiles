import gzip

filepath = r'D:\TMG\tasks\sample_data\MVT\tile_data.pbf'
with gzip.open(filepath, 'r') as fh:
    try:
        print(fh.read(1))
        print("gzip")
    except gzip.BadGzipFile:
        print('input_file is not a valid gzip file by BadGzipFile')
