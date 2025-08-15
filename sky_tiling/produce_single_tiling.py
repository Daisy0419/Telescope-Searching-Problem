name = "GW200216_220804"
telescope = "7dt"
confidence_interval = 0.99
max_rows = None

import numpy as np
import bin.rankedTilesGenerator as rankedTilesGenerator


# Creates a tiling using Shaon's method:
# https://github.com/shaonghosh/sky_tiling/
# Note that we modify the bin/rankedTilesGenerator.py with a new method, generateTable
def create_tiling(skymap_file, config_file, output_csv_file, confidence_interval, containment, max_rows=None) :
    tileObj = rankedTilesGenerator.RankedTileGenerator(skymap_file, config_file)
    [ranked_tile_indices, ranked_tile_probs] = tileObj.getRankedTiles(resolution=256)
    tbl = tileObj.generateTable(ranked_tile_indices, ranked_tile_probs, CI=confidence_interval, containment=containment, max_rows=max_rows)
    tbl.write(output_csv_file, format="csv", overwrite=True)

# Invoke create_tiling function for each file in list
def generate_tiles_from_maps(file_list, which, confidence_interval, containment, max_rows = None) :

    for idx, file in enumerate(file_list) :

        # Split to extract FITS file name and telescope name
        parts = file.split("_", 2)
        fits = parts[0] + "_" + parts[1]
        telescope = parts[2]

        skymap_file = "../data/ligo_healpix_flattened/" + fits + ".fits"
        config_file = "tile_pixel_maps/config_" + telescope + ".ini"
        output_file = "../data/"+ which + "/filtered_" + fits + "_" + telescope + ".csv"

        print "Processing", which, "HEALPix map", idx, "of", len(file_list), ":", output_file

        create_tiling(skymap_file, config_file, output_file, confidence_interval, containment, max_rows)


if __name__ == "__main__":

    name = "GW200216_220804"
    telescope = "7dt"

    skymap_file = "../"

    # Large instances: 99% containment
    generate_tiles_from_maps(large_files, "large", confidence_interval=0.99, containment=False)

    # Small instances: up to 99.5% containment or up to 100 tiles (smaller of two)
    generate_tiles_from_maps(small_files, "small", confidence_interval=0.995, containment=True, max_rows=100)
