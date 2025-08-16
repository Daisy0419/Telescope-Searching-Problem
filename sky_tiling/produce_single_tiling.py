name = "GW200216_220804" # Name of the GW event for which you would like to make tiles
telescope = "7dt" # Replace this with the name you selected for the telescope
confidence_interval = 0.99 # Include highest-probability tiles until
                           # cumulative probability exceeds this value
max_rows = None # Maximum number of tiles to include, None means no limit

import numpy as np
import os
import bin.rankedTilesGenerator as rankedTilesGenerator


# Creates a tiling using Shaon's method:
# https://github.com/shaonghosh/sky_tiling/
# Note that we modify the bin/rankedTilesGenerator.py with a new method, generateTable
def create_tiling(skymap_file, config_file, output_csv_file, confidence_interval, containment, max_rows=None) :
    tileObj = rankedTilesGenerator.RankedTileGenerator(skymap_file, config_file)
    [ranked_tile_indices, ranked_tile_probs] = tileObj.getRankedTiles(resolution=256)
    tbl = tileObj.generateTable(ranked_tile_indices, ranked_tile_probs, CI=confidence_interval, containment=containment, max_rows=max_rows)
    tbl.write(output_csv_file, format="csv", overwrite=True)


if __name__ == "__main__":


    # Create directory for output files
    skymap_dir = "../data/custom"
    if not os.path.exists(skymap_dir):
        os.makedirs(skymap_dir)

    skymap_file = "../data/ligo_healpix_flattened/" + name + ".fits"
    config_file = "tile_pixel_maps/config_" + telescope + ".ini"
    output_file = skymap_dir + "/" + name + "_" + telescope + ".csv"

    # Check if config_file exists, if not, create and write default content
    if not os.path.exists(config_file):
        default_config_content = """
[pixelTileMap]
preComputed_64 = ./tile_pixel_maps/preComputed_%s_pixel_indices_64.dat
preComputed_128 = ./tile_pixel_maps/preComputed_%s_pixel_indices_128.dat
preComputed_256 = ./tile_pixel_maps/preComputed_%s_pixel_indices_256.dat
preComputed_512 = ./tile_pixel_maps/preComputed_%s_pixel_indices_512.dat
preComputed_1024 = ./tile_pixel_maps/preComputed_%s_pixel_indices_1024.dat
preComputed_2048 = ./tile_pixel_maps/preComputed_%s_pixel_indices_2048.dat

[tileFiles]
tileFile = ./tile_center_files/%s_tiles_indexed.dat

[plot]
filenametag = %s
extension = png

[observation]
site = %s
time_magnitude = None
trigger_time = 

        """ % (telescope, telescope, telescope, telescope, telescope, telescope, telescope, telescope, telescope)
        with open(config_file, "w") as f:
            f.write(default_config_content)

    # Create the tiling for the custom skymap
    create_tiling(skymap_file, config_file, output_file, confidence_interval, containment=False, max_rows=max_rows)
    print "Created tiling for skymap:", skymap_file
    print "Output saved to:", output_file