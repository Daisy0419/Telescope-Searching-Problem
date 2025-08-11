
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

# Lists of files for large and small problem instances

large_files = [
    "GW191103_012549_7dt",
    "GW191109_010717_7dt",
    "GW191113_071753_7dt",
    "GW191126_115259_7dt",
    "GW191127_050227_7dt",
    "GW191204_110529_7dt",
    "GW191216_213338_7dt",
    "GW191219_163120_7dt",
    "GW191222_033537_7dt",
    "GW191230_180458_7dt",
    "GW200105_162426_7dt",
    "GW200112_155838_7dt",
    "GW200128_022011_7dt",
    "GW200208_222617_7dt",
    "GW200209_085452_7dt",
    "GW200210_092254_7dt",
    "GW200216_220804_7dt",
    "GW200220_061928_7dt",
    "GW200220_124850_7dt",
    "GW200302_015811_7dt",
    "GW200306_093714_7dt",
    "GW200308_173609_7dt_separate",
    "GW200322_091133_7dt_separate",
    "GW230529_181500_7dt_separate",
]

small_files = [
    "GW191105_143521_7dt_separate",
    "GW191129_134029_7dt_separate",
    "GW191204_171526_7dt_separate",
    "GW191215_223052_7dt_separate",
    "GW191216_213338_7dt_separate",
    "GW191230_180458_7dt_separate",
    "GW200115_042309_7dt_separate",
    "GW200129_065458_7dt_separate",
    "GW200202_154313_7dt_separate",
    "GW200209_085452_7dt_separate",
    "GW200219_094415_7dt_separate",
    "GW200225_060421_7dt_separate",
    "GW200316_215756_7dt_separate"
]

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

    # Large instances: 99% containment
    generate_tiles_from_maps(large_files, "large", confidence_interval=0.99, containment=False)

    # Small instances: up to 99.5% containment or up to 100 tiles (smaller of two)
    generate_tiles_from_maps(small_files, "small", confidence_interval=0.995, containment=True, max_rows=100)
