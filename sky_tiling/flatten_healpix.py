import ligo.skymap
import ligo.skymap.moc
import ligo.skymap.io.fits
import healpy as hp
import numpy as np
import os

# Function to flatten multi-order HEALPix
# LIGO HEALPix maps are transmitted in a multi-order hierarchical structure,
# allowing pixels with high probability to be further subdivided.
# The tiling software we use does not support this format, so we flatten the structure.
# This may use more space (low-probability pixels must also be subdivided),
# so we also reduce the resolution
def flatten(input_fits_file, output_fits_file) :
    map, header = ligo.skymap.io.fits.read_sky_map(input_fits_file, nest=True)

    # First reduce NSIDE (resolution)

    target_nside = 32 # Pixels = 6 * NSIDE ^ 2

    # Compute scaling factor

    nside_sqr = len(map)/12
    nside = np.sqrt(nside_sqr).astype(int)
    assert nside * nside == nside_sqr, f"Imported map is not a valid nside!"

    scale = nside / target_nside
    scale = (scale * scale).astype(int)

    level = np.log2(scale).astype(int)
    assert 2**level == scale, f"nside {nside} note a power of 2 multiple of target_nside {target_nside}"

    print(f"Scaling down {scale}X")

    # Apply scaling factor
    map_reduced = np.array(map).reshape(-1,scale).sum(axis=-1)

    # Convert from NEST to RING pixel order
    # (This just defines how the HEALPix is represented in the file)
    map_ring = hp.pixelfunc.reorder(map_reduced, inp="NESTED", out="RING", n2r=True)

    # Write flattened file
    hp.write_map(output_fits_file, np.array(map_ring), overwrite=True)

# Main function, loop over HEALPix (.fits) files in directory,
# flatten and place in new directory 
if __name__ == "__main__":

    input_directory = "../data/ligo_healpix"
    output_directory = "../data/ligo_healpix_flattened"

    for filename in os.listdir(input_directory) :
        input_fits_file = os.path.join(input_directory, filename)
        output_fits_file = os.path.join(output_directory, filename)
        flatten(input_fits_file, output_fits_file)