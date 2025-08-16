import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import os

# Function to convert HEALPix likelihood maps to PNG images
def visualize(input_fits_file, output_png_file):
    # Read the HEALPix map from the FITS file
    skymap = hp.read_map(input_fits_file)

    # Create a Mollweide projection of the HEALPix map
    name = os.path.splitext(os.path.basename(input_fits_file))[0]
    hp.mollview(skymap, title=name, unit='Probability', cmap='viridis')

    # Save the visualization as a PNG file
    hp.graticule()
    plt.savefig(output_png_file, dpi=300)
    print(f"Saved visualization to {output_png_file}")
    plt.close()

# Main function, loop over HEALPix (.fits) files in directory,
# generate png image and place in new directory 
if __name__ == "__main__":

    input_directory = "../data/ligo_healpix_flattened"
    output_directory = "../data/ligo_healpix_images"

    for filename in os.listdir(input_directory) :
        input_fits_file = os.path.join(input_directory, filename)
        output_png_file = os.path.join(output_directory, filename.replace(".fits", ".png"))
        visualize(input_fits_file, output_png_file)