import healpy as hp
import matplotlib.pyplot as plt

def visualize_converted_map(file):
    """Visualize a HEALPix map in Mollweide projection for paper use."""
    hmap = hp.read_map(file, verbose=False)

    # Configure figure
    hp.mollview(hmap,
                title='',  # remove title
                cmap='viridis',
                cbar=True,
                notext=True,
                norm='hist')  # histogram equalization for contrast

    fig = plt.gcf()
    fig.set_size_inches(6, 4)  # Adjust size for paper (can go larger for better resolution)

    # Optional: Customize colorbar
    cb = plt.gca().images[-1].colorbar
    if cb:
        cb.ax.tick_params(labelsize=8)
        cb.set_label("Posterior", fontsize=9)

    plt.savefig("ligo_map.png", dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

file="../../data_RTSS_update/GW200311_115853.fits_7dt.csv"
visualize_converted_map(file)