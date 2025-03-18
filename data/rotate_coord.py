import numpy as np
import healpy as hp
import pandas as pd

df = pd.read_csv("deep_slow.csv") 

def sph_to_cart(ra_deg, dec_deg):
    ra_rad = np.deg2rad(ra_deg)
    dec_rad = np.deg2rad(dec_deg)
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    return x, y, z

def cart_to_sph(x, y, z):
    dec_rad = np.arcsin(z)
    ra_rad = np.arctan2(y, x)
    ra_deg = np.rad2deg(ra_rad) % 360
    dec_deg = np.rad2deg(dec_rad)
    return ra_deg, dec_deg

def transform_coordinates(ra_deg, dec_deg, rotation_matrix):
    x, y, z = sph_to_cart(ra_deg, dec_deg)
    x_new, y_new, z_new = np.dot(rotation_matrix, [x, y, z])
    ra_rotated, dec_rotated = cart_to_sph(x_new, y_new, z_new)

    return ra_rotated, dec_rotated

# rotation_matrix = np.array([
#     [0, 0, 1],
#     [0, 1, 0],
#     [-1, 0, 0]
# ])

rotation_matrix = np.array([
    [1, 0, 0],
    [0, 0, 1],
    [0, -1, 0]
])

ra_transformed = []
dec_transformed = []
for ra, dec in zip(df["RA"], df["Dec"]):
    ra_t, dec_t = transform_coordinates(ra, dec, rotation_matrix)
    ra_transformed.append(ra_t)
    dec_transformed.append(dec_t)

df["RA"] = ra_transformed
df["Dec"] = dec_transformed


print("Transformed Data:")
print(df[["Rank", "Index", "RA", "Dec", "Probability"]])
df.to_csv("transformed_data.csv", index=False)
print("Data saved to 'transformed_data.csv'.")