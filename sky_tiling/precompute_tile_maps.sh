#!/bin/bash

python setup.py --telescope 7DT_recomputed --fov 2.5 --work ./tile_pixel_maps --site 7DT --timemag None --extension png
mv ./tile_pixel_maps/config.ini ./tile_pixel_maps/config_7DT_recomputed.ini
python setup.py --telescope 7DT_separated_recomputed --fov 10 --work ./tile_pixel_maps --site 7DT --timemag None --extension png 
mv ./tile_pixel_maps/config.ini ./tile_pixel_maps/config_7DT_separated_recomputed.ini