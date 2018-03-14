#!/bin/bash

python dem_density.py China.xyz > dem_China_density.txt

tessmodgen -s0.016667/0.016667 -z0 -v < dem_China_density.txt \
> dem-tess.txt