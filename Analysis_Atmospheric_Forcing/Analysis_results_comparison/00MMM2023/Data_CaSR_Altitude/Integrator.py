#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2024-12-24

@author: Nicolas Gasset
"""

import xarray as xr
import argparse

# Parametres 
parser = argparse.ArgumentParser(description="Combiner plusieurs fichiers NetCDF.")
parser.add_argument("-i","--input_files", nargs="+", help="Fichiers NetCDF en entree (au moins un).", required=True)
parser.add_argument("-o","--output_file", help="Fichier NetCDF de sortie.", required=True)
args = parser.parse_args()

input_files=args.input_files
output_file=args.output_file

# Liste des fichiers
tile_files = sorted(input_files)

print(f"Liste des fichiers d'entree :")
print(tile_files)

# Ouvre et combine les datasets pas coordonnees
combined = xr.open_mfdataset(tile_files, combine='by_coords', parallel=True)

# Definie correctement 'rotated_pole' variable
# Ne depends pas de rlat/rlon
ds = xr.open_dataset(tile_files[0])
combined['rotated_pole']=ds['rotated_pole']

# Sauve le dataset recombine dans un nouveau fichier NetCDF avec la dimension time UNLIMITED
combined.to_netcdf(output_file,unlimited_dims='time')

print(f"Fichier combine sauvegarde dans : {output_file}")
