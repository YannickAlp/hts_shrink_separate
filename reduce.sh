# run DBBAD
# hts_shrink_dbbad --chunk ../DEL_dataset_simpleAT/chunks/chunk_0.AP --active ../DEL_dataset_simpleAT/active_molecules.AP -np 12 -d 0.53

#!/bin/bash

# Spécifiez le chemin vers le dossier contenant les chunks
chemin_vers_chunks="del_dataset/data/chunks"

# Spécifiez le chemin vers le fichier contenant les molécules actives
chemin_vers_actives="del_dataset/data/active_molecules.AP"

# Nombre de processus parallèles
nombre_de_processus=12

# Distance
distance=0.53

# Boucle sur tous les fichiers du dossier "chunks"
for chunk in "$chemin_vers_chunks"/*; do
    if [ -f "$chunk" ]; then
        hts_shrink_dbbad --chunk "$chunk" --active "$chemin_vers_actives" -np "$nombre_de_processus" -d "$distance"
    fi
done