#!/bin/bash

OUTPUT_DIR_DOS="./OutputData/DOS_A1/"
OUTPUT_DIR_GAP="./OutputData/Gap_A1/"

mkdir -p $OUTPUT_DIR_DOS
mkdir -p $OUTPUT_DIR_GAP
# List of e_fermi values to loop through
e_fermi_values=(50 100)  # You can change or expand this list

python3 Runner/OutputMocker.py

# Backup the original input.nml file
cp ./input.nml ./input.nml.bak

for val in "${e_fermi_values[@]}"; do
    echo "Setting e_fermi = $val"

    # Copy the backup to restore original before each change
    cp ./input.nml.bak ./input.nml

    # Replace the e_fermi line with the new value
    sed -i "s/^\s*e_fermi\s*=.*/  e_fermi = $val/" ./input.nml

    # You can call your program here using the modified input.nml
    srun -c 48 ./bin/POST_LAO_STO.x
    #cp ./OutputData/DOS.dat ${OUTPUT_DIR_DOS}/DOS_${val}.dat
    cp ./OutputData/SuperconductingGap.dat ${OUTPUT_DIR_GAP}/SuperconductingGap_${val}.dat

done
