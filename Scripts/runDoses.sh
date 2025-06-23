#!/bin/bash

# List of e_fermi values to loop through
e_fermi_values=(-40 300)  # You can change or expand this list

# Backup the original input.nml file
cp ./input.nml ./input.nml.bak

for val in "${e_fermi_values[@]}"; do
    echo "Setting e_fermi = $val"

    # Copy the backup to restore original before each change
    cp ./input.nml.bak ./input.nml

    # Replace the e_fermi line with the new value
    sed -i "s/^\s*e_fermi\s*=.*/  e_fermi = $val/" ./input.nml

    # You can call your program here using the modified input.nml
    ./bin/POST_LAO_STO.x
    cp ./OutputData/SuperconductingGap.dat ./OutputData/SuperconductingGap_${val}.dat

done
