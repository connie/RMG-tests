#!/bin/bash

# This is a manual checking sript for the minimal target job.

target=$1
DIFF_DIR=/home/connieg/Diff
RMGPy=/home/connieg/Code/RMG-Py
TESTS_DIR=/home/connieg/Code/RMG-tests
RMGDatabase=/home/connieg/Code/RMG-database





echo 'Location of RMG-Py home directory: '$RMGPy
echo 'Target: ' $target

cd $RMGPy >/dev/null
echo 'RMG-Py Git Branch: ' git branch

echo 'Latest RMG-Py Git Commit: ' git log -n 1


echo 'Location of the RMG-database: ' $RMGDatabase

cd $RMGDatabase >/dev/null
echo 'RMG-database Git Branch: ' git branch
echo 'Latest RMG-database Git Commit: ' git commit



echo 'Making the target'
cd $RMGPy
make $target


SOURCE_FOLDER=$RMGPy/testing/$target/chemkin/


cd $DIFF_DIR
mkdir $target
cd $target

# check generated models:
# core:
python $TESTS_DIR/checkModels.py $target $SOURCE_FOLDER/chem_annotated.inp $SOURCE_FOLDER/species_dictionary.txt
# edge:
python $TESTS_DIR/checkModels.py $target $SOURCE_FOLDER/chem_edge_annotated.inp $SOURCE_FOLDER/species_edge_dictionary.txt