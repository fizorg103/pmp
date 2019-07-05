#!/bin/bash
for file in input/structures/*.xyz; 
do
  python3 mol_graph/xyz2mol.py file
done
