#! /bin/bash

module load PE-gnu/2.0

$ROSETTA_BIN/extract_pdbs.default.linuxgccrelease -in:file:silent $1 -in:file:tags $2 -overwrite
