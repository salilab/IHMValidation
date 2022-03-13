#!/bin/bash

# Simple script to get all .cif files from PDB-Dev and make the
# validation report for each one

MAX_ID=100  # Update as more structures are added to PDB-Dev

for id in $(seq 1 ${MAX_ID} ); do
  pdbid=$(printf "PDBDEV_%08d.cif" ${id})
  curl -fsLO https://pdb-dev.wwpdb.org/cif/${pdbid}
  if [ $? -eq 0 ]; then
    python Execute.py -f ${pdbid} >& ${id}.log
    if [ $? -ne 0 ]; then
      echo "${pdbid} FAILED; see ${id}.log"
    else
      echo "${pdbid} done"
    fi
  else
    echo "${pdbid} could not be downloaded"
  fi
done
