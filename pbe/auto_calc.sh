#!/bin/bash

surfaces=(100 110)

for i in "${surfaces[@]}"; do (
  ads_sites=("$i"/*)

  for j in "${ads_sites[@]}"; do (
    cd "$j"
    echo "----------------------------------------------------"
    echo "$j"

    poetry run ads_calcs.py reorder-atoms -f ./geometry.in -b Cu -c 6

    if [[ -d ./2_run_cont/ ]]; then
      cd ./2_run_cont/
    elif [[ -d ./run_cont/ ]]; then
      cd ./run_cont/
    fi

    poetry run ads_calcs.py reorder-atoms -f ./geometry.in.next_step -b Cu -c 6
    poetry run ads_calcs.py substrate-info -i ./geometry.in.next_step -b Cu -p -a C -l 8 | tee ./ads_info.txt
    echo
    poetry run ads_calcs.py parse-energies -f ./ -so /home/dylanmorgan/comp_chem/warwick/azupyrene_cu/refs/"$i"/ -ss /home/dylanmorgan/comp_chem/warwick/azupyrene_cu/refs/"$i"_sp/ -ao /home/dylanmorgan/comp_chem/warwick/azupyrene_cu/refs/azpy/ -as /home/dylanmorgan/comp_chem/warwick/azupyrene_cu/refs/azpy_sp/ -a C | tee ./energy_info.txt
    echo
    )
    done )
  done
