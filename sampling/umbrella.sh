#!/bin/bash

SAMPLING_PATH = 

cd $PATH

WINDOWS_close=$(seq 0 3 15)

mkdir wham

for d in $WINDOWS_close; do

	actual_distance=$( awk "BEGIN{ print $d / 100; exit}" )
	
	mkdir {$d}_close
	cd {$d}_close
	
	cp ../../../../SAM_sep/minim.mdp .
	cp ../../../../SAM_sep/md_pull_closecloseclose.mdp .
	sed -i s/StartingDistance/${actual_distance}/g NPT.mdp
	sed -i s/stepTime/${timeMove}/g NPT.mdp
	sed -i s/StartingDistance/${actual_distance}/g md_pull_closecloseclose.mdp


	gmx_mpi grompp -f md_pull_closecloseclose.mdp -c ../nvt_equil.gro -p ../../../HDT_AMD_025_B.top -r ../nvt_equil.gro -n ../index_index.ndx -o umbrella${d}_1.tpr  -maxwarn 3
	gmx_mpi mdrun -v -deffnm umbrella${d}_1 


	cp umbrella${d}_1_pullf.xvg ../wham/.
	cp umbrella${d}_1.tpr ../wham/.

	cd ..
done


