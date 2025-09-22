#!/bin/bash

source $HOME/local_installs/gromacs_plumed/bin/GMXRC

SAMPLING_PATH= 
INITIAL_CONFIGURATION_PATH= 

START=
INTERVAL=
END=

cd $PATH

WINDOWS_close=$(seq START INTERVAL END)

mkdir wham

for d in $WINDOWS_close; do

	actual_distance=$( awk "BEGIN{ print $d / 100; exit}" )
	
	mkdir {$d}_close
	cd {$d}_close
	
	cp ../../../../SAM_sep/minim.mdp .
	cp ../../../../SAM_sep/md_pull.mdp .
	sed -i s/StartingDistance/${actual_distance}/g md_pull.mdp


	gmx_mpi grompp -f md_pull_close.mdp -c $INITIAL_CONFIGURATION_PATH/nvt_equil.gro -p $INITIAL_CONFIGURATION_PATH/HDT_AMD_025_B.top -r $INITIAL_CONFIGURATION_PATH/nvt_equil.gro -n $INITIAL_CONFIGURATION_PATH/index_index.ndx -o umbrella${d}_1.tpr  -maxwarn 3
	gmx_mpi mdrun -v -deffnm umbrella${d}_1 


	cp umbrella${d}_1_pullf.xvg ../wham/.
	cp umbrella${d}_1.tpr ../wham/.

	cd ..
done


