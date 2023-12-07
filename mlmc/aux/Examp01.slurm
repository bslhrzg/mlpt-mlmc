#!/bin/bash

#SBATCH -p devtaras 
#SBATCH -N 1
#SBATCH --ntasks=4
#SBATCH --job-name=Examp01
#SBATCH --output=output-%x-%N.%j
#SBATCH --error=error-%x-%N.%j
#SBATCH --time=00:10:00

  # number of MC steps
  mctot=25

  # number of back report
  # it means that each 10 cycle 
  # the total average energy will be printed 
  # in the slurm output ! A good number should 
  # be each 1000 step do a report and copy if you 
  # are running mctot=50000-200000
  nwrite=5

  # temperature of the MC simulation
  temperature="300.00"

  # The vasp exec should be changed
  # if your not running a gamma point calculation
  # and also the vasp version 
  module load vasp/5.4.4
  export vasp=$(which vasp_gam)

########################################################
# you don't need to change anything from here; however 
# if you need to change a little the script, do it by 
# your risk. 
########################################################

  cd $SLURM_SUBMIT_DIR

  ulimit -s unlimited

  export OMP_NUM_THREADS=1
  export MKL_NUM_THREADS=1

  export TMPDIR=$(mktemp -d -p /tmp)

  echo '---- Began at ----'
  date
  echo ''

  echo '---- Cluster/Job Info ----'
  echo "    Job ID                     $SLURM_JOB_ID"
  echo "    Job Name                   $SLURM_JOB_NAME"
  echo "    Submit Directory           $SLURM_SUBMIT_DIR"
  echo "    Temporary Directory        $TMPDIR" 
  echo "    Host submitted from        $SLURM_SUBMIT_HOST"
  echo "    Number of nodes allocated  $SLURM_JOB_NUM_NODES"
  echo "    Number of cores/node       $SLURM_CPUS_ON_NODE"
  echo "    Total number of cores      $SLURM_NTASKS"
  echo "    Nodes assigned to job      $SLURM_JOB_NODELIST"
  echo ''

  cd $TMPDIR

  FOUT="$TMPDIR/MC-runall.txt"
  FXYZ="$TMPDIR/ACCEPTED.xyz"

  cp $SLURM_SUBMIT_DIR/aux/* .
  cp $SLURM_SUBMIT_DIR/runcalc/* .

##first iter
  iter=1

  my_training.py

  cp POSCAR.0 POSCAR
  mpirun -n $SLURM_NTASKS $vasp > vasp.out 2>&1 

  e0_prd=$( grep e_wo_entrp vasprun.xml | tail -1 | awk '{print $3}' )
  e0_tar=$( my_prediction.py $e0_prd )
  acc=1

  update_geom_accepted.py

  printf "%7d%12.6f%12.6f%3d\n" $iter $e0_prd $e0_tar $acc > $FOUT

  echo ""
  echo "  First Iteration:"
  echo ""
  printf "%7d%12.6f%12.6f%3d\n" $iter $e0_prd $e0_tar $acc 
  echo ""

  echo ""
  echo "  Reporting Averaged Target Energies at $nwrite steps:"
  echo ""

## mc iteraton
  while [ $iter -le $mctot ]; do

    mc_generator.py

    cp POSCAR POSCAR_
    cp hCONTCAR POSCAR

    mpirun -n $SLURM_NTASKS $vasp > vasp.out 2>&1 

    e1_prd=$( grep e_wo_entrp vasprun.xml | tail -1 | awk '{print $3}' )
    e1_tar=$( my_prediction.py $e1_prd )

    acc=$( detector.py $e0_tar $e1_tar $temperature )
    if [ $acc -eq 1 ]; then 
      mv hCONTCAR POSCAR
      e0_prd=$( echo $e1_prd )
      e0_tar=$( echo $e1_tar )
      update_geom_accepted.py
    else
      cp POSCAR_ POSCAR
    fi

    printf "%7d%12.6f%12.6f%3d\n" $iter $e0_prd $e0_tar $acc >> $FOUT

    nprt=$( echo $iter | awk -v n="$nwrite" '{print $1%n}' )
    if [ $nprt -eq 0 ] ; then 
      cp $FOUT $SLURM_SUBMIT_DIR
      cp $FXYZ $SLURM_SUBMIT_DIR
      avg=$( awk '{s+=$3}END{printf "%20.10e\n",s/NR}' ./MC-runall.txt )
      std=$( awk -v m="$avg" '{s+=($3-m)**2}END{printf "%20.10e\n",sqrt(s/(NR*(NR-1)))}' ./MC-runall.txt )
      printf "\t\tIter:%10d\t\tAvg:%20.10e\t\tStd:%20.10e\n" $iter $avg $std
    fi 

    iter=$((iter+1))

  done 

## copying output to home
  mv $FOUT $SLURM_SUBMIT_DIR
  mv $FXYZ $SLURM_SUBMIT_DIR

  cd $SLURM_SUBMIT_DIR

## cleaning the TMPDIR 
  rm $TMPDIR/__pycache__/*
  rmdir $TMPDIR/__pycache__
  rm $TMPDIR/*
  rmdir $TMPDIR

  echo ""
  echo '---- Finished at ----'
  date
  echo ''
