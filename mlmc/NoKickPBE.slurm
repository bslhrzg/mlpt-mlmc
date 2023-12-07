#!/bin/bash

#SBATCH -p tarasall
#SBATCH -N 1
#SBATCH --ntasks=16
#SBATCH --job-name=t-hchab
#SBATCH --output=output-%x-%N.%j
#SBATCH --error=error-%x-%N.%j

  # number of MC steps
  mctot=80000

  nwrite=10

  # temperature of the MC simulation
  temperature="300.00"


  module load vasp/5.4.4
  export vasp=$(which vasp_gam)

  cd $SLURM_SUBMIT_DIR

  ulimit -s unlimited

  export OMP_NUM_THREADS=1
  export MKL_NUM_THREADS=1

  export TMPDIR=$(mktemp -d -p /scratch/$USER)

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
  update_geom_accepted.py

  e0_prd=$( grep e_wo_entrp vasprun.xml | tail -1 | awk '{print $3}' )
  e0_tar=$( my_prediction.py $e0_prd )
  acc=1

  echo $iter $e0_prd $e0_tar $acc > $FOUT

  echo ""
  echo "  First Iteration:"
  echo ""
  echo $iter $e0_prd $e0_tar $acc 
  echo ""

  echo ""
  echo "  Reporting Averaged Target Energies at $nwrite steps:"
  echo ""

  iter=$((iter+1))

## mc iteraton
  while [ $iter -le $mctot ]; do

    ./mc_generator.py

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
      update_geom_refused.py
    fi

    echo $iter $e0_prd $e0_tar $acc >> $FOUT

    nprt=$( echo $iter | awk -v n="$nwrite" '{print $1%n}' )
    if [ $nprt -eq 0 ] ; then 
      cp $FOUT $SLURM_SUBMIT_DIR
      cp $FXYZ $SLURM_SUBMIT_DIR
      avg=$( awk '{s+=$3}END{print s/NR}' ./MC-runall.txt )
      std=$( awk -v m="$avg" '{s+=($3-m)**2}END{print sqrt(s/(NR*(NR-1)))}' ./MC-runall.txt )
      echo $iter $avg $std
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
