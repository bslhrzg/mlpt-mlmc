#!/bin/bash

#SBATCH -p tarasall
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --job-name=co2
#SBATCH --output=output-%x-%N.%j
#SBATCH --error=error-%x-%N.%j
#SBATCH --mem=60GB

ulimit -s unlimited

module purge
module use -a /auto/tms7/chagasda1/modulefiles
module load compiler/gnu-8.4.0 openmpi/gnu-8.4.0/4.1.1
module load OpenBLAS/gnu-8.4.0/0.3.13

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1

VASPBIN="mpirun -np $SLURM_NTASKS /auto/tms7/herzog1/build/ccsd_ml/tests/ccsdt-cc4s-fast/bin/vasp_gam"

CC4S="mpirun -np $SLURM_NTASKS /auto/tms7/herzog1/build/ccsd_ml/tests/ccsdt-cc4s-fast/bin/Cc4s"


enc=400
egwno=150
egw=150

cp cc4s_not.in cc4s.in

for pos in `ls POSCAR_*`
do


nb=`awk <OUTCAR.$k.HFT.$pos "/maximum number of plane-waves:/ { print \\$5*2-1 }"`
nocc=`awk <OUTCAR.$k.HFT.$pos "/NELEC/ { print \\$3/2 }"`
nfpa=`awk <OUTCAR.$k.HFT.$pos "/NELEC/ { print \\$3/2*6 }"`
nfpb=`awk <OUTCAR.$k.HFT.$pos "/NELEC/ { print \\$3/2*11 }"`
nfpc=`awk <OUTCAR.$k.HFT.$pos "/NELEC/ { print \\$3/2*21 }"`
nfpd=`awk <OUTCAR.$k.HFT.$pos "/NELEC/ { print \\$3/2*31 }"`

echo "going to use $nb bands"
echo "there are $nocc occupied bands"
echo "for fpa we need $nfpa  bands"
echo "for fpb we need $nfpb  bands"


#for nbno in 6000 6800 7600
for nbno in $nfpb
do


for egw in 300
do

$CC4S > cc4s.out.$a.enc.$enc.egw.$egw.nbno.$nbno.$pos.precfocka.1em4

done

done

done
