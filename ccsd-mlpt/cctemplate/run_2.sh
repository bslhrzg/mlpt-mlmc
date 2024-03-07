#!/bin/bash

#SBATCH -p tarasall
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --job-name=co2
#SBATCH --output=output-%x-%N.%j
#SBATCH --error=error-%x-%N.%j
#SBATCH --mem=60GB

ulimit -s unlimited

module load compiler/intel19u5 openmpi/4.1.1-intel19u5

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1


VASPBIN="mpirun -np $SLURM_NTASKS /auto/tms7/herzog1/build/ccsd_ml/tests/ccsdt-cc4s-fast/bin/vasp_gam"

CC4S="mpirun -np $SLURM_NTASKS /auto/tms7/herzog1/build/ccsd_ml/tests/ccsdt-cc4s-fast/bin/Cc4s"


enc=400
egwno=150
egw=150



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

cp WAVECAR.FNO WAVECAR

cat >INCAR <<!
LRHFATM=.TRUE.
NBANDS = $nbno
NBANDSHIGH = $nbno
ALGO = sub ; NELM = 1
ISMEAR = -1 ; SIGMA = 0.00001
EDIFF = 1E-4
LHFCALC=.TRUE.
AEXX=1.0
ENCUT = $enc
!
rm WAVEDER
cat INCAR
$VASPBIN
cp OUTCAR OUTCAR.HFTdiag.no.$nbno.$pos

for egw in 300
do

cat >INCAR <<!
LRHFATM=.TRUE.
NBANDS = $nbno
NBANDSHIGH = $nbno
ALGO = MP2 ; NELM = 1
ISMEAR = -1 ; SIGMA = 0.00001
EDIFF = 1E-5
LHFCALC=.TRUE.
AEXX=1.0
ENCUT = $enc
ENCUTGW = $egw
LMAXFOCKAE=4
!
rm WAVEDER
cat INCAR
$VASPBIN
cp OUTCAR OUTCAR.MP2.latt.$a.egw.$egw.nbno.$nbno.$pos

cat >INCAR <<!
LRHFATM=.TRUE.
NBANDS=$nbno
NBANDSHIGH = $nbno
ENCUT = $enc
ENCUTGW = $egw
SIGMA = 0.0001
EDIFF = 1E-3
PREC=A
PRECFOCK=A
NELM=100000
LHFCALC=.TRUE.
AEXX=1.0
ALGO=CC4S
ISYM=-1
#LFTODDUMP=.TRUE.
!
cat INCAR
$VASPBIN
cp OUTCAR OUTCAR.FTODDUMP.enc.$enc.nbno.$nbno.$pos

#cp cc4s_not.in cc4s.in


#$CC4S > cc4s.out.$a.enc.$enc.egw.$egw.nbno.$nbno.$pos.precfocka.1em4

done

done

done
