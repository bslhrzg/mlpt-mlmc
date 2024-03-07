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

cp $pos POSCAR
rm WAVECAR

cat >INCAR <<!
ENCUT=$enc
EDIFF  = 1.0E-05
NBANDS = 768 
NELMDL= -5
NPAR=4 
!
cat INCAR
$VASPBIN
cp OUTCAR OUTCAR.DFT.$pos


cat >INCAR <<!
LRHFATM=.TRUE.
ENCUT=$enc
EDIFF  = 1.0E-05
LHFCALC=.TRUE.
AEXX=1.0
ALGO=C
SIGMA=0.001
NPAR=4
!
cat INCAR
$VASPBIN
cp OUTCAR OUTCAR.$k.HFT.$pos


nb=`awk <OUTCAR "/maximum number of plane-waves:/ { print \\$5*2-1 }"`
nocc=`awk <OUTCAR "/NELEC/ { print \\$3/2 }"`
nfpa=`awk <OUTCAR "/NELEC/ { print \\$3/2*6 }"`
nfpb=`awk <OUTCAR "/NELEC/ { print \\$3/2*11 }"`
nfpc=`awk <OUTCAR "/NELEC/ { print \\$3/2*21 }"`
nfpd=`awk <OUTCAR "/NELEC/ { print \\$3/2*31 }"`

echo "going to use $nb bands"
echo "there are $nocc occupied bands"
echo "for fpa we need $nfpa  bands"
echo "for fpb we need $nfpb  bands"


cat >INCAR <<!
LRHFATM=.TRUE.
NBANDS = $nb
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
cp OUTCAR OUTCAR.HFTdiag.$pos

cat >INCAR <<!
LRHFATM=.TRUE.
NBANDS = $nb
ALGO = MP2 ; NELM = 1
ISMEAR = -1 ; SIGMA = 0.00001
EDIFF = 1E-4
LHFCALC=.TRUE.
AEXX=1.0
ENCUT = $enc
LSFACTOR=.TRUE.
LMAXFOCKAE=4
!
rm WAVEDER
cat INCAR
$VASPBIN
cp OUTCAR OUTCAR.MP2full.$pos

cp Mp2PairEnergies.dat Mp2PairEnergies.dat.$pos 
cp Mp2PairEnergies.yaml Mp2PairEnergies.yaml.$pos 

cat >INCAR <<!
LRHFATM=.TRUE.
NBANDS = $nb
ALGO = MP2NO ;
LAPPROX=.TRUE.
ISMEAR = -1 ; SIGMA = 0.00001
EDIFF = 1E-5
LHFCALC=.TRUE.
AEXX=1.0
ENCUT = $enc
ENCUTGW = $egwno
LMAXFOCKAE=4
!
rm WAVEDER
cat INCAR
$VASPBIN
cp OUTCAR OUTCAR.MP2NO.$pos

#for nbno in 6000 6800 7600
for nbno in $nfpa 
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

cp cc4s_t.in cc4s.in 

module purge
module use -a /auto/tms7/chagasda1/modulefiles
module load compiler/gnu-8.4.0 openmpi/gnu-8.4.0/4.1.1
module load OpenBLAS/gnu-8.4.0/0.3.13

$CC4S > cc4s.out.$a.enc.$enc.egw.$egw.nbno.$nbno.$pos.precfocka.1em4

module purge
module load compiler/intel19u5 openmpi/4.1.1-intel19u5

done

done

done
