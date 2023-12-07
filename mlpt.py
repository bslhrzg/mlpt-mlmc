#!/usr/bin/env python3


from ase.io import read
import numpy as np 
from sys import argv
from dscribe.descriptors import SOAP
from sklearn.preprocessing import normalize
from dscribe.kernels import REMatchKernel
from sklearn.kernel_ridge import KernelRidge
from joblib import Parallel, delayed
from os import system

re = REMatchKernel( metric="linear", gamma = 0.01, alpha=1, threshold=1e-6 )
KR = KernelRidge( alpha=1e-6, kernel="precomputed" )

conv = 23.061

trj   = []
e_low = []
x_train = []
y_train = []
y_shift = 0.0e0
e_rew = 0.0e0 
e_bia = 0.0e0 
Iw_index = 0.0e0

fsample = 'NULL'
ftrain = 'NULL'
ftest = 'NULL'
kT = 0.02585e0
ncpu = 1 

fmod = []

###############################################################################
def main() :

  # collect input variables from input line
  get_input_line()
  # getting the trajectories and source energies
  get_main_data()
  # getting the descriptor and the energy train 
  get_train_data()
  # training the model
  training()
  # predictions in parallel
  prediction()
  # printing the output
  report()


###############################################################################
def report():

  fo = open('Output.txt','w')
  fo.write("\t   Biased: %12.6f eV\n"%(e_bia) )
  fo.write("\t Reweight: %12.6f eV\n\n"%(e_rew) )
  fo.write("\t Iw Index:  %8.6f\n"%(Iw_index) )

  if ftest != 'NULL' :
    idx   = np.loadtxt(ftest+'.txt',usecols=0,dtype=int)
    e_hig = np.loadtxt(ftest+'.txt',usecols=1,dtype=np.float128)
    err = np.zeros(len(idx),dtype=np.float128)
    j = 0 
    for i in idx :
      err[j] = e_pred[i] - e_hig[j]
      j += 1
    rmse = np.sqrt( np.mean( np.square( err ) ) ) 
    merr = np.mean( err ) 
    fo.write( "\n<<<< TEST ERROR >>>>\n\n")
    fo.write( "\t RMS: %12.6f eV\t%12.6f kcal/mol\n"%( rmse, rmse*conv ) )
    fo.write( "\tMean: %12.6f eV\t%12.6f kcal/mol"%( merr, merr*conv ) )
  fo.flush()
  fo.close()

  fo = open('Energies.txt','w')
  fo.write('#E0  #E_pred\n')
  for i in range(len(trj)):
    fo.write(' %20.10f %20.10f\n'%(e_low[i],e_pred[i]))
  fo.flush()
  fo.close()


###############################################################################
def prediction():

  global e_pred, delta_E
  global e_rew, e_bia, Iw_index

  ntot = len( trj )

  dE = np.memmap( 'tmp-ediff', dtype=np.float, shape=ntot, mode='w+' )

  ndim = len( y_train )
  nblk = int( ntot / ndim )
  lint = nblk*ndim

  Parallel( n_jobs=ncpu, backend="multiprocessing" ) ( delayed(fpred) (i,nblk,dE) for i in range(0,lint,nblk) )

  if ( len( trj[lint:] ) != 0 ):
      xp = get_desc( trj[lint:] )
      KP = re.create( x_train, xp )
      dE[lint:] = fmod.predict(KP.T)

  delta_E = dE + y_shift
  e_pred  = e_low + delta_E 

  e_bia = np.mean( e_pred )
  e_rew, Iw_index = reweight( e_pred, delta_E )

  system('rm tmp-ediff')


###############################################################################
def fpred( ii, nsize, res ):
  xp = get_desc( trj[ii:ii+nsize] )
  KP = re.create( x_train, xp )
  res[ii:ii+nsize] = fmod.predict(KP.T)


################################################################################
def reweight( v, dv ):
  beta = 1.0e0 / kT
  y   = dv - min(dv)
  fy  = np.exp( -beta*y )
  sy  = 1.0e0/np.sum(fy)
  wy  = fy*sy
  erw = np.sum( v*wy )
  Wo = np.sort( fy )
  s1 = np.sum( Wo )
  s2 = 0.0e0
  ndim = len( Wo )
  for i in range( ndim ):
    s2  = s2 + Wo[i]
    if ( s2 / s1 ) >= 0.50e0:
      break
  mdim = i - 1
  Iw_idx = np.float( ( ndim - mdim ) / ndim )
  return erw,Iw_idx


################################################################################
def training():
  global fmod
  KT = re.create( x_train )
  fmod = KR.fit( KT, y_train )


###############################################################################
def get_train_data():

  global x_train, y_train, y_shift

  idx = np.loadtxt(ftrain+'.txt',usecols=0,dtype=int)
  e1 = np.loadtxt(ftrain+'.txt',usecols=1,dtype=np.float128)
  y_train = np.zeros(len(e1),dtype=np.float128)

  at = []
  j = 0 
  for i in idx:
    at.append(trj[i])
    y_train[j] = e1[j] - e_low[i]
    j += 1 

  y_shift  = np.mean( y_train )
  y_train -= y_shift

  x_train = get_desc ( at )


###############################################################################
def get_desc( at ) :
  anum = list( at[0].get_atomic_numbers() )
  dc = SOAP( species=anum, periodic=True,
             rcut=5.0, nmax=8, lmax=6, sigma=1,
             crossover=True, sparse=False )
  x = []
  for i in range(len(at)):
    xi = dc.create( at[i] )
    xi = normalize( xi )
    x.append(xi)
  return x


###############################################################################
def get_main_data() :
  global trj, e_low
  trj = read( fsample+'.xyz', index=':', format='extxyz' )
  e_low = np.loadtxt( fsample+'.txt', usecols=0, dtype=np.float128 )


###############################################################################
def get_input_line() :
  global fsample, ftrain, ftest
  global kT, ncpu
  fsample = argv[1]
  ftrain = argv[2]
  kT = float( argv[3] )
  ftest = argv[4]
  ncpu = int( argv[5] )


###############################################################################
if __name__ == '__main__' :
  main()
