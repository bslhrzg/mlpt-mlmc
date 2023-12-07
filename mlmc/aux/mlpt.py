import numpy as np
from joblib import dump, load
import ase
from ase.io import read, write
from dscribe.kernels import REMatchKernel
from sklearn.kernel_ridge import KernelRidge
from dscribe.descriptors import SOAP
from sklearn.preprocessing import normalize

re = REMatchKernel( metric="linear", gamma = 0.01, alpha=0.1, threshold=1e-6 )
KR = KernelRidge( alpha=1e-8, kernel="precomputed" )

class MLPT:

  def __init__( self ):
    pass 
#
# training routine
#
  def training( self ):
    idx = np.loadtxt( 'TRAIN.dat', usecols=[0], dtype=np.int )
    e1  = np.loadtxt( 'TRAIN.dat', usecols=[1], dtype=np.float )
    at = read('TRAIN.xyz',index=':',format='extxyz')
    ndim = len(at)
    e0 = np.zeros( ndim, dtype=np.float )
    for i in range( ndim ):
      e0[i] = at[i].get_potential_energy()
    Y_trn = e1 - e0
    mY_trn = np.mean( Y_trn )
    dump(mY_trn,'Ytrn_Mean.bin')
    Y_trn = Y_trn - mY_trn
    X_trn = self.build_soapN( at ) 
    dump(X_trn,'SOAP_desc.bin')
    KTrn  = re.create( X_trn )
    model = KR
    fmod  = model.fit( KTrn, Y_trn )
    dump(model,'Fitted_Model.bin')
    return 0
#
# prediction routine
#
  def prediction( self,fio,e0 ) :
    X_trn  = load('SOAP_desc.bin') 
    model = load('Fitted_Model.bin')
    mY_trn = load('Ytrn_Mean.bin')
    at = read(fio,index=':',format='vasp')
    X_prd = self.build_soapN( at )
    KPrd  = re.create( X_trn , X_prd )
    Y_prd = model.predict( KPrd.T )
    e1 = e0 + Y_prd + mY_trn
    return e1[0]
#
# SOAP N descriptor
#
  def build_soapN( self, atoms ) :
    at_num = list( atoms[0].get_atomic_numbers() )
    desc = SOAP( species=at_num, periodic=True, 
               rcut=5.227, nmax=8, lmax=6, sigma=0.0046, 
               crossover=True, sparse=False )
    soap_list = []
    for j in range(0,len(atoms)):
      soap = desc.create(atoms[j])
      soap = normalize(soap)
      soap_list.append(soap)
    return soap_list
