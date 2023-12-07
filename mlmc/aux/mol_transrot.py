#!/usr/bin/env python3

from ase.io import read, write
from ase.io.vasp import write_vasp 
import numpy as np 
import argparse

one = 1.0e0
zero = 0.0e0

center = -1
imol = [] 
max_ds = 0.1e0
fposcar = 'POSCAR' 
rot_lim = 15.0 
trans_lim = 0.5e0 
fout='out.POSCAR.vasp'

###############################################################################
def main() :

  parser()
  rnd = get_random()
  at = read(fposcar,index=':',format='vasp')
  dir_lat = at[0].get_cell()
  cry = at[0].get_scaled_positions(wrap=True)
  cry_uwr = unwrap_mol(cry)
  xyz = rec2car( cry_uwr, dir_lat )
  at_sym = at[0].get_chemical_symbols()
  xyz = translation( xyz, rnd[0], rnd[1], rnd[2] )
  xyz = rotations( xyz, rnd[3], rnd[4], rnd[5] )
  at[0].set_positions(xyz)
  write('MOL.xyz',at,format='extxyz',append=True)
  write_vasp(fout,at,label='System',vasp5=True)


###############################################################################
def myprt(a,x):
  nd = x.shape[0]
  print("{}".format(nd))
  print("")
  for i in range(nd):
   si = a[i]
   xi = x[i,0]
   yi = x[i,1]
   zi = x[i,2]
   print("%2s%12.6f%12.6f%12.6f"%(si,xi,yi,zi))


###############################################################################
def car2rec( v, u ) :
  nd = v.shape[0]
  u0 = u[0,:]
  u1 = u[1,:]
  u2 = u[2,:]
  w = np.zeros((nd,3),dtype=np.float)
  for i in range(nd) :
    vi = v[i,:]
    w[i,0] = np.dot( vi, u0 )
    w[i,1] = np.dot( vi, u1 )
    w[i,2] = np.dot( vi, u2 )
  return w


###############################################################################
def rec2car( v, u ) :
  nd = v.shape[0]
  u0 = u[:,0]
  u1 = u[:,1]
  u2 = u[:,2]
  w = np.zeros((nd,3),dtype=np.float)
  for i in range(nd) :
    vi = v[i,:]
    w[i,0] = np.dot( vi, u0 )
    w[i,1] = np.dot( vi, u1 )
    w[i,2] = np.dot( vi, u2 )
  return w


###############################################################################
def reciprocal(u) :
  v = np.zeros((3,3),dtype=np.float)
  v[0,:] = np.cross(u[1,:],u[2,:])
  v[1,:] = np.cross(u[2,:],u[0,:])
  v[2,:] = np.cross(u[0,:],u[1,:])
  omega = v[0,0]*u[0,0]+v[0,1]*u[0,1]+v[0,2]*u[0,2]
  v = v / omega
  return v, omega


################################################################################
def unwrap_mol( x_in ) :
  x_out = x_in
  nato = len(imol)
  d = np.zeros((nato,3),dtype=np.float)
  c = np.zeros((nato,3),dtype=np.float)
  cm = np.zeros(3,dtype=np.float)
  cm = x_out[center,:]
  for i in range(nato):
    for j in range(3) :
      c[i,j] = x_out[imol[i],j]
      d[i,j] = c[i,j] - cm[j]
      if np.fabs(d[i,j]) > 0.5:
        if d[i,j] > 0.0 : 
          d[i,j] = -1.0
        else :
          d[i,j] = 1.0
        x_out[imol[i],j] = c[i,j] + d[i,j]
  return x_out


###############################################################################
def parser():

  pr = argparse.ArgumentParser(description='Adding Rot-Translational movement')

  pr.add_argument('-p','--poscar', action='store', type=str, required=True,
                  help='POSCAR filename' )
  pr.add_argument('-r','--rot_lim', action='store', type=float, required=True, 
                  help='random limits for the rotation [-rot_lim,rot_lim]')
  pr.add_argument('-t','--trans_lim', action='store', type=float, required=True, 
                  help='random limits for the translation [-trans_lim,trans_lim]')
  pr.add_argument('-d','--max_ds', action='store', type=float, required=True,
                 help='max. translation accepted step dS=sqrt(x**2+y**2+z**2)' )
  pr.add_argument('-c','--center', action='store', type=int, required=True, 
                 help='central atom index for the rotations' )
  pr.add_argument('-i','--imol', action='store', type=int, nargs='+', required=True,
                 help='[at1,at2,at3,...] index vector of the molecular fragment' )
  pr.add_argument('-o','--out', action='store', type=str, required=True,
                 help='POSCAR output filename' )

  args = pr.parse_args()

  global fout, center, imol, max_ds, fposcar, rot_lim, trans_lim

  center=args.center
  imol=args.imol
  max_ds=args.max_ds
  fposcar=args.poscar
  rot_lim=args.rot_lim 
  trans_lim=args.trans_lim
  fout=args.out


###############################################################################
def translation( q, Tx, Ty, Tz ):
  for i in range(len(imol)):
    j = imol[i]
    q[j,0] += Tx
    q[j,1] += Ty
    q[j,2] += Tz
  return q 


###############################################################################
def rotations( cart, alpha, beta, gamma ):

  # Eulerian rotational matrix
  ca = np.cos(alpha)
  sa = np.sin(alpha)
  cb = np.cos(beta)
  sb = np.sin(beta)
  cg = np.cos(gamma)
  sg = np.sin(gamma)
  R = np.zeros((3,3),dtype=float)
  R[0,0]=ca*cb
  R[0,1]=ca*sb*sg-sa*cg
  R[0,2]=ca*sb*cg+sa*sg
  R[1,0]=sa*cb
  R[1,1]=sa*sb*sg+ca*cg
  R[1,2]=sa*sb*cg-ca*sg
  R[2,0]=-sb
  R[2,1]=cb*sg
  R[2,2]=cb*cg

  nd = len(imol)
  q = np.zeros((nd,3),dtype=np.float)
  for i in range(nd) :
    j = imol[i]
    q[i,0] = cart[j,0]
    q[i,1] = cart[j,1]
    q[i,2] = cart[j,2]

  xcm = cart[center,0]
  ycm = cart[center,1]
  zcm = cart[center,2]

  q[:,0] -= xcm
  q[:,1] -= ycm
  q[:,2] -= zcm

  for i in range(nd):
    x = q[i,0]*R[0,0] + q[i,1]*R[1,0] + q[i,2]*R[2,0]
    y = q[i,0]*R[0,1] + q[i,1]*R[1,1] + q[i,2]*R[2,1]
    z = q[i,0]*R[0,2] + q[i,1]*R[1,2] + q[i,2]*R[2,2]
    q[i,0] = x
    q[i,1] = y
    q[i,2] = z

  q[:,0] += xcm
  q[:,1] += ycm
  q[:,2] += zcm

  for i in range(nd) :
    j = imol[i]
    cart[j,0] = q[i,0]
    cart[j,1] = q[i,1]
    cart[j,2] = q[i,2]

  return cart 


###############################################################################
def get_random() :
  v=np.zeros(6,dtype=float)
  i = 0 
  ds = max_ds
  while ds >= max_ds :
    i+=1
    t1 = np.random.uniform(-trans_lim,trans_lim)
    t2 = np.random.uniform(-trans_lim,trans_lim)
    t3 = np.random.uniform(-trans_lim,trans_lim)
    ds = np.sqrt( (t1*t1) + (t2*t2) + (t3*t3) )
  v[0] = t1
  v[1] = t2
  v[2] = t3
  for i in range(3,6):
    vi = np.random.uniform(-rot_lim,rot_lim)
    vi = np.pi*vi/ 180.0
    v[i] = vi
  return v


###############################################################################
if __name__ == "__main__" :
  main()
