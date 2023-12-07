#!/usr/bin/env python2.7

from numpy import *
import sys,time

import string,random,math
import readandwrite,intcoord,takeinp,bmatrix,mymath
import datastruct,various,physconstants
import takeinpPOSCAR
import inputer
import int_values


def metropolis_criter(e1,e0,beta):
  t=1
  if e1>e0:
    b=(e0-e1)*beta
    b=e**b
    r=random.random()
    if r>b:
      t=0
  return t

def back_transf(toler,maxiter,xi,qi,qj,lattmat,primcoords,numofatoms,transmat_i,fname):
  # qj=q0=q(x0)=q(xj)

  maxstep=0.3
  rrr="   >Back-transformation"
  fname.write(rrr+"\n")
  rrr="            step       err   coord. num."
  fname.write(rrr+"\n")

  err=1000.
  qk=1*qi
  qk=sanity_operation_simple(qj,qk,primcoords)
  xk=1*xi
  count=0
  dq=qk-qj
  err=max(abs(dq))
  btscale=1.
  if err>0.: btscale=maxstep/err
  btscale=min(1.,btscale)
  #print btscale
  errarg=argmax(abs(dq))
  rrr="    b_tr>"+"   "+'%4d'%(count+1) + "    "+'%.3e'%(err) +'%4d'%(errarg+1)
  fname.write(rrr+"\n")
  while err> toler:
    count+=1
    if count>maxiter:
      sys.exit("ERROR: Back-transformation failed")
    b=bmatrix.bmatrix(xk,primcoords,numofatoms,lattmat,0)
    Bmat=b.Bmatrix
    Bmat=dot(Bmat,transmat_i)
    Gmatrix=dot(transpose(Bmat),Bmat)
    GmatrixInv=mymath.mygeneralized_inverse(Gmatrix)
    Amat=dot(GmatrixInv,transpose(Bmat))
    dx=btscale*dot(Amat,dq)
    xk=xk-dx
    x3=various.format_change(xk)
    qk=give_internal_coords(x3,primcoords,lattmat)
    qk=sanity_operation_simple(qj,qk,primcoords)
    dq=qk-qj
    err=max(abs(dq))
    errarg=argmax(abs(dq))
    rrr="    b_tr>"+"   "+'%4d'%(count+1) + "    "+'%.3e'%(err) +'%4d'%(errarg+1)
    fname.write(rrr+"\n")  
  fname.write("\n")
    #print "bt:",err
    
  return xk

def shake_algo(toler,maxiter,x,xNext,vNext,lattmat,primcoords,q0,masses,numofatoms,transmat_i,dt,fname):
  rrr="   >Shake-algorithm"
  fname.write(rrr+"\n")
  rrr="            step       err"
  fname.write(rrr+"\n")

  b=bmatrix.bmatrix(x,primcoords,numofatoms,lattmat,0)
  Bmat=b.Bmatrix
  #Bmat=setup_iconst(primcoords,Bmat)
  
  Bmat0m=dot(Bmat,transmat_i)
  for i in range(len(Bmat0m)):
    Bmat0m[i]/=masses

  gamma0=zeros(len(primcoords))
  cforce=zeros(3*numofatoms)
  dq=zeros(len(primcoords))
  err=10000.
  count=0

  while err>toler:
    count+=1
    if count>maxiter: 
      sys.exit("ERROR: shake algorithm failed")
    b=bmatrix.bmatrix(xNext,primcoords,numofatoms,lattmat,0)
    Bmat=b.Bmatrix
    Bmat=dot(Bmat,transmat_i)
    x3=various.format_change(xNext)
    qNext=give_internal_coords(x3,primcoords,lattmat)
    dq=qNext-q0
    for i in range(len(primcoords)):
      if primcoords[i].status!="Keep":
        dq[i]=0.
    err=max(abs(dq))
    rrr="    s_al>"+"   "+'%4d'%(count) + "    "+'%.3e'%(err)
    fname.write(rrr+"\n")
    #print err
    gamma=zeros(len(primcoords))
    for i in range(len(primcoords)):
      if primcoords[i].status=="Keep":
        gamma[i]=dq[i]/dot(Bmat0m[i],Bmat[i])/dt**2
        #gamma0[i]=gamma0[i]+dq[i]/dot(Bmat0m[i],Bmat)/dt**2

        cforce_=-gamma[i]*Bmat0m[i]
        cforce+=cforce_
        xNext+=cforce_*dt**2

  fname.write("\n")
  vNext+=cforce*dt

  return xNext,vNext


def sanity_operation(prims1,prims2,coords,Bmat,step,freport):
  maxval=3.10
  for m in range(len(coords)):
    if coords[m].dtyp=='simple':
      if coords[m].tag=='A' or coords[m].tag=='T':
        while (prims2[m]-prims1[m])>pi:
          prims2[m]=prims2[m]-2*pi
        while (prims2[m]-prims1[m])<-pi:
          prims2[m]=prims2[m]+2*pi
      if coords[m].tag=='A':
        if prims2[m]>maxval:
          freport.write("\n")
          rrr='   Warning(MD step'+str(step)+'): invalid value of coord. '+str(m)+', '+str(prims2[m])+', temporarily switched off'
          freport.write(rrr+"\n")
          freport.write("\n")
          print rrr
          prims2[m]=1*prims1[m]
          Bmat[m]=0.
        prims2[m]=abs(prims2[m])
    if coords[m].dtyp=='tsum':
      dmin=prims2[m]-prims1[m]
      for i in range(len(coords[m].tag)):
        if coords[m].tag[i]=='T':
          dminP=dmin+2*pi*coords[m].coefs[i]
          if abs(dminP)<abs(dmin): dmin=dminP
          dminP=min(dmin,dmin-2*pi*coords[m].coefs[i])
          if abs(dminP)<abs(dmin): dmin=dminP
      prims2[m]=prims1[m]+dmin

  return prims2

def sanity_operation_simple(prims1,prims2,coords):
  maxval=3.10
  for m in range(len(coords)):
    if coords[m].dtyp=='simple':
      if coords[m].tag=='A' or coords[m].tag=='T':
        while (prims2[m]-prims1[m])>pi:
          prims2[m]=prims2[m]-2*pi
        while (prims2[m]-prims1[m])<-pi:
          prims2[m]=prims2[m]+2*pi
    if coords[m].dtyp=='tsum':
      dmin=prims2[m]-prims1[m]
      #print "dmin1",dmin
      for i in range(len(coords[m].tag)):
        if coords[m].tag[i]=='T':
          dminP=dmin+2*pi*coords[m].coefs[i]
          if abs(dminP)<abs(dmin): dmin=dminP
          dminP=min(dmin,dmin-2*pi*coords[m].coefs[i])
          if abs(dminP)<abs(dmin): dmin=dminP
      #print "dmin2",dmin
      prims2[m]=prims1[m]+dmin
  return prims2


def give_internal_coords(coords_c,primcoords,lattmat):
  q=zeros(len(primcoords))
  deal=int_values.Int_values(coords_c,primcoords,lattmat,[])
  for i in range(len(primcoords)):
    q[i]=primcoords[i].value
  return q

def give_whattags(what,tags,atomictags):
  whattags=[None]
  if tags=='R' or tags=='A' or tags=='T' or tags=='tV' or tags=='IR1'\
  or tags=='IR6' or tags=='RatioR' or tags=='X' or \
  tags=='Y' or tags=='Z' or tags=='fX' or tags=='fY' or tags=='fZ':
    whattags=[]
    for i in range(len(what)):
      whattags.append(atomictags[what[i]])
  return whattags

def get_iconst(crt,lattmat,atomictags,fname):
  coords=[]
  try:
    consttags,constwhat,constwhere,constcoefs,conststat,complextype=\
    readandwrite.read_cwhatwhere(fname,crt,lattmat)

  except IOError:
    numconst=0
    constcoefs=[]
    conststat=[]
  else:
    numconst=len(constwhat)
  if numconst>0:
    if constcoefs==[]:
      for i in range(numconst):
        constwhattags=give_whattags(constwhat[i],consttags[i],atomictags)
        coords.append(datastruct.Complextype('simple',[1],consttags[i],constwhat[i],\
        constwhattags,constwhere[i],0.0,'Ti'))
    else:
      for i in range(len(constcoefs)):
        if len(constcoefs[i])!=len(consttags):
          print 'incorect definition of complex constrained coordinate!!!'
        else:
          if len(constcoefs[i])==1:
            constwhattags=give_whattags(constwhat[i],consttags[i],atomictags)
            coords.append(datastruct.Complextype('simple',constcoefs[i],consttags[i],constwhat[i],\
            constwhattags,constwhere[i],0.0,conststat[i]))
          else:
            coords.append(datastruct.Complextype(complextype[i],constcoefs[i],consttags,constwhat,\
            [None],constwhere,0.0,conststat[i]))
  return coords,conststat

def setup_iconst(primcoords,Bmatp):
  Cmatp=zeros((len(Bmatp),len(Bmatp[0])),float)
  constraints=[]
  for i in range(len(primcoords)):
    if primcoords[i].status=="Keep":
      constraints.append(i)
  for i in constraints:
    Cmatp[i]=Bmatp[i]
  Cmatp=mymath.orthonormalize_mat(Cmatp)
  oBmatp=dot(Bmatp,transpose(Cmatp))
  for i in range(len(oBmatp)):
    oBmatp[i][i]=0
  oBmatp=dot(oBmatp,Cmatp)
  Bmatp=Bmatp-oBmatp
  return Bmatp

def setup_iconst2(primcoords,Bmatp):
  for i in range(len(primcoords)):
    if primcoords[i].status!="Ti":
       Bmatp[i]=0.
  return Bmatp

def read_hessian(numofatoms,file):
  f=open(file,'r')
  energy0=0.
  xd0=zeros((numofatoms,3))
  fconst=zeros(3*numofatoms)
  modes=zeros((3*numofatoms,3*numofatoms))
  hessian=zeros((3*numofatoms,3*numofatoms))

  indx2=-1
  count=0
  for line in f.readlines():
    prefield=line.split()
    count+=1
    if count==1:
      energy0=float(prefield[0])
    elif count<=(numofatoms+1):
      xd0[count-2]=[float(prefield[0]),float(prefield[1]),float(prefield[2])]       
    else:
      indx1=count%(numofatoms+1)-1
      #indx2=count/(numofatoms+1)-1
      if indx1==-1: indx1=numofatoms #;indx2-=1
      if indx1==0:
        indx2+=1
        fconst[indx2]=float(prefield[0])
      else:
        #print indx1,indx2,prefield
        modes[indx2][3*(indx1-1)]=float(prefield[0])
        modes[indx2][3*(indx1-1)+1]=float(prefield[1])
        modes[indx2][3*(indx1-1)+2]=float(prefield[2])
  #print fconst
  f.close()
 
  hessian=zeros((3*numofatoms,3*numofatoms))
  for i in range(3*numofatoms):
    hessian[i][i]=fconst[i]

  hessian=dot(hessian,modes)
  hessian=dot(transpose(modes),hessian)
  
  #eigval,eigvect=mymath.eigh_tb(hessian)
  #print eigval

  return energy0,xd0,hessian

def boltzmannRandom(width):
  #random.seed()
  x=random.random()
  y=random.random()
  z=width*cos(2*pi*x)*(-2*log(y))**0.5
  return z

def calcMomentum(v,m):
  p=zeros(len(v),float)
  for i in range(len(v)):
    p[i]=v[i]*m[i/3]
  return p

def cmVelocity(vel,masses):
   """calculates velocity of the
   ce nter of masses
   """
   cmvel=zeros(3,float)
   totmass=zeros(3,float)
   for i in range(len(masses)):
     for j in range(3):
       cmvel[j]+=vel[3*i+j]*masses[i]
       totmass[j]+=masses[i]
   cmvel=cmvel/totmass
   return cmvel

def brandNewVelocity(meanT,allmasses,pC):
  kb=pC.bolk
  velo_=zeros(3*len(allmasses),float)
  for i in range(len(velo_)):
    velo_[i]=boltzmannRandom((meanT*kb/allmasses[i/3])**0.5)
  # m/s->A/fs:
  velo_*=1e-5
  return array(velo_)

def mdStep(x,v,a,dt,m,temperature,colprob,pC):
  #kb=pC.bolk
  #totalMass=sum(m)
  #veloCM=zeros(3,float)
  #veloCM[0]=boltzmannRandom((temperature*kb/totalMass)**0.5)
  #veloCM[1]=boltzmannRandom((temperature*kb/totalMass)**0.5)
  #veloCM[2]=boltzmannRandom((temperature*kb/totalMass)**0.5)
  # m/s->A/fs:
  #veloCM*=1e-5

  #for i in range(len(v)/3):
  #  v[3*i+0]+=veloCM[0]
  #  v[3*i+1]+=veloCM[1]
  #  v[3*i+2]+=veloCM[2]

  #v=v+a*dt
  #v_=brandNewVelocity(temperature,m,pC)
  #for i in range(len(v_)):
  #  random.seed()
  #  rn=random.random()
  #  if rn<colprob:
  #    v[i]=v_[i]
  cmvel=cmVelocity(v,m)
  v_=1*v
  for i in range(len(v)/3):
    v_[3*i+0]=v[3*i+0]-cmvel[0]
    v_[3*i+1]=v[3*i+1]-cmvel[1]
    v_[3*i+2]=v[3*i+2]-cmvel[2]
  x=x+v_*dt
  return x,v


def computeForce(x,x0,h):
  dx=x-x0
  f=-dot(h,dx)
  energy=-0.5*sum(dx*f)
  return f,energy

def computeAccel(f,m):
  a=zeros(len(f),float)
  for i in range(len(m)):
    a[3*i:3*i+3]=f[3*i:3*i+3]/m[i]
  return a

def writeStr(x,flags):
  print len(flags)
  print ''
  for i in range(len(flags)):
    print flags[i],x[3*i],x[3*i+1],x[3*i+2]

def writeHeaderPoscar(f,inpt):
  f.write(inpt.comment+'\n')
  f.write(str(1.000)+'\n')
  for i in range(len(inpt.lattmat)):
    row='   '+'%.12f'%(round(inpt.lattmat[i][0],12))+'   '+'%.12f'%(round(inpt.lattmat[i][1],12))+'   '+'%.12f'%(round(inpt.lattmat[i][2],12))+'\n'
    f.write(row)
  
  row="" 
  for i in range(inpt.ntypes):
    row+=inpt.atomicFlags[i]+" "
  row+='\n'
  f.write(row)
 
  row=""
  for i in range(inpt.ntypes):
    row+=str(inpt.types[i])+" "
  row+='\n'
  f.write(row)


def writeXDATCAR(step,x,latinv,f):
  rrr="Direct configuration= "+str(step)+'\n'
  f.write(rrr)
  rrr="  "
  for i in range(len(x)/3):
    dx=array([x[3*i],x[3*i+1],x[3*i+2]])
    dx=dot(dx,latinv)
    #print dx[0],dx[1],dx[2]
    rrr="   "+'%.12f'%(round(dx[0],12))+" "+'%.12f'%(round(dx[1],12))+" "+'%.12f'%(round(dx[2],12))+'\n'
    f.write(rrr)

def writeOSZICAR(step,temperature,energyT,energyP,energyK,f):
  rrr=str(step)+"\t"+"T=  "+'%.3f'%(round(temperature,3))+"  " + "E=  "+'%.6e'%(round(energyT,6))
  rrr+="  " + "F=  "+'%.6e'%(round(energyT,6))
  rrr+="  " + "E0=  "+'%.6e'%(round(energyP,6))+"  "+"EK= "+'%.6e'%(round(energyK,6))
  rrr+='\n'
  f.write(rrr)


def writeREPORT_step(step,fname):
  rrr="========================================="
  fname.write(rrr+"\n")
  rrr="         MD step No. "+str(step)
  fname.write(rrr+"\n")
  rrr="========================================"
  fname.write(rrr+"\n")
  fname.write("\n")

def writeREPORT_time(t0,t1,t2,t3,t4,t5,t6,fname):
  rrr="    Timing(s):"
  fname.write(rrr+"\n")
  rrr="      forces1: "+str(t1-t0)
  fname.write(rrr+"\n")
  rrr="      back-transf. : "+str(t2-t1)
  fname.write(rrr+"\n")
  rrr="      forces0: "+str(t3-t2)
  fname.write(rrr+"\n")
  rrr="      md step: "+str(t4-t3)
  fname.write(rrr+"\n")
  rrr="      shake : "+str(t5-t4)
  fname.write(rrr+"\n")
  rrr="      output: "+str(t6-t5)
  fname.write(rrr+"\n")
  rrr="      total: "+str(t6-t0)
  fname.write(rrr+"\n")
  fname.write("\n")

def writeREPORT(step,x,x0,coords,fname):

  rrr="   >Harmonic_coord_0"
  fname.write(rrr+"\n")
  for i in range(len(x)):
    tag="Z"
    if i%3==0: 
      tag="X"
    elif i%3==1:
      tag="Y" 
    indx=i/3+1
    rrr="    hc0> "+str(indx)+" "+tag+"\t"+'%.6f'%(round(x[i],6))+"\t"+'%.6f'%(round(x0[i],6))
    fname.write(rrr+"\n")

  fname.write("\n")



def same_cell(x,y):
  r=0.0
  for i in range(len(x)):
    for j in range(len(x[0])):
      while (x[i][j]-y[i][j])>0.5:
        x[i][j]-=1
      while (x[i][j]-y[i][j])<=-0.5:
        x[i][j]+=1
  return x


#random.seed('ahoj')
#random.seed('tono')

pC=physconstants.physConstants()

#c read basic structural data
#c NOTE that in the present version, we start from 
#c coordinates defined in HESSEMAT, not in POSCAR (the latter are ignored)
tknp=takeinpPOSCAR.TakeInput()
tknp.read("POSCAR")
lattmat=tknp.lattmat
numofatoms=tknp.numofatoms
atomquality=tknp.atomicFlags
katoms=tknp.types
coords_d=tknp.coords_d
coords_c=tknp.coords_c
xd0=coords_d
xc0=dot(xd0,lattmat)
#xd=tknp.coords_d

attags=[]
for i in range(len(katoms)):
  attags+=katoms[i]*[atomquality[i]]

coords_c=dot(coords_d,lattmat)


#c read in the simulation parameters
inps=inputer.inputer()

potim=inps.POTIM
nsw=inps.NSW
equilibration=inps.EQUIL
temperature=inps.TEBEG
colprob=inps.ANDERSEN_PROB
writeIncrem=inps.WINCREM
tiLambda=inps.TILAMBDA
shaketoler=inps.SHAKE_TOLER
bttoler=inps.BT_TOLER
shakemax=inps.SHAKE_MAXITER
btmax=inps.BT_MAXITER



try:
  f=open('MASSES','r')
  mass=[]
  for line in f.readlines():
    line=string.split(line)
    if len(line)>0:mass.append(float(line[0]))
  f.close()

except IOError:
  mass_=inps.POMASS
  mass=[]
  for i in range(len(katoms)):
    mass+=katoms[i]*[mass_[i]]

mass=array(mass)
mass*=pC.amutokg
masses=zeros(3*numofatoms)
for i in range(len(mass)):
  masses[3*i]=mass[i]
  masses[3*i+1]=mass[i]
  masses[3*i+2]=mass[i]


# convert into internals and back to cartesians in order 
# to get rid of translations and rotations
cartesian=various.change_format(coords_c)

primcoords,conststat=get_iconst(cartesian,lattmat,attags,'ICOORD')
q0=give_internal_coords(coords_c,primcoords,lattmat)
#print q0
#print koko

x=1*cartesian
#velocity=brandNewVelocity(temperature,mass,pC)
latinv=linalg.inv(lattmat)

#c output files 
fstruct=open('hCONTCAR','w')
freport=open('hREPORT','w')

rrr=  "           POTIM = "+str(potim)
freport.write(rrr+"\n")
#rrr=  "             NSW = "+str(nsw)
#freport.write(rrr+"\n")
#rrr=  "           EQUIL = "+str(equilibration)
#freport.write(rrr+"\n")
rrr=  "           TEBEG = "+str(temperature)
freport.write(rrr+"\n")
#rrr=  "   ANDERSEN_PROB = "+str(colprob)
#freport.write(rrr+"\n")
#rrr=  "         WINCREM = "+str(writeIncrem)
#freport.write(rrr+"\n")
rrr=  "          POMASS = "
for i in range(len(inps.POMASS)):
  rrr+=str(inps.POMASS[i])+" "
freport.write(rrr+"\n")
#rrr=  "        TILAMBDA = "+str(tiLambda)
#freport.write(rrr+"\n")
rrr=  "        BT_TOLER = "+str(bttoler)
freport.write(rrr+"\n")
rrr=  "      BT_MAXITER = "+str(btmax)
freport.write(rrr+"\n")
rrr=  "     SHAKE_TOLER = "+str(shaketoler)
freport.write(rrr+"\n")
rrr=  "   SHAKE_MAXITER = "+str(shakemax)
freport.write(rrr+"\n")
freport.write("\n")


beta=1./(pC.bolkEV*temperature)

writeHeaderPoscar(fstruct,tknp)

ftest=open('hPOSCAR','w')
writeHeaderPoscar(ftest,tknp)
writeXDATCAR(1,x,latinv,ftest)
ftest.close()
#print x
#print latinv

iconst0=0
iconst3=0
for i in range(len(primcoords)):
  if primcoords[i].status=="Keep": 
    iconst0+=1
  elif primcoords[i].status=="Ti":
    iconst3+=1
  

ndeg=3*numofatoms - iconst0 #-nZeros

rrr=  "   Total number of degrees of freedom: " + str(3*numofatoms)
freport.write(rrr+"\n")
rrr=  "   Number of active egrees of freedom: " + str(ndeg)
freport.write(rrr+"\n")
freport.write("\n")


# no acceleration used in MC
acceleration= zeros(3*numofatoms,float)

x=1*cartesian

t0=0;t1=0;t2=0;t3=0;t4=0;t5=0;t6=0

transmat=mymath.cd_transmatrix(lattmat,len(cartesian))
transmat_i=linalg.inv(transmat)

#main loop
step=1

x_=1*x

writeREPORT_step(step,freport)

# generate new velocities
velocity=brandNewVelocity(temperature,mass,pC)

x,velocity=mdStep(x_,velocity,acceleration,potim,mass,temperature,colprob,pC)

# call the shake algorithm only if there are constraints
if iconst0>0:
  x,velocity=shake_algo(shaketoler,shakemax,x_,x,velocity,lattmat,primcoords,q0,masses,numofatoms,transmat_i,potim,freport)


writeXDATCAR(step,x,latinv,fstruct)
  

fstruct.close()
freport.close()  

