import os,string
from numpy import *

import sys

class inputer:
 """handles the input information
 """
 EQUIL=0
 TEBEG=300.
 ANDERSEN_PROB=0.00
 WINCREM=100
 TILAMBDA=0.
 SHAKE_TOLER=1e-7
 BT_TOLER=1e-7
 SHAKE_MAXITER=100
 BT_MAXITER=100

 RELAX=0         # 0-realxation of atomic positions only, 1-relaxation of atomic pos.+lat. param.
 CART=0          # optimization in cartesians (1) or in delocalized internals (0)
 TORSIONS=1       # generate torsions (0-no, 1-yes)
 HESSIAN=1       # hessian initialised as a diag. matrix in cartesian (0),
                 # internal coord.(1) space, Lindhs model (2), Fischers model (3)
 HUPDATE=1       # hessian-update formula, 0-no update -GDIIS, 1-BFGS, 2-BFGS-TS, 3-SR1, 4-PSB,
                 # 5-SR1/PSB, 6-SR1/BFGS
 HREFRESH=0      # specifies how offen should be hessian re-initialised (0 - never,1-every step...)
 HUPDIM=2        # dimension of history for BFGS update, can differ from 2 only if HREFRESH=1
 OPTENGINE=0     # "engine" for optimization, 0 -DIIS, 1 -RFO, 3 -pRFO(TS), 4 -qNR
 LINESEARCH=0    # line-search Y/N 1/0
 LINEMAX=5       # max. number of line minimizations
 CONSTCONJ=0     # 0- constraints as defined, 1
 HSING=1e-04     # minimal allowed ratio of eigenvalue of Hessian to its maximal eigenvalue
 POTIM=1         # the value of diagonal component of initial component if self.hessian==0
                 # or if internal and cartesian components are mixed
 KMATAL=0          # calculate K-matrix (0-no, 1-yes)
 NFREE=10        # number of history-steps involved in DIIS, 1(steepest descent)-10
 NSW=1000         # maximal number of relaxation steps
 STEPLIM=0.3     # maximal allowed step in internal coodrs (bohr for bonds, rad. for angles and torsions)
 GCRITER=0.05    # convergence criteria for gradients (eV/A)
 ECRITER=1e-02   # convergence criteria for energy (eV)
 SCRITER=1e-01   # convergence criteria for geometry step (A)
 ASCALE=1.0      # scaling factor for covalent atomic radii. if set to zero, only cartesian coords
                 # will be detected
 BSCALE=2.0      # as ASCALE multiplied bz this if more fragments is present in the cell
 CSCALE=1.0
 ANGLECRIT=6     # criterion for internal-coordinates (bonding angles) detection scheme 
 TORSIONCRIT=4   # criterion for internal-coordinates (torsions) detection scheme
 VDWRAD=0.0      # radius for empirical vdw force field
 FRAGCOORD=1     # if more fragments this determine what to do (0-add cartesians for all but the alrgest
                 # fragments, 1-add longer distances, 2-add inverse power distances (1/R),3-add
                 # 1/R^6). in cases of 1, 2 and 3 the distances are generated using BSCALE*ASCALE. new
                 # coordinates are not used for generation of angles, torsions etc...
 PRIMRESET=0     # reset primitive internals (0/1 --->N/Y)
 TS=0
 IEIGVAL=0       # which eigenvalue should be made negative in TS==0 simulations
 HOMEADDRESS=os.getcwd()
 COMPADDRESS=os.getcwd()
 PULAYGUESS=zeros(6,float)  # this will be subtracted from stress tensor, makes sense only
                            # for full relaxations (RELAX=1)
 RIGIDSTRESS=0              # 0/1 --->actual/given stress
 RSTRESS=zeros(6,float)     # if RIGIDSTRESS=1: only for test purposes (replaces the actual stress by this line)
 SCALEBT=1.0
 SUBST=100
 SUBSTINDEX=10000
 CARTSUBST=0  # use Cartesian coordinates for substrate atoms
# beg ew
 IFCON=0	#forces full connectivity in the cell; IFCON=0 - proceed normally; IFCON=1 - ASCALE is increased until full connectivity is reached
# IFCON=2 - long-range bonds are added to the short-range bonds
 LRVDW=0
# end ew
 ATRADII=[]
 iteration=0
 
 POMASS=[]

 def __init__(self):
  cexist=0

  try:   
    ifile=open('INPDAT','r')
  except IOError:
    print "problem reading input file (INPDAT), I'll try to use defaults"
  else:
    for line in ifile.readlines():
     line=string.split(line,'=')
     if len(line)!=2:continue
     line[0]=string.strip(line[0])
     line[1]=string.strip(line[1])
     line[1]=string.split(line[1])
     line_long=[]
     if len(line[1])>1:
       line_long=line[1][1:]
     line[1]=line[1][0]
     if line[0]=='RELAX':
       self.RELAX=int(line[1])
     elif line[0]=='EQUIL':
       self.EQUIL=int(line[1])
     elif line[0]=='WINCREM':
       self.WINCREM=int(line[1])
     elif line[0]=='BT_MAXITER':
       self.BT_MAXITER=int(line[1])
       if self.BT_MAXITER<1: self.BT_MAXITER=100
     elif line[0]=='SHAKE_MAXITER':
       self.SHAKE_MAXITER=int(line[1])
       if self.SHAKE_MAXITER<1: self.SHAKE_MAXITER=100
     elif line[0]=='TEBEG':
       self.TEBEG=float(line[1])
     elif line[0]=='ANDERSEN_PROB':
       self.ANDERSEN_PROB=float(line[1])
       if self.ANDERSEN_PROB<0.: self.ANDERSEN_PROB=0.
       if self.ANDERSEN_PROB>1.: self.ANDERSEN_PROB=0.
     elif line[0]=='TILAMBDA':
       self.TILAMBDA=float(line[1])
       if self.TILAMBDA>1.:self.TILAMBDA=0.
       if self.TILAMBDA<0.:self.TILAMBDA=0.
     elif line[0]=='SHAKE_TOLER':
       self.SHAKE_TOLER=float(line[1])
       if self.SHAKE_TOLER<0.:self.SHAKE_TOLER=1e-7
     elif line[0]=='BT_TOLER':
       self.BT_TOLER=float(line[1])
       if self.BT_TOLER<0.:self.BT_TOLER=1e-7
     elif line[0]=='CART':
       self.CART=int(line[1])
     elif line[0]=='TORSIONS':
       self.TORSIONS=int(line[1])
       if abs(self.TORSIONS)>1:self.TORSIONS=1
     elif line[0]=='HESSIAN':
       self.HESSIAN=int(line[1])
       if self.HESSIAN>3 or self.HESSIAN<0:self.HESSIAN=1
     elif line[0]=='HUPDATE':
       self.HUPDATE=int(line[1])
     elif line[0]=='HREFRESH':
       self.HREFRESH=int(line[1])
     elif line[0]=='PRIMRESET':
       self.PRIMRESET=int(line[1])
     elif line[0]=='TS':
       self.TS=int(line[1])
     elif line[0]=='IEIGVAL':
       self.IEIGVAL=int(line[1])
       if self.IEIGVAL<0: self.IEIGVAL=0
     elif line[0]=='HUPDIM':
       self.HUPDIM=int(line[1])
     elif line[0]=='OPTENGINE':
       self.OPTENGINE=int(line[1])
     elif line[0]=='LINESEARCH':
       self.LINESEARCH=int(line[1])
     elif line[0]=='LINEMAX':
       self.LINEMAX=int(line[1])
     elif line[0]=='CONSTCONJ':
       self.CONSTCONJ=int(line[1])
     elif line[0]=='HSING':
       self.HSING=float(line[1])
     elif line[0]=='POTIM':
       self.POTIM=float(line[1])
     elif line[0]=='KMATAL':
       self.KMATAL=int(line[1])
       if self.KMATAL>1 or self.KMATAL<0: self.KMATAL=0
     elif line[0]=='NFREE':
       self.NFREE=int(line[1])
       if self.NFREE>10:self.NFREE==10
     elif line[0]=='NSW':
       self.NSW=int(line[1])
     elif line[0]=='STEPLIM':
       self.STEPLIM=float(line[1])
     elif line[0]=='GCRITER':
       self.GCRITER=float(line[1])
     elif line[0]=='ECRITER':
       self.ECRITER=float(line[1])
     elif line[0]=='SCRITER':
       self.SCRITER=float(line[1])
     elif line[0]=='ASCALE':
       self.ASCALE=float(line[1])
     elif line[0]=='BSCALE':
       self.BSCALE=float(line[1])
     elif line[0]=='CSCALE':
       cexist=1
       self.CSCALE=float(line[1])
     elif line[0]=='ANGLECRIT':
       self.ANGLECRIT=int(line[1])
     elif line[0]=='TORSIONCRIT':
       self.TORSIONCRIT=int(line[1])
     elif line[0]=='SCALEBT':
       self.SCALEBT=float(line[1])
     #elif line[0]=='VDWRAD':
     #  self.VDWRAD=float(line[1])
     elif line[0]=='SUBST':
       self.SUBST=int(line[1])
     elif line[0]=='SUBSTINDEX':
       self.SUBSTINDEX=int(line[1])
     elif line[0]=='CARTSUBST':
       self.CARTSUBST=int(line[1])
       if abs(self.CARTSUBST)>1: self.CARTSUBST=0
# beg ew
     elif line[0]=='IFCON':
       self.IFCON=int(line[1])
     elif line[0]=='LRVDW':
       self.LRVDW=int(line[1])
# end ew
     elif line[0]=='FRAGCOORD':
       self.FRAGCOORD=int(line[1])
     elif line[0]=='HOMEADDRESS':
       self.HOMEADDRESS=line[1]
     elif line[0]=='COMPADDRESS':
       self.COMPADDRESS=line[1]
     elif line[0]=='PULAYGUESS':
       self.PULAYGUESS[0]=float(line[1])
       for i in range(1,6):
         self.PULAYGUESS[i]=float(line_long[i-1])
     elif line[0]=='RIGIDSTRESS':
       self.RIGIDSTRESS=1
       self.RSTRESS[0]=float(line[1])
       for i in range(1,6):
         self.RSTRESS[i]=float(line_long[i-1])
     elif line[0]=='ATRADII':
       self.ATRADII.append(float(line[1]))
       for i in range(len(line_long)):
	 try:
	   ll=float(line_long[i])
         except ValueError:
	   break
	 else:
	   self.ATRADII.append(ll)
     elif line[0]=='POMASS':
       self.POMASS.append(float(line[1]))
       for i in range(len(line_long)):
         try:
           ll=float(line_long[i])
         except ValueError:
           break
         else:
           self.POMASS.append(ll)


    if self.EQUIL>=self.NSW: self.EQUIL=0	   
    if self.HUPDIM<2:self.HUPDIM=2
    if (self.TS==1 and (self.HUPDATE==1 or self.HUPDATE==6)):
      self.HUPDATE=5
    if cexist==0:self.CSCALE=self.ASCALE
    if self.CART==1:self.ASCALE=0.0
# beg ew
#    if self.IFCON!=0:
#       print "IFCON forces full connectivity in the unit cell - parameter default values SUBST=100 and BSCALE=2.0 are used"
#       self.SUBST=100
#      # self.BSCALE=2.0
# end ew
    ifile.close()

