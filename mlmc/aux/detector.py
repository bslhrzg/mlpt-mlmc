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

e0=float(sys.argv[1])
e1=float(sys.argv[2])
temperature=float(sys.argv[3])

pC=physconstants.physConstants()

beta=1./(pC.bolkEV*temperature)


verdict=metropolis_criter(e1,e0,beta)
print(verdict)

