#!/usr/bin/env python3.6

import numpy as np
from ase.io import read
from ase.io.vasp import write_vasp

def main():

  at = read('ACCEPTED.xyz',index='-1')
  write_vasp('POSCAR.last',at,label='System',vasp5=True)

if __name__ == '__main__':
  main()

