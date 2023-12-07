#!/usr/bin/env python3.6

from ase.io import read, write

def main():
  at = read('vasprun.xml',index=':',format='vasp-xml')
  write('ACCEPTED.xyz',at[0],format='extxyz',append=True)

if __name__ == '__main__':
  main() 
