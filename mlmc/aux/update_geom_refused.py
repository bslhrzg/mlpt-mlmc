#!/usr/bin/env python3.6

from ase.io import read, write

def main():
  at = read('ACCEPTED.xyz',index='-1',format='extxyz')
  write('ACCEPTED.xyz',at,format='extxyz',append=True)

if __name__ == '__main__':
  main() 
