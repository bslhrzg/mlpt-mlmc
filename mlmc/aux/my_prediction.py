#!/usr/bin/env python3.6

from mlpt import MLPT
from sys import argv

def main():
  e0=MLPT().prediction('POSCAR',float(argv[1]))
  print('{0}'.format(e0))

if __name__ == '__main__':
  main() 
