#!/usr/local/bin/ python3
# -*- coding: utf-8 -*-


import os


def main():
  '''Copy a set of directories in its entirety from Bochum to home'''
  ip = '192.18.0.27'
  dir = '/home/peter/FAT_Main/Test_Analysis/GDL_fits/'

  #get the directories to copy.
  with open(f'{dir}Galaxies.txt') as file:
      tmp = file.readlines()
  dir_requested = []
  for line in tmp:
    line_split = line.split()
    directory_requested = line_split[0].strip()
    os.system(f'scp -r {ip}:./FAT/Database/{directory_requested}/ {dir}')


if __name__ == '__main__':
    main()
