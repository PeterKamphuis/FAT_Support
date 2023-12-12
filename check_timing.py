#!/usr/local/bin/ python3
# -*- coding: utf-8 -*-
import numpy as np

def main():
    file= '/home/peter/FAT_Main/HPC_Results/Timing_Result.txt'
    galaxies  = read_galaxies_in_timing_result(file)

    occurence = [[x, galaxies.count(x)] for x in set(galaxies)]
    for pew in occurence:
        if pew[1] > 1:
            print(f'The galaxy {pew[0]} has been fitted {pew[1]} times')

def read_galaxies_in_timing_result(filename):

    with open(filename) as file:
        all = file.readlines()
    galaxies = []
    for line in all:
        tmp = line.split()
        if tmp[0] == 'The' and tmp[1] == 'galaxy':
            galaxies.append(tmp[4].split('/')[-2])
    return galaxies
if __name__ == '__main__':
    main()
