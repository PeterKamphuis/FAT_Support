#!/usr/local/bin/ python3


# This program reads the system configurations of your machine
# determines the threads available divides up the amount of threads and runs FAT at 95% capacity by dividing up the input list but it keeps using 4 cores per galaxy.
FATInput='/Users/Peter/WALLABY/LVHIS-26_v2.3/Input.config'
#FATInput='/data/users/kamphuis/Artificial/FAT_INPUT.config'


import os
import numpy as np
import multiprocessing
import subprocess
from collections import OrderedDict
import time
import copy
import sys
# Code for creating a proper dictionary instead of this python unordered nonsense.
# This also includes an insert class.
class MyOrderedDict(OrderedDict):
     def insert(self,existing_key,new_key,key_value):
        done = False
        if new_key in self:
            self[new_key]=key_value
            done = True
        else:
            new_orderded_dict=self.__class__()
            for key, value in self.items():
                new_orderded_dict[key]=value
                if key==existing_key:
                    new_orderded_dict[new_key]=key_value
                    done = True
            if not done:
                new_orderded_dict[new_key]=key_value
                done = True
                print("----!!!!!!!! YOUR new key was appended at the end as you provided a non-existing key to add it after!!!!!!---------")
            self.clear()
            self.update(new_orderded_dict)

        if not done:
            print("----!!!!!!!!We were unable to add your key!!!!!!---------")



# first get the amount of available threads
ncores = multiprocessing.cpu_count()
#ncores = 40

start_time = time.time()

tmp = open(FATInput,'r')
Template_in = MyOrderedDict({})
unarranged = tmp.readlines()
# Separate the keyword names
for tmp in unarranged:
    # python is really annoying with needing endlines. Let's strip them here and add them when writing
    Template_in[tmp.split('=',1)[0].strip().upper()]=tmp.rstrip()
    if  tmp.split('=')[0] == 'catalogue':
        filefound=  tmp.split('=')[1].strip()
    if  tmp.split('=')[0] == 'startgalaxy':
        start= int(tmp.split('=')[1].strip())
    if  tmp.split('=')[0] == 'endgalaxy':
        end= int(tmp.split('=')[1].strip())

#replace space
filefound = "\ ".join(filefound.split())
print(filefound)
if end == -1:
    proc = subprocess.Popen(["wc", filefound], stdout=subprocess.PIPE)
    (out, err) = proc.communicate()
    end = int(out.split()[0]) -1

# So we need to fit this many galaxies
no_gal=end-start
# and we use four cores per galaxy so we can run
no_procs = int(ncores/4)
# if this is more than the number of galaxies we use  1 run per galaxy
if no_procs > no_gal: no_procs = no_gal
#we need so many galaxies per run
gal_step_1= np.ceil(no_gal/no_procs)
gal_step_2= np.floor(no_gal/no_procs)
count=0
if gal_step_1 != gal_step_2:
    times1 = 1
    times2 = 1
    total=times1*gal_step_1+times2*gal_step_2
    while no_gal-total > gal_step_1:
        times1 += 1
        times2 += 1
        total=times1*gal_step_1+times2*gal_step_2
        print(total)
        count += 1
#now write input files and spawn the processes








print(filefound,start,end)
outfile=FATInput.split("/")[-1].split(".config")[0]
outdir="/".join(FATInput.split("/")[0:-1])
counts=np.zeros(no_procs)
proc= dict()
print(outfile,outdir)
filegh= Template_in["OUTPUTCATALOGUE"].split("=")[1]
for i in range(no_procs):
    Template=copy.deepcopy(Template_in)
    Template["STARTGALAXY"]="startgalaxy = {}".format(int(start))
    if i <= count:
        Template["ENDGALAXY"]="endgalaxy = {}".format(int(start+gal_step_1-1))
        start +=gal_step_1
    else:
        Template["ENDGALAXY"]="endgalaxy = {}".format(int(start+gal_step_2-1))
        start +=gal_step_2
    if i == no_procs-1:
        Template["ENDGALAXY"]="endgalaxy = {}".format(-1)
    Template["OUTPUTCATALOGUE"]=Template["OUTPUTCATALOGUE"].split(".txt")[0]+"_auto{:d}".format(i)+".txt"
    tri = open(outdir+"/"+outfile+'_auto{:d}.config'.format(i), 'w')
    tri.writelines([Template[key]+"\n" for key in Template])
    tri.close()
    proc[i] = subprocess.Popen(["idl", "-e fat.pro,configuration_file='"+outdir+"/"+outfile+'_auto{:d}.config'.format(i)+"'"], stdout=subprocess.PIPE)
    #test = subprocess.Popen("cp "+filegh.split(".txt")[0]+"_cauto{:d}\ copy".format(i)+".txt "+filegh.split(".txt")[0]+"_auto{:d}".format(i)+".txt " , stdout=subprocess.PIPE, shell=True)
print("yes we finished starting")
# We want to wait for all the fitting
fullout= ['                                                        Name         AC1         AC2\n']
for i in range(no_procs):
    proc[i].wait()
    print("Process {} finished".format(i))
    counter = 0
#    print(filegh.split(".txt")[0]+"_auto{:d}\ copy".format(i)+".txt "+filegh.split(".txt")[0]+"_auto{:d}".format(i)+".txt ")

    with open(filegh.split(".txt")[0]+"_auto{:d}".format(i)+".txt",'r') as f:
        for line in f:
            if (counter != 0) & (not line.isspace()) :
                fullout.append(line)
                print("we doing this")
            counter += 1
    #cleanup
    subprocess
tri = open(Template_in["OUTPUTCATALOGUE"].split("=")[1] , 'w')
tri.writelines([f for f in fullout])
tri.close()
#And the we want to read the results to combine them into the single output classes
print(fullout)
print("It took {} seconds to run this".format(time.time() - start_time))
#need to clean up
#proc = subprocess.Popen("rm -f "+outdir+"/"+outfile+'_auto*.config', stdout=subprocess.PIPE, shell=True)
proc = subprocess.Popen("rm -f "+filegh.split(".txt")[0]+"_auto*"+".txt", stdout=subprocess.PIPE, shell=True)


#for i in range(5):
#    for func in [mapcount, simplecount, bufcount, opcount]:
#        start_time = time.time()
#        assert func("big_file.txt") == 1209138
#        counts[func].append(time.time() - start_time)
