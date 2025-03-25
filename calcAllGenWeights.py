#! /usr/bin/python3
import os

#version - if wanted
ver='w46'

# PTG and HT binned samples
pthtbin_list= ['summer2024P8_PTG10to100-HT40to100', 'summer2024P8_PTG10to100-HT100to200', 
            'summer2024P8_PTG10to100-HT200to400', 'summer2024P8_PTG10to100-HT400to600', 
            'summer2024P8_PTG10to100-HT600to1000', 'summer2024P8_PTG10to100-HT1000toInf', 
            'summer2024P8_PTG100to200-HT40to200', 'summer2024P8_PTG100to200-HT200to400', 
            'summer2024P8_PTG100to200-HT400to600', 'summer2024P8_PTG100to200-HT600to1000', 
            'summer2024P8_PTG100to200-HT1000toInf', 'summer2024P8_PTG200toInf-HT40to400', 
            'summer2024P8_PTG200toInf-HT400to600', 'summer2024P8_PTG200toInf-HT600to1000', 
            'summer2024P8_PTG200toInf-HT1000toInf']

#os.system("rm *.so *.d *.pcm")
#os.system("root -l -b -q mk_CondFormats.C") #could add any prerequisites here
os.system("root -l -b -q mk_GenWeightLibrary.C")

for pthtbin in pthtbin_list:
    print("Process CalcGenWeight.C+g for bin: "+pthtbin)
    #os.system("ls -ltrh files/genweight_"+pthtbin+".txt")
    #os.system("ls -ltrh log_CalcGenWeight_"+pthtbin+"_"+ver+".txt")
    os.system("root -l -b -q 'mk_CalcGenWeight.C(\""+pthtbin+"\")' > log/log_CalcGenWeight_"+pthtbin+"_"+ver+".txt &")
    #os.system("root -l -b -q 'mk_CalcGenWeight.C(\""+pthtbin+"\")' &")
#    os.system("fs flush")
#    wait()
#    time.sleep(sleep_time)
