#! /usr/bin/python
import os

# How to merge files into a bigger one. First one is the target
#IOV_list= ['2018P8','2018A','2018B','2018C','2018D1','2018D2']
IOV_list_of_lists = [
#    ['2016BCDEF','2016BCD','2016EF'],
#    ['2017BCDEF','2017B','2017C','2017D','2017E','2017F'],
#    ['2018ABCD','2018A1','2018A2','2018B','2018C',
#     '2018D1','2018D2','2018D3','2018D4'],
#    ['Run2','2016BCDEF','2016FGH','2017BCDEF','2018ABCD']
#    ['2022CD','2022C','2022D'],
#    ['2022CDE','2022C','2022D','2022E'],
#    ['2022FG','2022F','2022G'],
#    ['2023Cv4D','2023Cv4','2023D'],
#    ['Run3','2022C','2022D','2022E','2022F','2022G',
#     '2023Cv123','2023Cv4','2023D']
##    ['Run3Summer23', '2023Cv123', '2023Cv4','2023D'],
    ['2024only', '2024B', '2024C', '2024D']
#    ['Run3', '2023Cv123', '2023Cv4','2023D']
    ]
MC_list_of_lists = [
#    ['Run2P8','2016P8','2016APVP8','2017P8','2018P8'],
#    ['Run2QCD','2016QCD','2016QCDAPV','2017QCD','2018QCD'],
#    ['Run3P8','2022P8','2022EEP8']
##    ['Run3Summer23', '2023P8', '2023P8-BPix'] # !! does not make sense to merge these together? - double check!!
    ]

version = 'w17'

os.system("ls rootfiles/GamHistosFill_data_*_"+version+".root")
for IOV_list in IOV_list_of_lists:
    command = "hadd "
    for iov in IOV_list:
        command = command + "rootfiles/GamHistosFill_data_"+iov+"_"+version+".root "
    print("\""+command+"\"...")
    os.system(command)

os.system("ls rootfiles/GamHistosFill_mc_*_"+version+".root")
for MC_list in MC_list_of_lists:
    command = "hadd "
    for mc in MC_list:
        command = command + "rootfiles/GamHistosFill_mc_"+mc+"_"+version+".root "
    print("\""+command+"\"...")
    os.system(command)

#for iov in IOV_list:
#    print "Process GamHistFill.C for IOV "+iov
#    os.system("ls -ltrh files/GamHistosFill_mc_"+iov+".root")
#    os.system("ls -ltrh files/GamHistosFill_data_"+iov+".root")
#    os.system("ls -ltrh log_"+iov+"_"+version+".txt")
#    os.system("root -l -b -q 'mk_GamHistosFill.C(\""+iov+"\")' > log_"+iov+"_"+version+".txt &")
#    os.system("fs flush")
#    wait()
#    time.sleep(sleep_time)
