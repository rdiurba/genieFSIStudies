import os
import numpy as np
path=""
tgt=[1000060120,1000260560,1000822080]
ke=[0.125,0.25,0.5]
num=1000000
probe=[211,-211]
seed=np.random.randint(0,10000)
seed=400
for particle in probe:
    probeStr="piPlus"
    if( particle==2212): probeStr="protonPlus"
    if (particle==-2212): probeStr="protonMinus"
    if(particle==-211): probeStr="piMinus"
    if (particle==111): probeStr="pi0"
    if (particle==2112): probeStr="neutron"
    for t in tgt:
        for k in ke:
            seed=seed+1
            ghepFile=path+""+str(probeStr)+"_"+str(t)+"_"+str(k)+"GeV_HG4BertCasc_1M"
            gstFile=path+""+str(probeStr)+"_"+str(t)+"_"+str(k)+"GeV_HG4BertCasc_1M.ginuke.root"
            os.system("gevgen_hadron --seed "+str(seed)+" -n "+str(num)+" -p "+str(particle)+" -t "+str(t)+" -k "+str(k)+" -m HG4BertCasc -o "+str(ghepFile))
            os.system("gntpc -f ginuke -i "+str(ghepFile)+".0.ghep.root -o "+str(gstFile))
            seed=seed+1
            """
            ghepFile="/dune/data/users/rdiurba/genieFSIStudies/"+str(probeStr)+"_"+str(t)+"_"+str(k)+"GeV_hA2018_1M"
            gstFile="/dune/data/users/rdiurba/genieFSIStudies/"+str(probeStr)+"_"+str(t)+"_"+str(k)+"GeV_hA018_1M.ginuke.root"
            os.system("gevgen_hadron -n "+str(num)+" -p "+str(particle)+" -t "+str(t)+" -k "+str(k)+" -m hA2018 -o "+str(ghepFile)+" --seed "+str(seed))
            os.system("gntpc -f ginuke -i "+str(ghepFile)+".0.ghep.root -o "+str(gstFile))
            """


