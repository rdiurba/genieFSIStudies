import os
#import numpy as np
#os.system("source ../setupNuSystematics.sh")
tgt=[1000060120,1000260560,1000822080,1000180400]
ke=[0.125,0.25,0.5]
num=1000000
probe=[211]
#seed=np.random.randint(0,10000)
seed=0
os.system("rm pionAbs*.txt")
with open("pionAbsDifference.txt", "w") as f:
    f.write("probe,ke,target,meanhA,meanErrhA,stdhA,stdErrhA,meannhN,meanErrhN,stdhN,stdErrhN,meanINCL,meanErrINCL,stdINCL,stdErrINCL,meanGeant4,meanErrGeant4,stdGeant4,stdErrGeant4")
    f.write("\n")
with open("pionAbsSum.txt", "w") as f:
    f.write("probe,ke,target,meanhA,meanErrhA,stdhA,stdErrhA,meannhN,meanErrhN,stdhN,stdErrhN,meanINCL,meanErrINCL,stdINCL,stdErrINCL,meanGeant4,meanErrGeant4,stdGeant4,stdErrGeant4")
    f.write("\n")

with open("pionAbsDifferenceThresh.txt", "w") as f:
    f.write("probe,ke,target,meanhA,meanErrhA,stdhA,stdErrhA,meannhN,meanErrhN,stdhN,stdErrhN,meanINCL,meanErrINCL,stdINCL,stdErrINCL,meanGeant4,meanErrGeant4,stdGeant4,stdErrGeant4")
    f.write("\n")
with open("pionAbsSumThresh.txt", "w") as f:
    f.write("probe,ke,target,meanhA,meanErrhA,stdhA,stdErrhA,meannhN,meanErrhN,stdhN,stdErrhN,meanINCL,meanErrINCL,stdINCL,stdErrINCL,meanGeant4,meanErrGeant4,stdGeant4,stdErrGeant4")
    f.write("\n")
for particle in probe:
    for t in tgt:
        for k in ke:
            seed=seed+1
            os.system("root -l -b -q 'plotAbsDiff.C("+str(particle)+","+str(k)+","+str(t)+")'")
            os.system("root -l -b -q 'plotAbsSum.C("+str(particle)+","+str(k)+","+str(t)+")'")
            os.system("root -l -b -q 'plotAbsDiffThresh.C("+str(particle)+","+str(k)+","+str(t)+")'")
            os.system("root -l -b -q 'plotAbsSumThresh.C("+str(particle)+","+str(k)+","+str(t)+")'")





