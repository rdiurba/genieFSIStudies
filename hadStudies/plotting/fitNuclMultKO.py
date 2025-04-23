import os
#import numpy as np
#os.system("source ../setupNuSystematics.sh")
tgt=[1000060120,1000260560,1000822080,1000180400]
ke=[0.125,0.25,0.5]
num=1000000
probe=[2212]
#seed=np.random.randint(0,10000)
seed=0
os.system("rm nucleonKO*.txt")
with open("nucleonKODifference.txt", "w") as f:
    f.write("probe,ke,target,meanhA,meanErrhA,stdhA,stdErrhA,meannhN,meanErrhN,stdhN,stdErrhN,meanINCL,meanErrINCL,stdINCL,stdErrINCL,meanGeant4,meanErrGeant4,stdGeant4,stdErrGeant4")
    f.write("\n")
with open("nucleonKOSum.txt", "w") as f:
    f.write("probe,ke,target,consthA,constErrhA,slopehA,slopeErrhA,consthN,constErrhN,slopehN,slopeErrhN,constINCL,constErrINCL,slopeINCL,slopeErrINCL,constGeant4,constErrGeant4,slopeGeant4,slopeErrGeant4")
    f.write("\n")
with open("nucleonKODifferenceThresh.txt", "w") as f:
    f.write("probe,ke,target,meanhA,meanErrhA,stdhA,stdErrhA,meannhN,meanErrhN,stdhN,stdErrhN,meanINCL,meanErrINCL,stdINCL,stdErrINCL,meanGeant4,meanErrGeant4,stdGeant4,stdErrGeant4")
    f.write("\n")
with open("nucleonKOSumThresh.txt", "w") as f:
    f.write("probe,ke,target,consthA,constErrhA,slopehA,slopeErrhA,consthN,constErrhN,slopehN,slopeErrhN,constINCL,constErrINCL,slopeINCL,slopeErrINCL,constGeant4,constErrGeant4,slopeGeant4,slopeErrGeant4")
    f.write("\n")





for particle in probe:
    for t in tgt:
        for k in ke:
            seed=seed+1
            os.system("root -l -b -q 'plotKODiff.C("+str(particle)+","+str(k)+","+str(t)+")'")
            os.system("root -l -b -q 'plotKOSum.C("+str(particle)+","+str(k)+","+str(t)+")'")
            os.system("root -l -b -q 'plotKODiffThresh.C("+str(particle)+","+str(k)+","+str(t)+")'")
            os.system("root -l -b -q 'plotKOSumThresh.C("+str(particle)+","+str(k)+","+str(t)+")'")




