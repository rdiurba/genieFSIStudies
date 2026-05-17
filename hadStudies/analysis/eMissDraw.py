import os
#import numpy as np
#os.system("source ../setupNuSystematics.sh")
tgt=[1000180400]
ke = [round(0.075 + 0.1 * i, 4) for i in range(10)]
models=["hA2018","HINCL"]
num=1000000
probe=[211,2212]
#seed=np.random.randint(0,10000)
seed=0


for particle in probe:
    for t in tgt:
        for model in models:
            for k in ke:
                seed=seed+1
    
                os.system(f'root -l -b -q \'plotVisibleE.C({particle},{k},{t},1,1,"{model}")\'')
