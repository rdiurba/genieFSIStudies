import os
import subprocess
import multiprocessing as mp
from itertools import product

diff = 0.05
start = 0.075

# -----------------------
# Configuration
# -----------------------
tgt = [1000060120, 1000180400, 1000080160, 1000260560]
ke = [round(start + diff * i, 4) for i in range(40)]
num = 1000000
probe = [2212, 211, 2112, -211,111]
BASE_SEED = 0
N_WORKERS = 40
os.system("rm nucleonKO*.txt")
os.system("rm pionAbs*.txt")
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

with open("pionAbsDifferenceCompoundThresh.txt", "w") as f:
    f.write("probe,ke,target,meanhA,meanErrhA,stdhA,stdErrhA,meannhN,meanErrhN,stdhN,stdErrhN,meanINCL,meanErrINCL,stdINCL,stdErrINCL,meanGeant4,meanErrGeant4,stdGeant4,stdErrGeant4")
    f.write("\n")
with open("pionAbsSumCompoundThresh.txt", "w") as f:
    f.write("probe,ke,target,meanhA,meanErrhA,stdhA,stdErrhA,meannhN,meanErrhN,stdhN,stdErrhN,meanINCL,meanErrINCL,stdINCL,stdErrINCL,meanGeant4,meanErrGeant4,stdGeant4,stdErrGeant4")

with open("nucleonKODifferenceCompoundThresh.txt", "w") as f:
    f.write("probe,ke,target,meanhA,meanErrhA,stdhA,stdErrhA,meannhN,meanErrhN,stdhN,stdErrhN,meanINCL,meanErrINCL,stdINCL,stdErrINCL,meanGeant4,meanErrGeant4,stdGeant4,stdErrGeant4")
    f.write("\n")
with open("nucleonKOSumCompoundThresh.txt", "w") as f:
    f.write("probe,ke,target,consthA,constErrhA,slopehA,slopeErrhA,consthN,constErrhN,slopehN,slopeErrhN,constINCL,constErrINCL,slopeINCL,slopeErrINCL,constGeant4,constErrGeant4,slopeGeant4,slopeErrGeant4")
    f.write("\n")

"""
def run_root(cmd):
    os.system(cmd)

def main():
    diff = 0.05
    start = 0.075

    tgt = [1000060120, 1000180400, 1000080160, 1000260560]
    ke = [round(start + diff * i, 4) for i in range(40)]
    probe = [2212, 211]

    seed = 0

    with ProcessPoolExecutor(max_workers=5) as executor:
        for particle in probe:
            for t in tgt:
                for k in ke:
                    seed += 1

                    futures = [
                        executor.submit(run_root, f"root -l -b -q 'plotVisibleE.C({particle},{k},{t},1,0)' > /dev/null 2>&1"),
                         executor.submit(run_root, f"root -l -b -q 'plotVisibleEStack.C({particle},{k},{t},1,0)' > /dev/null 2>&1"),

                        executor.submit(run_root, f"root -l -b -q 'plotTotalNucl.C({particle},{k},{t},1,0)' > /dev/null 2>&1"),
                        executor.submit(run_root, f"root -l -b -q 'plotDiffNucl.C({particle},{k},{t},1,0)' > /dev/null 2>&1"),
                    ]

                    wait(futures)

if __name__ == "__main__":
    main()
"""
# -----------------------
# Helpers
# -----------------------


def run_command(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def run_job(args):
    particle, target, k, job_id = args

    run_command([
        "root", "-l", "-b", "-q",
        f"plotVisibleE.C({particle},{k},{target},1,0)"
    ])
    run_command([
        "root", "-l", "-b", "-q",
        f"plotVisibleEStack.C({particle},{k},{target},1,0)"
    ])
    run_command([
        "root", "-l", "-b", "-q",
        f"plotTotalNucl.C({particle},{k},{target},1,0)"
    ])
    run_command([
        "root", "-l", "-b", "-q",
        f"plotDiffNucl.C({particle},{k},{target},1,0)"
    ])
    run_command([
        "root", "-l", "-b", "-q",
        f"plotVisibleE.C({particle},{k},{target},1,1)"
    ])
    run_command([
        "root", "-l", "-b", "-q",
        f"plotVisibleEStack.C({particle},{k},{target},1,1)"
    ])
    run_command([
        "root", "-l", "-b", "-q",
        f"plotTotalNucl.C({particle},{k},{target},1,1)"
    ])
    run_command([
        "root", "-l", "-b", "-q",
        f"plotDiffNucl.C({particle},{k},{target},1,1)"
    ])
# -----------------------
# Main
# -----------------------
if __name__ == "__main__":

    jobs = []
    job_id = 0
    for particle, target, k in product(probe, tgt, ke):
        jobs.append((particle, target, k, job_id))
        job_id += 1



    with mp.Pool(processes=N_WORKERS) as pool:
        pool.map(run_job, jobs)
