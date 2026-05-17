import os
import subprocess
import multiprocessing as mp
from itertools import product

diff = 0.1
start = 0.125

# -----------------------
# Configuration
# -----------------------
tgt = [1000180400]
ke = [round(start + diff * i, 4) for i in range(5)]
probe = [2212, 211]
rw = ["hN2018", "INCL++", "Geant4"]
N_WORKERS = 40

# -----------------------
# Helpers
# -----------------------
def run_command(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def run_job(args):
    particle, target, k, r, job_id = args

    run_command([
        "root", "-l", "-b", "-q",
        f"plotTotNuclRW.C({particle},{k},{target},1,\"{r}\")"
    ])
    run_command([
        "root", "-l", "-b", "-q", f"plotVisibleERW.C({particle},{k},{target},1,\"{r}\")"
    ])
    run_command([
        "root", "-l", "-b", "-q", f"plotVisibleEwAllRW.C({particle},{k},{target},1,\"{r}\")"
    ])
    run_command([
        "root", "-l", "-b", "-q",
        f"plotTotNuclwAllRW.C({particle},{k},{target},1,\"{r}\")"
    ])
# -----------------------
# Main
# -----------------------
if __name__ == "__main__":
    jobs = []
    job_id = 0
    for particle, target, k, r in product(probe, tgt, ke, rw):
        jobs.append((particle, target, k, r, job_id))
        job_id += 1

    with mp.Pool(processes=N_WORKERS) as pool:
        pool.map(run_job, jobs)