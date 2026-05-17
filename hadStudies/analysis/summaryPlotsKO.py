# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 15:13:03 2025

@author: Richie-LHEP
"""

import csv
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from iminuit import Minuit
from iminuit.cost import LeastSquares

# ── Global tight-layout settings ──────────────────────────────
plt.rcParams['figure.autolayout'] = True        # equivalent to tight_layout() on every fig
plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.bbox'] = 'tight'          # bbox_inches='tight' on every savefig
plt.rcParams['savefig.pad_inches'] = 0.05       # minimal padding on save
# ──────────────────────────────────────────────────────────────
minKE=0.125;
maxKE=2.00;
pionCarbonhA=[]
ke=[]
probe=211
tgt=1000060120
index=0
typ="Difference" # "Sum"

def linear(x, slope, intercept):
    return intercept + slope * x

def expo(x, slope, intercept,floor,exp):
    return intercept*np.exp((x**exp)*slope)+floor
def hAIntranuke2018Calculator(probe,A,ke,typ):
    val=-99
    ke=ke*1000.0
    c1 = 0.041 + ke * 0.0001525;
    c2 = -0.003444 - ke * 0.00002324;
    
    c3 = 0.064 - ke * 0.000015;
    gam_ns = c1 * np.exp(c2*A) + c3;
    if(gam_ns<0.002): gam_ns = 0.002;
    if (typ=="Sum"): val=gam_ns
    Z=6
    if (A==40): Z=18 
    if (A==56): Z=26
    if (A==208): Z=82 

    if (typ=="Difference"):
        if (A - Z > Z):
            nd0 = 135.227 * np.exp(-7.124 * (A - Z)/A) - 2.762
        else:
            nd0 = -135.227 * np.exp(-7.124 * Z / A) + 4.914
        val=nd0

    return val
def hAIntranuke2018Std(probe,A,ke,typ):
    val=-99
    ke=ke*1000.0
    c1 = 0.041 + ke * 0.0001525;
    c2 = -0.003444 - ke * 0.00002324;
    
    c3 = 0.064 - ke * 0.000015;
    gam_ns = c1 * np.exp(c2*A) + c3;
    if(gam_ns<0.002): gam_ns = 0.002;
    if (typ=="Sum"): val=gam_ns
    Z=6
    if (A==40): Z=18 
    if (A==56): Z=26
    if (A==208): Z=82 
    RemnA = A    # approximate remnant nucleons
    if (typ=="Difference"):
        val=2.034 + RemnA * 0.007846
    return val
def linearPlotSigmaEnergy(probe,tgt,index,typ,thresh=0):
    threshStr=""
    if (thresh==1):
        threshStr="Thresh"
    if (thresh==2):
        threshStr="CompoundThresh"
    pionCarbonhA=[]
    pionCarbonhAlt=[]
    pionCarbonhAErr=[]
    pionCarbonhAltErr=[]
    ke=[]
    A=int(str(tgt)[6:9])
    with open('nucleonKO'+str(typ)+str(threshStr)+'.txt', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            if(row[0]=="probe"): continue
            if(row[0]!=str(probe)): continue
            if(row[2]!=str(tgt)): continue
            if (float(row[1])<minKE or float(row[1])>maxKE): continue;

            pionCarbonhA.append(float(row[5]))
            pionCarbonhAErr.append(float(row[6]))
            pionCarbonhAlt.append(float(row[5 + 4 * index]))
            pionCarbonhAltErr.append(float(row[6 + 4 * index]))

            ke.append(float(row[1]))
    indexStr="INCL++"
    if (index==1): indexStr="hN2018"
    if (index==3): indexStr="Geant4"
    plt.plot(ke,pionCarbonhA,"o",label="Mean Value Measured for hA2018")
    plt.plot(ke,pionCarbonhAlt,"o",label="Mean Value Measured for "+str(indexStr))

    plt.xlabel("Kinetic Energy [GeV]")
    probeLabel=r"$p$"

    plt.title("Proton Knockout Final State for "+str(probeLabel)+" for A="+str(int(str(int(tgt/10-1E8))[-3:]))+" Target")
    ylabel="Std. of "+typ+" of Nucleons"
    if (thresh==1): ylabel="Std. of "+typ+" of Nucleons above Threshold"
    cost_hA = LeastSquares(ke, pionCarbonhA, pionCarbonhAErr, linear)
    m_hA = Minuit(cost_hA, slope=1.0, intercept=0.0)
    m_hA.migrad()   # find minimum
    m_hA.hesse()    # compute uncertainties

    slope      = m_hA.values["slope"]
    slopeErr   = m_hA.errors["slope"]
    yIntercept    = m_hA.values["intercept"]
    yInterceptErr = m_hA.errors["intercept"]

    cost_alt = LeastSquares(ke, pionCarbonhAlt, pionCarbonhAltErr, linear)
    m_alt = Minuit(cost_alt, slope=1.0, intercept=0.0)
    m_alt.migrad()
    m_alt.hesse()

    slopeAlt      = m_alt.values["slope"]
    slopeAltErr   = m_alt.errors["slope"]
    yInterceptAlt    = m_alt.values["intercept"]
    yInterceptAltErr = m_alt.errors["intercept"]
    fittedVal=[linear(k,slope,yIntercept) for k in ke]
    fittedValAlt=[linear(k,slopeAlt,yInterceptAlt) for k in ke]
    simVal=[hAIntranuke2018Calculator(probe,A,k,typ) for k in ke]
    plt.plot(ke,fittedVal,"-", label="Fitted Slope for hA2018")
    plt.plot(ke,fittedValAlt,":", label="Fitted Slope for "+str(indexStr))

    #plt.plot(ke,simVal,"--", label="hA2018 Underlying Model")
    plt.legend(loc="lower right")
    plt.ylim(0,5)


    plt.ylabel(ylabel)
    plt.savefig("../plotting/linearPlotKEValueStd_"+str(probe)+"_"+str(tgt)+"_"+str(typ)+"_"+str(threshStr)+".pdf")

    plt.show()
    plt.close()
    return {
        "category": "Difference",
        "quantity": "Std",
        "hA2018_slope": slope,
        "hA2018_intercept": yIntercept,
        "alt_slope": slopeAlt,
        "alt_intercept": yInterceptAlt
    }
def linearPlotGammaEnergy(probe,tgt,index,typ,thresh=0):
    threshStr=""
    if (thresh==1):
        threshStr="Thresh"
    if (thresh==2):
        threshStr="CompoundThresh"
    pionCarbonhA      = []
    pionCarbonhAlt    = []
    pionCarbonhAErr   = []
    pionCarbonhAltErr = []
    ke=[]
    A=int(str(tgt)[6:9])
    with open('nucleonKO'+str(typ)+str(threshStr)+'.txt', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            if(row[0]=="probe"): continue
            if(row[0]!=str(probe)): continue
            if(row[2]!=str(tgt)): continue
            if (float(row[1])<minKE or float(row[1])>maxKE): continue;

            pionCarbonhA.append(float(row[5]))
            pionCarbonhAErr.append(float(row[6]))
            pionCarbonhAlt.append(float(row[5 + 4 * index]))
            pionCarbonhAltErr.append(float(row[6 + 4 * index]))

            ke.append(float(row[1]))
    indexStr="INCL++"
    if (index==1): indexStr="hN2018"
    if (index==3): indexStr="Geant4"
    plt.plot(ke,pionCarbonhA,"o",label="Mean Value Measured for hA2018")
    plt.plot(ke,pionCarbonhAlt,"o",label="Mean Value Measured for "+str(indexStr))

    plt.xlabel("Kinetic Energy [GeV]")
    probeLabel=r"$p$"

    plt.title("Proton Knockout Final State for "+str(probeLabel)+" for A="+str(int(str(int(tgt/10-1E8))[-3:]))+" Target")
    ylabel="Decay Constant of Nucleon Mult."
    if (thresh==1): ylabel="Decay Constant of Nucleon Mult. above Threshold"

    cost_hA = LeastSquares(ke, pionCarbonhA, pionCarbonhAErr, expo)
    m_hA = Minuit(cost_hA, slope=-1.0, intercept=5, floor=0,exp=1)
    m_hA.limits["intercept"] = (0, 10)
    m_hA.limits["slope"]     = (-50, 0)   # must be negative for decay
    m_hA.limits["exp"]       = (0.1, 5)   # positive, reasonable range
    m_hA.limits["floor"]     = (0, 2)
    m_hA.migrad()   # find minimum
    m_hA.hesse()    # compute uncertainties

    slope      = m_hA.values["slope"]
    slopeErr   = m_hA.errors["slope"]
    yIntercept    = m_hA.values["intercept"]
    yInterceptErr = m_hA.errors["intercept"]

    cost_alt = LeastSquares(ke, pionCarbonhAlt, pionCarbonhAltErr, expo)
    m_alt = Minuit(cost_alt, slope=-1.0, intercept=5, floor=0,exp=1)
    m_alt.limits["intercept"] = (0, 10)
    m_alt.limits["slope"]     = (-50, 0)   # must be negative for decay
    m_alt.limits["exp"]       = (0.1, 5)   # positive, reasonable range
    m_alt.limits["floor"]     = (0, 2)
    m_alt.migrad()
    m_alt.hesse()

    slopeAlt      = m_alt.values["slope"]
    slopeAltErr   = m_alt.errors["slope"]
    yInterceptAlt    = m_alt.values["intercept"]
    yInterceptAltErr = m_alt.errors["intercept"]
    floorAlt    = m_alt.values["floor"]
    floor = m_hA.values["floor"]   
    expAlt    = m_alt.values["exp"]
    exp = m_hA.values["exp"]    
    print(expAlt, exp, slopeAlt, slope)
    fittedVal=[expo(k,slope,yIntercept,floor,exp) for k in ke]
    fittedValAlt=[expo(k,slopeAlt,yInterceptAlt,floorAlt,expAlt) for k in ke]
    plt.plot(ke,fittedVal,"-", label="Fitted Slope for hA2018")
    plt.plot(ke,fittedValAlt,":", label="Fitted Slope for "+str(indexStr))
    plt.ylim(0,3)
    #plt.plot(ke,simVal,"--", label="hA2018 Underlying Model")

    plt.legend(loc="upper right")
    plt.ylabel(ylabel)
    plt.ylim(0,5)
    plt.savefig("../plotting/linearPlotKEValueGamma_"+str(probe)+"_"+str(tgt)+"_"+str(typ)+"_"+str(threshStr)+".pdf")

    plt.show()
    plt.close()
    return {
        "category":         "Sum",
        "quantity":         "Gamma",
        "hA2018_slope":     slope,
        "hA2018_intercept": yIntercept,
        "hA2018_floor":     floor,
        "hA2018_exp":       exp,
        "alt_slope":        slopeAlt,
        "alt_intercept":    yInterceptAlt,
        "alt_floor":        floorAlt,
        "alt_exp":          expAlt
    }
def linearPlotEnergy(probe,tgt,index,typ,thresh=0):
    threshStr=""
    if (thresh==1):
        threshStr="Thresh"
    if (thresh==2):
        threshStr="CompoundThresh"
    pionCarbonhA      = []
    pionCarbonhAlt    = []
    pionCarbonhAErr   = []
    pionCarbonhAltErr = []
    ke=[]
    A=int(str(tgt)[6:9])
    with open('nucleonKO'+str(typ)+str(threshStr)+'.txt', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            if(row[0]=="probe"): continue
            if(row[0]!=str(probe)): continue
            if(row[2]!=str(tgt)): continue
            if (float(row[1])<minKE or float(row[1])>maxKE): continue;

            pionCarbonhA.append(float(row[3]))
            pionCarbonhAErr.append(float(row[4]))
            pionCarbonhAlt.append(float(row[3 + 4 * index]))
            pionCarbonhAltErr.append(float(row[4 + 4 * index]))

            ke.append(float(row[1]))
    indexStr="INCL++"
    if (index==1): indexStr="hN2018"
    if (index==3): indexStr="Geant4"
    plt.plot(ke,pionCarbonhA,"o",label="Mean Value Measured for hA2018")
    plt.plot(ke,pionCarbonhAlt,"o",label="Mean Value Measured for "+str(indexStr))

    plt.xlabel("Kinetic Energy [GeV]")
    probeLabel=r"$p$"

    plt.title("Proton Knockout Final State for "+str(probeLabel)+" for A="+str(int(str(int(tgt/10-1E8))[-3:]))+" Target")
    ylabel="Mean "+typ+" of Nucleons"
    if (thresh==1): ylabel="Mean "+typ+" of Nucleons above Threshold"
    cost_hA = LeastSquares(ke, pionCarbonhA, pionCarbonhAErr, linear)
    m_hA = Minuit(cost_hA, slope=1.0, intercept=0.0)
    m_hA.migrad()   # find minimum
    m_hA.hesse()    # compute uncertainties

    slope      = m_hA.values["slope"]
    slopeErr   = m_hA.errors["slope"]
    yIntercept    = m_hA.values["intercept"]
    yInterceptErr = m_hA.errors["intercept"]

    cost_alt = LeastSquares(ke, pionCarbonhAlt, pionCarbonhAltErr, linear)
    m_alt = Minuit(cost_alt, slope=1.0, intercept=0.0)
    m_alt.migrad()
    m_alt.hesse()

    slopeAlt      = m_alt.values["slope"]
    slopeAltErr   = m_alt.errors["slope"]
    yInterceptAlt    = m_alt.values["intercept"]
    yInterceptAltErr = m_alt.errors["intercept"]

    
    
    
    fittedVal=[yIntercept+slope*k for k in ke]

    fittedValAlt=[yInterceptAlt+slopeAlt*k for k in ke]
    simVal=[hAIntranuke2018Calculator(probe,A,k,typ) for k in ke]
    plt.plot(ke,fittedVal,"-", label="Fitted Slope for hA2018")
    plt.plot(ke,fittedValAlt,":", label="Fitted Slope for "+str(indexStr))

    #plt.plot(ke,simVal,"--", label="hA2018 Underlying Model")
    
    plt.legend(loc="lower right")
    plt.ylim(-7,7)
    if (typ=="Difference"): plt.ylim(-5,5)
    plt.ylabel(ylabel)
    plt.savefig("../plotting/linearPlotKEValue_"+str(indexStr)+"_"+str(probe)+"_"+str(tgt)+"_"+str(typ)+"_"+str(threshStr)+".pdf")

    plt.show()
    plt.close()
    return {
        "category": "Difference",
        "quantity": "Mean",
        "hA2018_slope": slope,
        "hA2018_intercept": yIntercept,
        "alt_slope": slopeAlt,
        "alt_intercept": yInterceptAlt
    }
def linearPlotTarget(probe,ke,index,typ,thresh=0):
    threshStr=""
    if (thresh==1):
        threshStr="Thresh"
    if (thresh==2):
        threshStr="CompoundThresh"
    pionCarbonhA=[]
    pionCarbonhAlt=[]
    tgt=[]
    with open('nucleonKO'+str(typ)+str(threshStr)+'.txt', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            if(row[0]=="probe"): continue
            if(row[0]!=str(probe)): continue
            if(row[1]!=str(ke)): continue
            pionCarbonhA.append(float(row[3+4*index]))
            pionCarbonhAlt.append(float(row[3]))
            tgt.append(float(row[2][6:9]))
    probeLabel=r"$p$"

    plt.title("Proton Knockout Final State for "+str(probeLabel)+" for T="+str(1000*ke)+" MeV")
    indexStr="INCL++"
    if (index==1): indexStr="hN2018"
    if (index==3): indexStr="Geant4"
    plt.plot(tgt,pionCarbonhA,"o",label="Mean Value Measured for hA2018")
    plt.plot(tgt,pionCarbonhAlt,"o",label="Mean Value Measured for "+str(indexStr))
    (slope, yIntercept) = np.polyfit(tgt, pionCarbonhA, 1)
    fittedVal=[yIntercept+slope*A for A in tgt]
    (slopeAlt, yInterceptAlt) = np.polyfit(tgt, pionCarbonhAlt, 1)
    fittedValAlt=[yInterceptAlt+slopeAlt*A for A in tgt]
    simVal=[hAIntranuke2018Calculator(probe,A,ke,typ) for A in tgt]
    plt.plot(tgt,fittedVal,"-", label="Fitted Line for hA2018")
    plt.plot(tgt,fittedValAlt,":", label="Fitted Line for "+str(indexStr))
    #plt.plot(tgt,simVal,"--",label="hA2018 Underlying Model")
    plt.xlabel("A-value of Target")
    ylabel="Mean "+typ+" of Nucleons"
    if (thresh==1): ylabel="Mean "+typ+" of Nucleons above Threshold"
    plt.ylabel(ylabel)
    plt.ylim(-7,7)
    if (typ=="Difference"): plt.ylim(-5,5)
    plt.legend(loc="lower right")
    plt.savefig("../plotting/linearPlotAValue_"+str(probe)+"_"+str(ke)+"GeV_"+str(typ)+"_"+str(threshStr)+".pdf")
    plt.show()
    plt.close()
def linearPlotSigmaTarget(probe,ke,index,thresh=0):
    typ="Difference"
    threshStr=""
    if (thresh==1):
        threshStr="Thresh"
    if (thresh==2):
        threshStr="CompoundThresh"
    pionCarbonhA=[]
    pionCarbonhAlt=[]
    tgt=[]
    with open('nucleonKO'+str(typ)+str(threshStr)+'.txt', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            if(row[0]=="probe"): continue
            if(row[0]!=str(probe)): continue
            if(row[1]!=str(ke)): continue
            pionCarbonhA.append(float(row[5+4*index]))
            pionCarbonhAlt.append(float(row[5]))
            tgt.append(float(row[2][6:9]))
    probeLabel=r"$p$"

    plt.title("Proton Knockout Final State for "+str(probeLabel)+" for T="+str(1000*ke)+" MeV")
    indexStr="INCL++"
    if (index==1): indexStr="hN2018"
    if (index==3): indexStr="Geant4"
    plt.plot(tgt,pionCarbonhA,"o",label="Mean Value Measured for hA2018")
    plt.plot(tgt,pionCarbonhAlt,"o",label="Mean Value Measured for "+str(indexStr))
    (slope, yIntercept) = np.polyfit(tgt, pionCarbonhA, 1)
    fittedVal=[yIntercept+slope*A for A in tgt]
    (slopeAlt, yInterceptAlt) = np.polyfit(tgt, pionCarbonhAlt, 1)
    fittedValAlt=[yInterceptAlt+slopeAlt*A for A in tgt]
    simVal=[hAIntranuke2018Calculator(probe,A,ke,typ) for A in tgt]
    plt.plot(tgt,fittedVal,"-", label="Fitted Line for hA2018")
    plt.plot(tgt,fittedValAlt,":", label="Fitted Line for "+str(indexStr))
    #plt.plot(tgt,simVal,"--",label="hA2018 Underlying Model")
    plt.xlabel("A-value of Target")
    ylabel="Std. of "+typ+" of Nucleons"
    if (thresh==1): ylabel="Std. of "+typ+" of Nucleons above Threshold"
    plt.ylabel(ylabel)
    plt.ylim(0,5)
    if (typ=="Difference"): plt.ylim(-5,5)

    plt.savefig("../plotting/linearPlotAValueStd_"+str(indexStr)+"_"+str(probe)+"_"+str(ke)+"GeV_"+str(typ)+"_"+str(threshStr)+".pdf")
    plt.show()
    plt.close()
def linearPlotGammaTarget(probe,ke,index,tgt,thresh=0):
    typ="Sum"
    threshStr=""
    if (thresh==1):
        threshStr="Thresh"
    if (thresh==2):
        threshStr="CompoundThresh"
    pionCarbonhA=[]
    pionCarbonhAlt=[]
    tgt=[]
    with open('nucleonKO'+str(typ)+str(threshStr)+'.txt', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            if(row[0]=="probe"): continue
            if(row[0]!=str(probe)): continue
            if(row[1]!=str(ke)): continue
            pionCarbonhA.append(float(row[5+4*index]))
            pionCarbonhAlt.append(float(row[5]))
            tgt.append(float(row[2][6:9]))
    probeLabel=r"$p$"

    plt.title("Proton Knockout Final State for "+str(probeLabel)+" for T="+str(1000*ke)+" MeV")
    indexStr="INCL++"
    if (index==1): indexStr="hN2018"
    if (index==3): indexStr="Geant4"
    plt.plot(tgt,pionCarbonhA,"o",label="Mean Value Measured for hA2018")
    plt.plot(tgt,pionCarbonhAlt,"o",label="Mean Value Measured for "+str(indexStr))
    
    (slope, yIntercept) = np.polyfit(tgt, pionCarbonhA, 1)
    fittedVal=[yIntercept+slope*A for A in tgt]
    (slopeAlt, yInterceptAlt) = np.polyfit(tgt, pionCarbonhAlt, 1)
    fittedValAlt=[yInterceptAlt+slopeAlt*A for A in tgt]
    simVal=[hAIntranuke2018Std(probe,A,ke,typ) for A in tgt]
    plt.plot(tgt,fittedVal,"-", label="Fitted Line for hA2018")
    plt.plot(tgt,fittedValAlt,":", label="Fitted Line for "+str(indexStr))
    #plt.plot(tgt,simVal,"--",label="hA2018 Underlying Model")
    
    plt.xlabel("A-value of Target")
    ylabel="Decay Constant of Nucleon Mult."
    if (thresh==1): ylabel="Decay Constant of Nucleon Mult. above Threshold"
    plt.ylabel(ylabel)
    plt.ylim(0,5)
    plt.legend(loc="upper right")
    plt.savefig("../plotting/linearPlotAValueGamma_"+str(indexStr)+"_"+str(probe)+"_"+str(ke)+"GeV_"+str(typ)+"_"+str(threshStr)+".pdf")
    plt.show()
    plt.close()
fit_rows = {}

probes=[2212,2112]
targets=[1000060120, 1000180400, 1000080160, 1000260560]
kes = [round(0.225 + 0.1 * i, 4) for i in range(10)]

types=["Sum","Difference"]

indexes=[1,2,3]
indexStrings=["hN2018","INCL","Geant4"]

for probe in probes:
    for index in indexes:
        indexStr=indexStrings[index-1]
        for tgt in targets:
            res_list = [
                linearPlotGammaEnergy(probe, tgt, index, "Sum", thresh=1),
                linearPlotSigmaEnergy(probe, tgt, index, "Difference", thresh=1),
                linearPlotEnergy(probe, tgt, index, "Difference", thresh=1)
            ]
            for res in res_list:
                key = (tgt, probe, res["category"], res["quantity"])
                if key not in fit_rows:
                    fit_rows[key] = {
                        "target":           tgt,
                        "particle":       probe,
                        "category":         res["category"],
                        "quantity":         res["quantity"],
                        "hA2018_slope":     res["hA2018_slope"],
                        "hA2018_intercept": res["hA2018_intercept"],
                        "hA2018_floor":     res.get("hA2018_floor", ""),
                        "hA2018_exp":       res.get("hA2018_exp",   ""),
                    }

                fit_rows[key][f"{indexStr}_slope"]     = res["alt_slope"]
                fit_rows[key][f"{indexStr}_intercept"] = res["alt_intercept"]
                if "alt_floor" in res:
                    fit_rows[key][f"{indexStr}_floor"] = res["alt_floor"]
                if "alt_exp" in res:
                    fit_rows[key][f"{indexStr}_exp"]   = res["alt_exp"]

with open("nucleonKO_keLinearFitParameters.txt", "w", newline="") as f:
    writer = csv.writer(f)
    header = [
        "category", "quantity", "particle","target",
        "hA2018_slope", "hA2018_intercept", "hA2018_floor", "hA2018_exp",
        "hN2018_slope", "hN2018_intercept", "hN2018_floor", "hN2018_exp",
        "INCL_slope",   "INCL_intercept",   "INCL_floor",   "INCL_exp",
        "Geant4_slope", "Geant4_intercept", "Geant4_floor", "Geant4_exp"
    ]
    writer.writerow(header)
    for row in fit_rows.values():
        writer.writerow([
            row["category"], row["quantity"], row["particle"], row["target"],
            row.get("hA2018_slope", ""), row.get("hA2018_intercept", ""), row.get("hA2018_floor", ""), row.get("hA2018_exp", ""),
            row.get("hN2018_slope", ""), row.get("hN2018_intercept", ""), row.get("hN2018_floor", ""), row.get("hN2018_exp", ""),
            row.get("INCL_slope",   ""), row.get("INCL_intercept",   ""), row.get("INCL_floor",   ""), row.get("INCL_exp",   ""),
            row.get("Geant4_slope", ""), row.get("Geant4_intercept", ""), row.get("Geant4_floor", ""), row.get("Geant4_exp", "")
        ])