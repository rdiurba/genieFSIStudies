# -*- coding: utf-8 -*-
"""
Script to plot multiplicity of nucleons as a sum and difference as a function of KE and target
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
minKE=0.325;
maxKE=1.25;
pionCarbonhA=[]
ke=[]
probe=211
tgt=1000060120
index=0
typ="Difference" # "Sum"

def linear(x, slope, intercept):
    return intercept + slope * x
def hAIntranuke2018Calculator(probe,A,ke,typ):
    val=-99

    if (typ=="Sum"):
        val=.0001*(1.+(ke*1000)/250.) * (A-10)*(A-10) + 3.5
    if (typ=="Difference"):
        val=(1.+ke*1000/250.) - ((A/200.)*(1. + 2.*ke*1000/250.))
    return val
def hAIntranuke2018Std(probe,A,ke,typ):
    val=-99
    if (typ=="Sum"):
        val=(10. + 4. * (ke*1000.)/250.)*np.power(A/250.,0.9); 
    if (typ=="Difference"):
        val=4*(1 - np.exp(-0.03*ke*1000.))
    return val
def linearPlotSigmaEnergy(probe,tgt,index,typ,thresh=0):
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
    with open('pionAbs'+str(typ)+str(threshStr)+'.txt', newline='') as csvfile:
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
    probeLabel=r"$\pi^{+}$"
    if (probe==111): probeLabel=r"$\pi^{0}$"
    if (probe==-211): probeLabel=r"$\pi^{-}$"
    plt.title("Pion Absorption Final State for "+str(probeLabel)+" for A="+str(int(str(int(tgt/10-1E8))[-3:]))+" Target")
    ylabel="Std. of "+typ+" of Nucleons"
    if (threshStr==1): ylabel="Std. of "+typ+" of Nucleons above Threshold"
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
    plt.legend()
    plt.ylim(0,5)
    if (max(pionCarbonhAlt)>5): plt.ylim(0,10)

    plt.ylabel(ylabel)
    plt.savefig("../plotting/linearPlotKEValueStd_"+str(indexStr)+"_"+str(probe)+"_"+str(tgt)+"_"+str(typ)+"_"+str(threshStr)+".pdf")

    plt.show()
    plt.close()
    return {
        "category": typ,
        "quantity": "Std",
        "hA2018_slope": slope,
        "hA2018_intercept": yIntercept,
        "alt_slope": slopeAlt,
        "alt_intercept": yInterceptAlt
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
    with open('pionAbs'+str(typ)+str(threshStr)+'.txt', newline='') as csvfile:
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
    probeLabel=r"$\pi^{+}$"
    if (probe==111): probeLabel=r"$\pi^{0}$"
    if (probe==-211): probeLabel=r"$\pi^{-}$"
    plt.title("Pion Absorption Final State for "+str(probeLabel)+" for A="+str(int(str(int(tgt/10-1E8))[-3:]))+" Target")
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
    plt.legend()
    plt.ylim(0,10)
    if (typ=="Difference"): plt.ylim(-5,5)

    plt.ylabel(ylabel)
    plt.savefig("../plotting/linearPlotKEValue_"+str(indexStr)+"_"+str(probe)+"_"+str(tgt)+"_"+str(typ)+"_"+str(threshStr)+".pdf")

    plt.show()
    plt.close()
    return {
        "category": typ,
        "quantity": "Mean",
        "hA2018_slope": slope,
        "hA2018_intercept": yIntercept,
        "alt_slope": slopeAlt,
        "alt_intercept": yInterceptAlt
    }
def linearPlotTarget(probe,ke,index,tgt,typ,thresh=0):
    threshStr=""
    if (thresh==1):
        threshStr="Thresh"
    if (thresh==2):
        threshStr="CompoundThresh"
    pionCarbonhA      = []
    pionCarbonhAlt    = []
    pionCarbonhAErr   = []
    pionCarbonhAltErr = []
    tgt=[]
    with open('pionAbs'+str(typ)+str(threshStr)+'.txt', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            if(row[0]=="probe"): continue
            if(row[0]!=str(probe)): continue
            if(row[1]!=str(ke)): continue

            pionCarbonhA.append(float(row[3]))
            pionCarbonhAErr.append(float(row[4]))
            pionCarbonhAlt.append(float(row[3 + 4 * index]))
            pionCarbonhAltErr.append(float(row[4 + 4 * index]))
            tgt.append(float(row[2][6:9]))
    probeLabel=r"$\pi^{+}$"
    if (probe==111): probeLabel=r"$\pi^{0}$"
    if (probe==-211): probeLabel=r"$\pi^{-}$"
    plt.title("Pion Absorption Final State for "+str(probeLabel)+" for T="+str(1000*ke)+" MeV")
    indexStr="INCL++"
    if (index==1): indexStr="hN2018"
    if (index==3): indexStr="Geant4"
    plt.plot(tgt,pionCarbonhA,"o",label="Mean Value Measured for hA2018")
    plt.plot(tgt,pionCarbonhAlt,"o",label="Mean Value Measured for "+str(indexStr))
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
    
    
    fittedVal=[yIntercept+slope*A for A in tgt]

    fittedValAlt=[yInterceptAlt+slopeAlt*A for A in tgt]
    simVal=[hAIntranuke2018Calculator(probe,A,ke,typ) for A in tgt]
    plt.plot(tgt,fittedVal,"-", label="Fitted Line for hA2018")
    plt.plot(tgt,fittedValAlt,":", label="Fitted Line for "+str(indexStr))
    #plt.plot(ke,simVal,"--", label="hA2018 Underlying Model")

    plt.xlabel("A-value of Target")
    ylabel="Mean "+typ+" of Nucleons"
    if (thresh==1): ylabel="Mean "+typ+" of Nucleons above Threshold"
    plt.ylabel(ylabel)
    plt.legend()
    plt.ylim(0,5)
    if (typ=="Difference"): plt.ylim(-5,5)
    plt.savefig("../plotting/linearPlotAValue_"+str(indexStr)+"_"+str(probe)+"_"+str(ke)+"GeV_"+str(typ)+"_"+str(threshStr)+".pdf")
    plt.show()
    plt.close()
def linearPlotSigmaTarget(probe,ke,index,tgt,typ,thresh=0):
    threshStr=""
    if (thresh==1):
        threshStr="Thresh"
    if (thresh==2):
        threshStr="CompoundThresh"
    pionCarbonhA=[]
    pionCarbonhAlt=[]
    tgt=[]
    with open('pionAbs'+str(typ)+str(threshStr)+'.txt', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            if(row[0]=="probe"): continue
            if(row[0]!=str(probe)): continue
            if(row[1]!=str(ke)): continue
            pionCarbonhA.append(float(row[5+4*index]))
            pionCarbonhAlt.append(float(row[5]))
            tgt.append(float(row[2][6:9]))
    probeLabel=r"$\pi^{+}$"
    if (probe==111): probeLabel=r"$\pi^{0}$"
    if (probe==-211): probeLabel=r"$\pi^{-}$"
    plt.title("Pion Absorption Final State for "+str(probeLabel)+" for T="+str(1000*ke)+" MeV")
    indexStr="INCL++"
    if (index==1): indexStr="hN2018"
    if (index==3): indexStr="Geant4"
    plt.plot(tgt,pionCarbonhA,"o",label="Mean Value Measured for hA2018")
    plt.plot(tgt,pionCarbonhAlt,"o",label="Mean Value Measured for "+str(indexStr))
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
    
    
    
    fittedVal=[yIntercept+slope*A for A in tgt]

    fittedValAlt=[yInterceptAlt+slopeAlt*A for A in tgt]
    simVal=[hAIntranuke2018Calculator(probe,A,ke,typ) for A in tgt]

    plt.plot(tgt,fittedVal,"-", label="Fitted Line for hA2018")
    plt.plot(tgt,fittedValAlt,":", label="Fitted Line for "+str(indexStr))
    plt.plot(tgt,simVal,"--",label="hA2018 Underlying Model")
    plt.xlabel("A-value of Target")
    ylabel="Std. of "+typ+" of Nucleons"
    if (thresh==1): ylabel="Std. of "+typ+" of Nucleons above Threshold"
    plt.ylabel(ylabel)
    plt.ylim(0,5)
    plt.legend()
    plt.savefig("../plotting/linearPlotAValueStd_"+str(indexStr)+"_"+str(probe)+"_"+str(ke)+"GeV_"+str(typ)+"_"+str(threshStr)+".pdf")
    plt.show()
    plt.close()
probes=[211,111,-211]
targets=[1000060120, 1000180400, 1000080160, 1000260560]

fit_rows = {}
types=["Sum","Difference"]

indexes=[1,2,3]
indexStrings=["hN2018","INCL","Geant4"]

for probe in probes:
    for index in indexes:
        indexStr=indexStrings[index-1]
        for typ in types:
            for tgt in targets:
                res_list = [linearPlotEnergy(probe,tgt,index,typ,thresh=1),
                linearPlotSigmaEnergy(probe,tgt,index,typ,thresh=1)]
                for res in res_list:
                    key = (tgt,probe, res["category"], res["quantity"])
                
                    if key not in fit_rows:
                        fit_rows[key] = {
                            "target": tgt,
                            "particle": probe,
                            "category": res["category"],
                            "quantity": res["quantity"],
                            "hA2018_slope": res["hA2018_slope"],
                            "hA2018_intercept": res["hA2018_intercept"]
                        }
    
                    fit_rows[key][f"{indexStr}_slope"] = res["alt_slope"]
                    fit_rows[key][f"{indexStr}_intercept"] = res["alt_intercept"]


with open("pionAbs_keLinearFitParameters.txt", "w", newline="") as f:
    writer = csv.writer(f)

    header = [
        "category",
        "quantity",
        "particle",
        "target",
        "hA2018_slope", "hA2018_intercept",
        "hN2018_slope", "hN2018_intercept",
        "INCL_slope", "INCL_intercept",
        "Geant4_slope", "Geant4_intercept"
    ]
    writer.writerow(header)

    for row in fit_rows.values():
        writer.writerow([
            row["category"],
            row["quantity"],
            row["particle"],
            row["target"],
            row.get("hA2018_slope", ""),
            row.get("hA2018_intercept", ""),
            row.get("hN2018_slope", ""),
            row.get("hN2018_intercept", ""),
            row.get("INCL_slope", ""),
            row.get("INCL_intercept", ""),
            row.get("Geant4_slope", ""),
            row.get("Geant4_intercept", "")
        ])

