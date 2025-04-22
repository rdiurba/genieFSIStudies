# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 15:13:03 2025

@author: Richie-LHEP
"""

import csv
import matplotlib.pyplot as plt
import numpy as np
pionCarbonhA=[]
ke=[]
probe=211
tgt=1000060120
index=0
typ="Difference" # "Sum"
def hAIntranuke2018Calculator(probe,A,ke,typ):
    val=-99
    if (typ=="Sum"):
        val=-9999
    if (typ=="Difference"):
          Z=6
          if (A==40): Z=18 
          if (A==56): Z=26
          if (A==208): Z=82 

          if (A-Z > Z):
            nd0 =  135.227 * np.exp(-7.124*(A-Z)/(A)) - 2.762;
          else:
            nd0 = -135.227 * np.exp(-7.124*Z /A) + 4.914;
          val=nd0
    return val
def hAIntranuke2018Std(probe,A,ke,typ):
    val=-99
    ke=ke*1000.0
    if (typ=="Sum"):
          c1 = 0.041 + ke * 0.0001525;
          c2 = -0.003444 - ke * 0.00002324;

          c3 = 0.064 - ke * 0.000015;
          gam_ns = c1 * np.exp(c2*A) + c3;
          if(gam_ns<0.002): gam_ns = 0.002;
          val=gam_ns
    if (typ=="Difference"):
        val=2.034 + A * 0.007846
    return val
def linearPlotSigmaEnergy(probe,tgt,index,typ,thresh=0):
    threshStr=""
    if (thresh==1):
        threshStr="Thresh"
    pionCarbonhA=[]
    pionCarbonhAlt=[]
    ke=[]
    A=int(str(tgt)[6:9])
    with open('nucleonKO'+str(typ)+str(threshStr)+'.txt', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            if(row[0]=="probe"): continue
            if(row[0]!=str(probe)): continue
            if(row[2]!=str(tgt)): continue
            pionCarbonhA.append(float(row[5]))
            pionCarbonhAlt.append(float(row[5+4*index]))

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
    (slope, yIntercept) = np.polyfit(ke, pionCarbonhA, 1)
    fittedVal=[yIntercept+slope*k for k in ke]
    (slopeAlt, yInterceptAlt) = np.polyfit(ke, pionCarbonhAlt, 1)
    fittedValAlt=[yInterceptAlt+slopeAlt*k for k in ke]
    simVal=[hAIntranuke2018Calculator(probe,A,k,typ) for k in ke]
    plt.plot(ke,fittedVal,"-", label="Fitted Slope for hA2018")
    plt.plot(ke,fittedValAlt,":", label="Fitted Slope for "+str(indexStr))
    plt.ylim(-7,7)
    plt.plot(ke,simVal,"--", label="hA2018 Underlying Model")
    plt.legend(loc="lower right")
    plt.ylabel(ylabel)
    plt.savefig("linearPlotKEValue_"+str(probe)+"_"+str(tgt)+"_"+str(typ)+"_"+str(threshStr)+".pdf")

    plt.show()
def linearPlotGammaEnergy(probe,tgt,index,typ,thresh=0):
    threshStr=""
    if (thresh==1):
        threshStr="Thresh"
    pionCarbonhA=[]
    pionCarbonhAlt=[]
    ke=[]
    A=int(str(tgt)[6:9])
    with open('nucleonKO'+str(typ)+str(threshStr)+'.txt', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            if(row[0]=="probe"): continue
            if(row[0]!=str(probe)): continue
            if(row[2]!=str(tgt)): continue
            pionCarbonhA.append(float(row[5]))
            pionCarbonhAlt.append(float(row[5+4*index]))

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
    (slope, yIntercept) = np.polyfit(ke, pionCarbonhA, 1)
    fittedVal=[yIntercept+slope*k for k in ke]
    (slopeAlt, yInterceptAlt) = np.polyfit(ke, pionCarbonhAlt, 1)
    fittedValAlt=[yInterceptAlt+slopeAlt*k for k in ke]
    simVal=[hAIntranuke2018Calculator(probe,A,k,typ) for k in ke]
    plt.plot(ke,fittedVal,"-", label="Fitted Slope for hA2018")
    plt.plot(ke,fittedValAlt,":", label="Fitted Slope for "+str(indexStr))
    plt.ylim(-7,7)
    plt.plot(ke,simVal,"--", label="hA2018 Underlying Model")
    plt.legend(loc="lower right")
    plt.ylabel(ylabel)
    plt.savefig("linearPlotKEValueGamma_"+str(probe)+"_"+str(tgt)+"_"+str(typ)+"_"+str(threshStr)+".pdf")

    plt.show()
def linearPlotEnergy(probe,tgt,index,typ,thresh=0):
    threshStr=""
    if (thresh==1):
        threshStr="Thresh"
    pionCarbonhA=[]
    pionCarbonhAlt=[]
    ke=[]
    A=int(str(tgt)[6:9])
    with open('nucleonKO'+str(typ)+str(threshStr)+'.txt', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            if(row[0]=="probe"): continue
            if(row[0]!=str(probe)): continue
            if(row[2]!=str(tgt)): continue
            pionCarbonhA.append(float(row[3]))
            pionCarbonhAlt.append(float(row[3+4*index]))

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
    (slope, yIntercept) = np.polyfit(ke, pionCarbonhA, 1)
    fittedVal=[yIntercept+slope*k for k in ke]
    (slopeAlt, yInterceptAlt) = np.polyfit(ke, pionCarbonhAlt, 1)
    fittedValAlt=[yInterceptAlt+slopeAlt*k for k in ke]
    simVal=[hAIntranuke2018Calculator(probe,A,k,typ) for k in ke]
    plt.plot(ke,fittedVal,"-", label="Fitted Slope for hA2018")
    plt.plot(ke,fittedValAlt,":", label="Fitted Slope for "+str(indexStr))

    plt.plot(ke,simVal,"--", label="hA2018 Underlying Model")
    plt.legend(loc="lower right")
    plt.ylim(-7,7)
    plt.ylabel(ylabel)
    plt.savefig("linearPlotKEValue_"+str(probe)+"_"+str(tgt)+"_"+str(typ)+"_"+str(threshStr)+".pdf")

    plt.show()
def linearPlotTarget(probe,ke,index,tgt,typ,thresh=0):
    threshStr=""
    if (thresh==1):
        threshStr="Thresh"
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
    plt.plot(tgt,simVal,"--",label="hA2018 Underlying Model")
    plt.xlabel("A-value of Target")
    ylabel="Mean "+typ+" of Nucleons"
    if (thresh==1): ylabel="Mean "+typ+" of Nucleons above Threshold"
    plt.ylabel(ylabel)
    plt.ylim(-7,7)
    plt.legend(loc="lower right")
    plt.savefig("linearPlotAValue_"+str(probe)+"_"+str(ke)+"GeV_"+str(typ)+"_"+str(threshStr)+".pdf")
    plt.show()
def linearPlotSigmaTarget(probe,ke,index,tgt,typ,thresh=0):
    threshStr=""
    if (thresh==1):
        threshStr="Thresh"
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
    plt.plot(tgt,simVal,"--",label="hA2018 Underlying Model")
    plt.xlabel("A-value of Target")
    ylabel="Std. of "+typ+" of Nucleons"
    if (thresh==1): ylabel="Std. of "+typ+" of Nucleons above Threshold"
    plt.ylabel(ylabel)
    plt.ylim(-7,7)
    plt.legend(loc="lower right")
    plt.savefig("linearPlotAValue_"+str(probe)+"_"+str(ke)+"GeV_"+str(typ)+"_"+str(threshStr)+".pdf")
    plt.show()
    
def linearPlotGammaTarget(probe,ke,index,tgt,typ,thresh=0):
    threshStr=""
    if (thresh==1):
        threshStr="Thresh"
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
    plt.plot(tgt,simVal,"--",label="hA2018 Underlying Model")
    plt.xlabel("A-value of Target")
    ylabel="Decay Constant of Nucleon Mult."
    if (thresh==1): ylabel="Decay Constant of Nucleon Mult. above Threshold"
    plt.ylabel(ylabel)
    plt.ylim(-7,7)
    plt.legend(loc="lower right")
    plt.savefig("linearPlotAValueGamma_"+str(probe)+"_"+str(ke)+"GeV_"+str(typ)+"_"+str(threshStr)+".pdf")
    plt.show()
probes=[2212]
indexes=[0,1,2,3]
targets=[1000060120,1000180400,1000260560,1000822080]
kes=[0.125,0.25,0.5]
types=["Sum","Difference"]
#kes=[0.125]
#targets=[1000060120]
indexes=[1,2,3]
for probe in probes:
    for index in indexes:
        for tgt in targets:
            linearPlotEnergy(probe, tgt, index, "Difference",thresh=1)
            linearPlotEnergy(probe, tgt, index, "Difference",thresh=0)
            linearPlotSigmaEnergy(probe, tgt, index, "Difference",thresh=1)
            linearPlotSigmaEnergy(probe, tgt, index, "Difference",thresh=0)
            linearPlotGammaEnergy(probe, tgt, index, "Sum",thresh=1)
            linearPlotGammaEnergy(probe, tgt, index, "Sum",thresh=0)
        for ke in kes:
            linearPlotTarget(probe, ke, index, tgt, "Difference",thresh=1)
            linearPlotTarget(probe, ke, index, tgt, "Difference",thresh=0)
            linearPlotSigmaTarget(probe, ke, index, "Difference", typ,thresh=1)
            linearPlotSigmaTarget(probe, ke, index, "Difference", typ,thresh=0)
            linearPlotGammaTarget(probe, ke, index, "Sum", typ,thresh=1)
            linearPlotGammaTarget(probe, ke, index, "Sum", typ,thresh=0)


