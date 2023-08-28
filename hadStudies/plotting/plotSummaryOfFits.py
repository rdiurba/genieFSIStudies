import csv
import matplotlib.pyplot as plt
import numpy as np
pionCarbonhA=[]
ke=[]
probe=211
tgt=1000060120
index=0
typ="Difference" # "Sum"
def hAIntranuke2018CalculatorStd(probe,A,ke,typ):
    std=0
    if (typ=="Sum"):
        std = (10. + 4. * ke/250.)*np.power(A/250.,0.9);  
    if (typ=="Difference"):
        std = 4*(1 - np.exp(-0.03*ke));
    return std
def hAIntranuke2018Calculator(probe,A,ke,typ):
    val=-99
    if (typ=="Sum"):
        val=.0001*(1.+(ke*1000)/250.) * (A-10)*(A-10) + 3.5
    if (typ=="Difference"):
        val=(1.+ke*1000/250.) - ((A/200.)*(1. + 2.*ke*1000/250.))
    """
    if (probe!=211 and probe!=111):
        if (typ=="Sum"):
            ke=ke*1000
            c1 = 0.041 + ke * 0.0001525;
            c2 = -0.003444 - ke * 0.00002324;
            c3 = 0.064 - ke * 0.000015;
            val = c1 * np.exp(c2*A) + c3;
            print(val)
            if(val<0.002): val = 0.002;
        if (typ=="Difference"):
            val=(1.+ke*1000/250.) - ((A/200.)*(1. + 2.*ke*1000/250.))
      """      

    return val
def linearPlotEnergy(probe,tgt,index,typ):
    pionCarbonhA=[]
    ke=[]
    A=int(str(tgt)[6:9])
    fileTitle="pionAbs"
    if (probe==2112 or probe==2212): fileTitle="nucleonKO"
    with open(fileTitle+str(typ)+'.txt', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            if(row[0]=="probe"): continue
            if(row[0]!=str(probe)): continue
            if(row[2]!=str(tgt)): continue
            if ((typ=="Difference" and (abs(probe)==2212 or probe==2112)) or (abs(probe)==211 or probe==111)):
                pionCarbonhA.append(float(row[3+4*index]))
            else:
                pionCarbonhA.append(-float(row[5+4*index]))
            ke.append(float(row[1]))
    indexStr="INCL++"
    if (index==3): indexStr="Geant4"
    plt.plot(ke,pionCarbonhA,"o",label="Simulated Points for "+str(indexStr))
    plt.xlabel("Kinetic Energy [GeV]")
    probeLabel=r"$\pi^{+}$"
    probeString="piPlus"
    if (probe==111): 
        probeLabel=r"$\pi^{0}$"
        probeString="pi0"
    if (probe==-211):
        probeLabel=r"$\pi^{-}$"
        probeString="piMinus"
    if (probe==2112):
        probeLabel=r"$n$"
        probeString="neutron"
    if (probe==2212):
        probeLabel=r"$p^{+}$"
        probeString="proton"
    scatter="Pion Absorption"
    if (probe==2112 or probe==2212): scatter="Nucleon Knockout"
    plt.title(scatter+" Final State for "+str(probeLabel)+" for A="+str(int(str(int(tgt/10-1E8))[-3:]))+" Target")
    variableType="Mean"
    if ((probe==2212 or probe==2112) and typ=="Sum"):
        variableType="Gamma"
    ylabel=variableType+" "+typ+" of Nucleons from "+scatter
    (slope, yIntercept) = np.polyfit(ke, pionCarbonhA, 1)
    fittedVal=[yIntercept+slope*k for k in ke]
    simVal=[hAIntranuke2018Calculator(probe,A,k,typ) for k in ke]
    plt.plot(ke,fittedVal,"-", label="Fitted Slope for "+str(indexStr))
    plt.plot(ke,simVal,"--", label="hA2018 Model Slope")
    plt.legend()
    plt.ylabel(ylabel)
    plt.savefig(typ+probeString+str(indexStr)+str(tgt)+"MeV.png")
    plt.show()
def linearPlotTarget(probe,ke,index,tgt,typ):
    pionCarbonhA=[]
    tgt=[]
    fileTitle="pionAbs"
    if (probe==2112 or probe==2212): fileTitle="nucleonKO"
    with open(fileTitle+str(typ)+'.txt', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            if(row[0]=="probe"): continue
            if(row[0]!=str(probe)): continue
            if(row[1]!=str(ke)): continue
            if ((typ=="Difference" and (abs(probe)==2212 or probe==2112)) or (abs(probe)==211 or probe==111)):
                pionCarbonhA.append(float(row[3+4*index]))
            else:
                pionCarbonhA.append(-float(row[5+4*index]))
            tgt.append(float(row[2][6:9]))
    probeLabel=r"$\pi^{+}$"
    probeString="piPlus"
    if (probe==111): 
        probeLabel=r"$\pi^{0}$"
        probeString="pi0"
    if (probe==-211):
        probeLabel=r"$\pi^{-}$"
        probeString="piMinus"
    if (probe==2112):
        probeLabel=r"$n$"
        probeString="neutron"
    if (probe==2212):
        probeLabel=r"$p^{+}$"
        probeString="proton"
    scatter="Pion Absorption"
    if (probe==2112 or probe==2212): scatter="Nucleon Knockout"
    variableType="Mean"
    if ((probe==2212 or probe==2112) and typ=="Sum"):
        variableType="Gamma"
    ylabel=variableType+" "+typ+" of Nucleons from "+scatter
    plt.title(scatter+" Final State for "+str(probeLabel)+" for T="+str(1000*ke)+" MeV")
    indexStr="INCL++"
    if (index==3): indexStr="Geant4"
    plt.plot(tgt,pionCarbonhA,"o",label="Simulated Points for "+str(indexStr))
    if(typ=="Sum"):
        (yIntercept,slope,quadratic) = np.polyfit(tgt, pionCarbonhA, 2)
        fittedVal=[yIntercept+slope*A+quadratic*A*A for A in tgt]
    else:
        (yIntercept,slope) = np.polyfit(tgt, pionCarbonhA, 1)
        fittedVal=[yIntercept+slope*A for A in tgt]
    simVal=[hAIntranuke2018Calculator(probe,A,ke,typ) for A in tgt]
    plt.plot(tgt,fittedVal,"-", label="Fitted Slope for "+str(indexStr))
    plt.plot(tgt,simVal,"--",label="hA2018 Model Slope")
    plt.xlabel("A-value of Target")
    ylabel="Mean "+typ+" of Nucleons from "+scatter
    plt.ylabel(ylabel)
    plt.legend()
    plt.savefig(typ+probeString+str(indexStr)+str(1000*ke)+"MeV.png")
    plt.show()

probes=[211,-211,111,2212,2112]
indexes=[0,1,2,3]
targets=[1000060120,1000260560,1000822080]
kes=[0.125,0.25,0.5]
types=["Sum","Difference"]
probes=[211]
#kes=[0.125]
#targets=[1000060120]
indexes=[2,3]
for probe in probes:
    for index in indexes:
        for tgt in targets:
            for typ in types: 
                linearPlotEnergy(probe, tgt, index, typ)
            #for ke in kes:
            #    linearPlotTarget(probe, ke, index, tgt, typ)
print(hAIntranuke2018Calculator(211, 208, 0.5, "Sum"))
