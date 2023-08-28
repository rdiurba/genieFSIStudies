#!/bin/bash

export events="2000000"
export version="v3_4_0"
export tune="AR23_20i_00_000"
export probe="-14"
export target="1000180400"
export interaction="Default"
export minE="0."
export maxE="50."
export fluxfile="/dune/app/users/kunzmann/BuildEventGenerators/BuildEventGenerators/jobcards/g4numiv6_minervame_me000z-200i_0_0001.dk2nu"
export fluxhisto="hEnumu_cv"
export GEOM="/dune/app/users/kunzmann/BuildEventGenerators/BuildEventGenerators/jobcards/Merged2x2MINERVA_noRock.gdml"
export DET="ProtoDUNE-ND"
export GXMLPATH="/dune/app/users/kunzmann/BuildEventGenerators/BuildEventGenerators/jobcards/GNuMIFlux.xml:$GXMLPATH"
export GDK2NUFLUXXML="/dune/app/users/kunzmann/BuildEventGenerators/BuildEventGenerators/jobcards/GNuMIFlux.xml"
export GENIEXSECFILE2="/pnfs/dune/persistent/users/rdiurba/gxspl-NUsmall.xml"
# Produce the GENIE splines
#gmkspl -p 14,-14,12,-12 -t ${target} -e ${maxE} -o ${target}_${interaction}_${version}_${tune}_DUNE.xml --tune ${tune} --event-generator-list ${interaction}

# Convert the xml splines to root format
#gspl2root -f /dune/app/users/kunzmann/BuildEventGenerators/BuildEventGenerators/jobcards/splines/${probe}_${target}_${interaction}_${version}_${tune}_DUNE.xml --event-generator-list ${interaction} -p ${probe} -t ${target} -o splines/${probe}_${target}_${interaction}_${version}_${tune}_DUNE.root --tune ${tune}

# Generate GENIE events
gevgen_fnal -o /dune/data/users/rdiurba/test -n $events -g ${target} --event-generator-list ${interaction} --tune ${tune} --cross-sections ${GENIEXSECFILE2} -f ${fluxfile},${DET} --seed 114101 
#old with maxpath file
#gevgen_fnal -n $events -g ${target} --event-generator-list ${interaction} --tune ${tune} --cross-sections splines/${probe}_${target}_${interaction}_${version}_${tune}_DUNE.xml -f ${fluxfile},${DET},14 -m ${MAXPATH_FILE} -o samples/${probe}_${target}_${interaction}_${version}_${tune}_DUNE.ghep.root

# Convert file from ghep to gst
gntpc -f gst -i /dune/data/users/rdiurba/test.0.ghep.root -o /dune/data/users/rdiurba/test.gst.root
# Convert file from ghep to nuisance format
#PrepareGENIE -i samples/${probe}_${target}_${interaction}_${version}_${tune}_DUNE.ghep.root -t ${target}[1] -o samples/${probe}_${target}_${interaction}_${version}_${tune}_DUNE.gprep.root -f ${fluxfile},${fluxhisto}

# Convert to nuisance flat tree format
#nuisflat -i GENIE:samples/${probe}_${target}_${interaction}_${version}_${tune}_DUNE.gprep.root -o samples/${probe}_${target}_${interaction}_${version}_${tune}_DUNE.flat.root

# Remove all unnecessary files
#rm *.status 
#rm input-flux.root
