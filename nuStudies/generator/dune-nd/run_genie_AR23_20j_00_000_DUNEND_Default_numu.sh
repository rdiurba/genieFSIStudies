#!/bin/bash
export events="2000000"
export version="v3_4_0"
export tune="AR23_20j_00_000"
export probe="14"
export target="1000180400"
export interaction="Default"
export minE="0."
export maxE="50."
export fluxfile="/pnfs/dune/persistent/users/rdiurba/DUNE_OptimizedEngineeredNov2017_REGULAR.root"
export fluxhisto="numu_NDFHC_flux"
export GENIEXSECFILE2="/pnfs/dune/persistent/users/rdiurba/gxspl-NUsmall_AR23_20j_00_000.xml"
# Produce the GENIE splines
##gmkspl -p ${probe} -t ${target} -e ${maxE} -o ${probe}_${target}_${interaction}_${version}_${tune}.xml --tune ${tune} --event-generator-list ${interaction}

# Convert the xml splines to root format
##gspl2root -f ${probe}_${target}_${interaction}_${version}_${tune}.xml --event-generator-list ${interaction} -p ${probe} -t ${target} -o ${probe}_${target}_${interaction}_${version}_${tune}.xml.root --tune ${tune}

# Generate GENIE events
gevgen -n $events -p ${probe} -t ${target} -e ${minE},${maxE} --message_thresholds $GENIE/config/Messenger_whisper.xml  --event-generator-list ${interaction} --tune ${tune} --cross-sections ${GENIEXSECFILE2} -f ${fluxfile},${fluxhisto} -o /dune/data/users/rdiurba/genieFSIStudies/DUNEND_${probe}_${target}_${interaction}_${tune}.ghep.root --seed 121111
# Convert file from ghep to gst
gntpc -f gst -i /dune/data/users/rdiurba/genieFSIStudies/DUNEND_${probe}_${target}_${interaction}_${tune}.ghep.root -o /dune/data/users/rdiurba/genieFSIStudies/DUNEND_${probe}_${target}_${interaction}_${tune}.gst.root

# Convert file from ghep to nuisance format
#PrepareGENIE -i ${probe}_${target}_${interaction}_${version}_${tune}.ghep.root -t ${target}[1] -o ${probe}_${target}_${interaction}_${version}_${tune}.gprep.root -f ${fluxfile},${fluxhisto}

# Convert to nuisance flat tree format
#nuisflat -i GENIE:${probe}_${target}_${interaction}_${version}_${tune}.gprep.root -o ${probe}_${target}_${interaction}_${version}_${tune}.flat.root

# Remove all unnecessary files
rm *.status 
rm input-flux.root
