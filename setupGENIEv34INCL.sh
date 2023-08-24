#!/bin/bash

source /grid/fermiapp/products/dune/setup_dune.sh
#setup root v6_22_08d -q e20:p392:prof
#setup pythia v6_4_28r -q gcc930:prof
#setup log4cpp v1_1_3c -q e20:prof
#setup lhapdf v6_3_0 -q e20:p392:prof
setup dk2nugenie v01_10_01k -q "e20:prof"
unsetup genie
setup genie v3_04_00d -q "e20:+inclxx:prof"
export GXMLPATH=`pwd`/genieTunes:$GXMLPATH

#setup genie v3_04_00 -q "e20:prof"
# GENIE
export GENIE_VERSION=v3_04_00d


export LD_LIBRARY_PATH=${GENIE_FQ_DIR}/lib/:${LD_LIBRARY_PATH}
export PATH=${GENIE_FQ_DIR}/bin/:${PATH}
export ROOT_INCLUDE_PATH=${GENIE_FQ_DIR}/include/GENIE/:${ROOT_INCLUDE_PATH}
export CMAKE_PREFIX_PATH=${GENIE_FQ_DIR}:${CMAKE_PREFIX_PATH}
export PKG_CONFIG_PATH=${GENIE_FQ_DIR}:${PKG_CONFIG_PATH}

export NuSystDir=`pwd`
export NUSYST_INC=${NuSystDir}/build/Linux/include
export NUSYST_LIB=${NuSystDir}/build/Linux/lib

export NUSYST_DQ_DIR=${NuSystDir}/build/Linux/

export LD_LIBRARY_PATH=${NUSYST_DQ_DIR}/lib/:$LD_LIBRARY_PATH
export PATH=${NUSYST_DQ_DIR}/bin/:$PATH
export ROOT_INCLUDE_PATH=${NUSYST_DQ_DIR}/include/${ROOT_INCLUDE_PATH}
export CMAKE_PREFIX_PATH=${NUSYST_DQ_DIR}:${CMAKE_PREFIX_PATH}
export PKG_CONFIG_PATH=${NUSYST_DQ_DIR}:${PKG_CONFIG_PATH}
