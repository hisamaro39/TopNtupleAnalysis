#!/bin/bash

mkdir installed

. /afs/cern.ch/sw/lcg/external/gcc/4.6/x86_64-slc6/setup.sh
. /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.07_python2.7/x86_64-slc6-gcc46-opt/root/bin/thisroot.sh

PY=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc6-gcc46-opt

export PATH=${PY}/bin:$PATH
export LD_LIBRARY_PATH=${PY}/lib:$LD_LIBRARY_PATH

export PYTHONDIR=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc6-gcc46-opt
export PYTHONPATH=$PYTHONPATH:$ROOTSYS/lib
export LD_LIBRARY_PATH=$ROOTSYS/lib:$PYTHONDIR/lib:$LD_LIBRARY_PATH:/opt/rh/python27/root/usr/lib64

export PATH=$PWD/installed/bin:$PATH
export PYTHONPATH=$PWD/installed/lib/python2.7/site-packages:$PYTHONPATH
