#!/bin/bash
#### Description: This is a wrapper for excavator2 scripts
#. /venv/bin/activate
if [[ "$1" = "TargetPerla.pl" ]]
then
  TargetPerla.pl "${@:2}"
elif [[ "$1" = "EXCAVATORDataPrepare.pl" ]]
then
  EXCAVATORDataPrepare.pl "${@:2}"
elif [[ "$1" = "EXCAVATORDataAnalysis.pl" ]]
then
  EXCAVATORDataAnalysis.pl "${@:2}"
else
  echo "This is excavator2. Either run TargetPerla.pl or EXCAVATORDataPrepare.pl or EXCAVATORDataAnalysis.pl "
fi

