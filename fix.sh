#!/bin/bash
PROJDIR=/projects/kg98/kristina/GenofCog/datadir/derivatives
SCRIPTDIR=/projects/kg98/kristina/code/
SUBJ=sub-000
FEATDIR=$PROJDIR/$SUBJ/prepro.feat

module load fsl/5.0.9
module load fix/1.064
module load R/3.3.1
module load matlab/r2016a



echo -e "\nRunning FIX: $SUBJ \n"

/usr/local/fix/1.064/bin/fix $FEATDIR /projects/kg98/kristina/GenofCog/training/Training.RData 20 -m

echo -e "\nDone! \n"


