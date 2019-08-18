#!/bin/bash

basedir=PPWWDD

source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530
cd /user/jdeclerc/CMSSW_8_0_30/src ; eval `scram runtime -sh` ; cd - >/dev/null
export X509_USER_PROXY=/user/$USER/x509up_u$(id -u $USER)
export X509_USER_PROXY=/tmp/x509up_u20641
voms-proxy-init -pwstdin < /user/jdeclerc/CMSSW_8_0_30/src/STEP1_Sexaq/password
#tmpdir=/user/lowette/`date +"%s"`_$RANDOM
tmpdir=$TMPDIR
mkdir -p $tmpdir
cp $basedir/SUS-RunIISummer16DR80Premix-00068_IDDD_cfg.py $tmpdir
cd $tmpdir
cmsRun SUS-RunIISummer16DR80Premix-00068_IDDD_cfg.py > JobIDDD.log 2>&1

# move job log back
mv JobIDDD.log $basedir
#mv *.std* $basedir # probably useless and done by job scheduler
# make sure folder on storage element exists
gfal-mkdir srm://maite.iihe.ac.be:8443/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Step1/`basename $basedir` 2>/dev/null
# put root output on storage element
for f in `ls *root` ; do
  gfal-copy file://$tmpdir/$f srm://maite.iihe.ac.be:8443/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Step1/`basename $basedir`/$f
  #rm $f
done
