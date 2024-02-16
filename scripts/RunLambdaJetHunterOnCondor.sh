#!/bin/bash
# -----------------------------------------------------------------------------
# 'RunLambdaJetHunterCondor.sh'
# Derek Anderson
# 02.13.2024
#
# Script to run the lambda-tagged jet finder
# module via condor.
# -----------------------------------------------------------------------------

# grab home directory and set up environment
export HOME=/sphenix/u/${LOGNAME}
source /opt/sphenix/core/bin/sphenix_setup.sh

# initialize no. of events to run and verbosity
nEvents=0
verbose=0

# create array of file lists
inputFiles="{"
for fileList in $@; do
  inputFiles+="\"${fileList}\","
done
inputFiles=${inputFiles::-1}
inputFiles+="}"

# run script
echo running: RunLambdaJetHunterOnCondor.sh $*
root -b -q "Fun4All_RunLambdaJetHunterOnCondor.C(${inputFiles}, $nEvents, $verbose)"
echo Script done

# end -------------------------------------------------------------------------
