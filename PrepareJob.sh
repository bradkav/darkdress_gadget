#!/bin/bash

BASEDIR=/home/kavanagh/EMRI/code
TARGETDIR=/home/kavanagh/EMRI/code/sims/$1

RUNDIR=$TARGETDIR/run
OUTDIR=$TARGETDIR/out

echo "Setting up folder: $TARGETDIR"

mkdir -p $TARGETDIR
mkdir -p $RUNDIR
mkdir -p $OUTDIR

cp $BASEDIR/run/Gadget2 $RUNDIR

cp $BASEDIR/run/EMRI_HPC.param $RUNDIR/EMRI.param

cp $BASEDIR/run/EMRI1.dat $RUNDIR
cp $BASEDIR/RunEMRI_template.sh $TARGETDIR/RunEMRI_$1.sh

cd $TARGETDIR
#sed -i -- "s@RUN_DIRECTORY@$RUNDIR@g" RunPBH_$1.sh
sed -i -- "s@JOB_ID@$1@g" RunEMRI_$1.sh
