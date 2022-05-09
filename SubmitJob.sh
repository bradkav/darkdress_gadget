#!/bin/bash

cd /home/kavanagh/EMRI/code

./PrepareJob.sh $1

cd /home/kavanagh/EMRI/code/sims/$1
qsub RunEMRI_$1.sh
