#!/bin/bash
MAX_JOBS=20

FROM_RUNS=$1
TO_RUNS=$2

START=100
END=140

seq $FROM_RUNS $TO_RUNS | xargs -I {} -P $MAX_JOBS /home/paval/papers/CORREDORAS/pablo_python/.conda/bin/python basa_transfer_parallel.py {} $START $END