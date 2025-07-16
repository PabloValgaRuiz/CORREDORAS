#!/bin/bash
MAX_JOBS=12

FROM_RUNS=$1
TO_RUNS=$2

START=$3
END=$4

seq $FROM_RUNS $TO_RUNS | xargs -I {} -P $MAX_JOBS /home/paval/papers/CORREDORAS/pablo_python/.conda/bin/python transfer_parallel.py {} $START $END