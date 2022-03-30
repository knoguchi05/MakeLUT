#!/bin/sh

bsub -q 4h -n 2 -o Log.Fit root betaBatch.C 
