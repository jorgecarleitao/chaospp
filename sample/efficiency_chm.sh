#!/usr/bin/env bash

export TARGET=efficiency_chm

for SAMPLES in 10 50 100 1000
do
export SAMPLES

for DIMENSION in 2 4 6 8 10
do
export DIMENSION

for MAX_ESCAPE_TIME in 2 4 8 16 20 24 28 32 36 40 48 56 64
do
export MAX_ESCAPE_TIME

### execute
if [ $(uname) == Linux ]; then
qsub -N ${TARGET}_${MAX_ESCAPE_TIME}_${MAP_DIMENSION} ./run.sh
else
cmake . && make ${TARGET} && ./${TARGET}
fi

done
done
done
