#!/bin/bash

TOTAL=`ls -d Au*Pd*O* | wc -l`
CALC=`ls Au*Pd*O*/lasp.out | wc -l`

echo "TOTAL JOBS: $TOTAL"
echo "CALC-ED JOBS: $CALC"
