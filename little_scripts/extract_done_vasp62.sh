#!/bin/bash

mkdir DFT_DONE
mkdir DFT_EXCEED
DIRTYPE=para
LIMIT=200
for i in ${DIRTYPE}*
do
	if [[ ! -z $(grep E0 ${i}/OSZICAR) ]]; then
		NELMs=$(cat ${i}/OSZICAR | wc -l)
		if [[ $NELMs -le $LIMIT ]]
		then
			mv ${i} DFT_DONE
		else
			mv ${i} DFT_EXCEED
		fi
	fi
done


