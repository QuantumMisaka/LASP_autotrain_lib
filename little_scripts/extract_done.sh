#!/bin/bash

mkdir LASP_DONE
for i in AuPdO*
do
	status=`grep elapse ${i}/lasp.out`
	if [[ ! -z $status ]]; then
		mv ${i} LASP_DONE
	fi
done


