#!/bin/bash
ROOTDIR=$(pwd)
ALL_JOBS=0
FIN_JOBS=0
BAD_JOBS=0
for i in para*
do
	if [[ -a ${i}/OSZICAR ]]
	then
		INFO=$(tail -n 2 ${i}/OSZICAR)
		RESULT=$(echo $INFO | grep "E0")
        	if [[ "$RESULT" != "" ]]
		then
			let FIN_JOBS+=1
		else
			let BAD_JOBS+=1
		fi
	else
		let BAD_JOBS+=1
	fi
	let ALL_JOBS+=1
done

echo "FINISHED/ALL: " $FIN_JOBS "/" $ALL_JOBS
if (( $FIN_JOBS == $ALL_JOBS  ))
then
	echo "Single_Calc TOTALLY DONE in " $ROOTDIR "!!!!"
else
	echo "Single_Calc NOT Fully done in " $ROOTDIR
fi






