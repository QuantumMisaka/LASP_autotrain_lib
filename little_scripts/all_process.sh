#!/bin/bash
rm process.txt
touch process.txt
ROOTDIR=$(pwd)
for i in V*
do
	cd ${i}

	pwd >> ../process.txt
	../forset_process.sh >> ../process.txt
        cd $ROOTDIR
done
echo "Done!"

