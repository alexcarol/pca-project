#!/bin/bash

molecules=( 1hba.parsed 2kai.parsed 2pka.parsed 4hhb.parsed 5pti.parsed )

for static in ${molecules[@]}; do
	for mobile in ${molecules[@]}; do
		runtime=`./ftdock -static $static -mobile $mobile | grep "Total time:" | cut -f2 -d ":"`
		echo -static $static -mobile $mobile: ${runtime}s
	done
done

