#!/bin/bash

#make shell script
python Transformations_p1.py

#cd to temporary directory to run script
chmod +x transform.sh
./transform.sh

#compile results
python Transformations_p2.py