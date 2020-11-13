#!/bin/bash

java -Xmx200g -jar beagle-18May20.d20.jar gt=WF1.vcf out=WF1Bout.vcf.gz

#accuracy test

java -Xmx50g -jar beagle-18May20.d20.jar gt=AcTestIn.vcf out=AcTestOut