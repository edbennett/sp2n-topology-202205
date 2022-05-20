#!/usr/bin/env bash

set -eu

SCALES_FILE=$1
T0_FILE=$2
W0_FILE=$3

awk '{print "$",$3,"$&$",$4,"$&$",$5,"$&$",$6,"$ \\\\"}' ${SCALES_FILE} > ${T0_FILE}
sed -i -e '$s/\\\\//' ${T0_FILE}
awk '{print "$",$3,"$&$",$7,"$&$",$8,"$&$",$9,"$ \\\\"}' ${SCALES_FILE} > ${W0_FILE}
sed -i -e '$s/\\\\//' ${W0_FILE}
