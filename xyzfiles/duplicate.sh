#!/bin/bash

cat 100_PBE_pure.xyz | awk '{ if (NF == 1){ print $0*2;} else if (NF == 4){ print $0; print $1"  "$2"  "$3"  "$4 + 14.41; } else{print "";}}' 
