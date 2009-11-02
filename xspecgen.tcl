# Script for generating xspec code
# Takes as input the output from sphvol.tcl

set filenamein "volumes.out"
set filenameout "fakemekal.xspec"

set fin [open $filenamein r]

close $fin

set fout [open $filenameout w]

close $fout