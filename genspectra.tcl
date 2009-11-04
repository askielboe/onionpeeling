proc genspectra { args } {
	
	# Read-in arguments RIGHT NOW ONLY ALPHA IS ALLOWED
	set alpha $args
	
	# FUNCTION & PARAMETER DEFINITIONS
	proc func_temp {R} {
		set temp [expr -16.5*$R+10.]
		return $temp
	}
	proc func_rho {R i} {
		if {$i < 4} {
			set rho 0.008
		} else {
			set rho [expr pow($R*100,-2.)+0.003]
		}
		return $rho
	}
	set abundance 0.3
	set redshift 0.18

	# Script for generating xspec code
	# Takes as input the output from sphvol.tcl and the chosen cylinder-shell we would like to observe.

	set filenameradii "shells.dat"

	# Read radii and define parameters
	set e 2.71828183

	set fradii [open $filenameradii r]
	set i 1
	while {![eof $fradii]} {
		gets $fradii line ; # Read line until we encounter a space (" ") and put into n
		for {set count 0} {[string index $line $count]!=" "} {incr count} {
			append R($i) [string index $line $count]
		}
		# Define temperature as a function of radius T(R)
		set t($i) [func_temp $R($i)]
	
		# Define density as a function of radius rho(R)
		set rho($i) [func_rho $R($i) $i]
	
		incr i
	}
	close $fradii
	
	# Set number of observations = number of shells = number of lines in filenameradii
	set nobs [expr $i-1]

	# Generate volumes with sphvol
	set newvolumes [sphvol $alpha]
	array unset volumes
	array set volumes $newvolumes

	# Define procedure to write addcomp to models as a function of volume, shell and temperature
	# Syntax: addcomp [Model Name:]<new component number> componentName
	for {set iobs 1} {$iobs <= $nobs} {incr iobs} {
		# Initiate the model
		model 1:fakemek mekal & $t(1) & $rho(1) & $abundance & $redshift & & $volumes($iobs,$iobs) &
		if {[expr $nobs - $iobs] > 0} {
			set modelno 2
			for {set i [expr $iobs+1]} {$i <= $nobs} {incr i} {
				addcomp fakemek:$modelno mekal & $t($i) & $rho($i) & $abundance & $redshift & & $volumes($i,$iobs) &
				incr modelno
			}
		}
		data 1666_3.pi
		fakeit & y & & fm_shell$iobs.fak & &
	}
}