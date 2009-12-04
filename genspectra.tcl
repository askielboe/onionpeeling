proc genspectra { args } {
	
	if {$args == ""} {
		puts "Syntax: genspectra <alpha-value>"
		break
	}
	
	# Read-in arguments RIGHT NOW ONLY ALPHA IS ALLOWED
	set alpha $args
	
	# FUNCTION & PARAMETER DEFINITIONS
	proc func_temp {R} {
		set temp [expr 1.35*pow($R/0.045,1.9)+0.45)/(pow($R/0.045,1.9)+1.)/pow(1+pow($R/0.6,2.),0.45)]
		#set temp [expr -16.5*$R+10.]
		# Martina temperature function:
		#set temp [expr -10.*$R+10.]
		return $temp
	}
	proc func_norm {R i} {
		if {$i < 4} {
			set norm 0.008
		} else {
			set norm [expr pow($R*100,-2.)+0.003]
		}
		# Martina density function:
		# set norm [expr pow(1/(1+pow($R/0.15,2.)),2.)]
		set upscalednorm [expr $norm*1000.] ; # UP-SCALING THE HYDROGEN NUMBER DENSITY OR NORM
		return $upscalednorm
	}
	
	set nH 1
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
		set norm($i) [func_norm $R($i) $i]
	
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
		model 1:fakemek mekal*constant & $t($iobs) & $nH & $abundance & $redshift & & $norm($iobs) & $volumes($iobs,$iobs) &
		if {[expr $nobs - $iobs] > 0} {
			set modelno 3
			for {set i [expr $iobs+1]} {$i <= $nobs} {incr i} {
				addcomp fakemek:$modelno mekal & $t($i) & $nH & $abundance & $redshift & & $norm($i) &
				incr modelno
				addcomp fakemek:$modelno constant & $volumes($i,$iobs) &
				incr modelno
			}
		}
		data 1666_3.pi
		fakeit & y & & fm_shell$iobs.fak & &
		puts "DATA FAKED!"
		puts "OUTPUT: fm_shell$iobs.fak"
	}
	# for {set i 1} {$i < 8} {incr i} {
	# 	puts $norm($i)
	# }
}