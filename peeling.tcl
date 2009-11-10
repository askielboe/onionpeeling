proc peeling { args } {
	
	# Define parameters
	set t(1) 1.
	set norm(1) 1.
	set nH 1.
	set abundance 0.3
	set redshift 0.18
	
	# First we select a step-size and fitting range for alpha
	set alpha 1.
	set alphamin 0.1
	set alphamax 1.
	set alphastep 0.1
	
	# Calculate number of shells from lines in shells.dat
	set filenameradii "shells.dat"
	set fradii [open $filenameradii r]
	set i 1
	while {![eof $fradii]} {
		gets $fradii line 
		incr i
	}
	close $fradii
	set nshells [expr $i-1]
	
	# Fit outermost shell
	
	data fm_shell$nshells.fak
	
	# Then we loop through the fit for each value of alpha recording the reduced chi-square.
	set i 1
	for {set alpha $alphamin} {$alpha <= $alphamax} {set alpha [expr $alpha + $alphastep]} {
		
		# Calculate volumes for given alpha-value using sphvol
		set newvolumes [sphvol $alpha]
		array unset volumes
		array set volumes $newvolumes
		
		# Define mekal model from volume
		model 1:fakemek mekal*constant & $t(1) & $nH & $abundance & $redshift & & $norm(1) & $volumes(7,7) &
		freeze fakemek:7 ; # Freeze the volume
		fit 100 ; # Do the fit
		
		# Read parameters and errors from xspec output (using the tcloutr command)
		
		# TEMPERATURE
		for {set count 0} {[string index [tcloutr param fakemek:1] $count]!=" "} {incr count} {
			append param_temp($i) [string index [tcloutr param fakemek:1] $count]
			set n $count
		}
		for {set count [expr $count+1]} {[string index [tcloutr param fakemek:1] $count]!=" "} {incr count} {
			append param_temp_error($i) [string index [tcloutr param fakemek:1] $count]
		}
		# NORM (DENSITY?)
		for {set count 0} {[string index [tcloutr param fakemek:6] $count]!=" "} {incr count} {
			append param_norm($i) [string index [tcloutr param fakemek:6] $count]
		}
		for {set count [expr $count+1]} {[string index [tcloutr param fakemek:6] $count]!=" "} {incr count} {
			append param_norm_error($i) [string index [tcloutr param fakemek:6] $count]
		}
		# VOLUME
		for {set count 0} {[string index [tcloutr param fakemek:7] $count]!=" "} {incr count} {
			append param_vol($i) [string index [tcloutr param fakemek:7] $count]
		}
		for {set count [expr $count+1]} {[string index [tcloutr param fakemek:7] $count]!=" "} {incr count} {
			append param_vol_error($i) [string index [tcloutr param fakemek:7] $count]
		}
		
		# Get the reduced chi-square from the fit
		set prompt "XSPEC12>"
		
		# Get degrees of freedom
		set dof 1
		unset dof
		for {set count 0} {[string index [tcloutr dof] $count]!=" "} {incr count} {
			append dof [string index [tcloutr dof] $count]
		}
		
		# Get chi-square and compute reduced chi-square
		set chisquare($i) [expr [tcloutr stat]/$dof]
		
		set nalpha($i) $alpha
		
		puts "ALPHA = $alpha"
		
		incr i
	}
	
	# Find the alpha value corresponding to the best fit
	set minindex 1
	set chimin $chisquare(1)
	for {set n 2} {$n < $i} {incr n} {
		puts "Chi-square = $chisquare($n), chimin = $chimin"
		if {$chisquare($n) < $chimin} {
			set chimin $chisquare($n)
			set imin $n
			puts "NEW MINIMUM FOUND!"
		}
	}
	
	# Make a file with chisquare as a function of alpha for external plotting
	set fout [open chi_vs_alpha.txt w]
	for {set n 1} {$n < $i} {incr n} {
		puts $fout "$nalpha($n) $chisquare($n)"
	}
	close $fout
	
	set fout [open vol_vs_alpha.txt w]
	for {set n 1} {$n < $i} {incr n} {
		puts $fout "$nalpha($n) $param_vol($n) $param_vol_error($n)"
	}
	close $fout
	
	set fout [open temp_vs_alpha.txt w]
	for {set n 1} {$n < $i} {incr n} {
		puts $fout "$nalpha($n) $param_temp($n) $param_temp_error($n)"
	}
	close $fout
	
	set fout [open norm_vs_alpha.txt w]
	for {set n 1} {$n < $i} {incr n} {
		puts $fout "$nalpha($n) $param_norm($n) $param_norm_error($n)"
	}
	close $fout

	# Print the results
	puts "Best fit achieved for alpha = $nalpha($imin). Here the reduced chi-square was = $chisquare($imin)."
	puts "Parameters for best fit:"
	puts "Temperature kT = $param_temp($imin)"
	puts "Density nH = $param_rho($imin)"
	puts "Normalization (volume) HVAD ER ENHEDEN??? = $param_norm($imin)"
	
	# Wait for user to press any key
	puts "Press any key to start the onionpeeling..."
	gets stdin ans
	
	# # Print out results of the fits TO THE SCREEN
	# for {set n 1} {$n < $i} {incr n} {
	# 	puts "$n: Alpha = $nalpha($n), Reduced Chi-Square = $chisquare($n)"
	# }
	# 
	# # Print out results of the fits TO A FILE
	# set fout [open alphafits.txt w]
	# for {set n 1} {$n < $i} {incr n} {
	# 	puts $fout "$n: Alpha = $nalpha($n), Reduced Chi-Square = $chisquare($n)"
	# }
}



