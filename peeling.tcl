proc peeling { args } {
	
	# Define parameters
	set t(1) 1.
	set rho(1) 1.
	set abundance 0.3
	set redshift 0.18
	
	# First we select a step-size and fitting range for alpha
	set alpha 1.
	set alphamin 0.1
	set alphamax 1.
	set alphastep 0.1
	
	# Fit outermost shell
	
	data fm_shell7.fak
	
	# Then we loop through the fit for each shell recording the reduced chi-square
	set i 1
	for {set alpha $alphamin} {$alpha <= $alphamax} {set alpha [expr $alpha + $alphastep]} {
		# Calculate volumes for given alpha-value
		set newvolumes [sphvol $alpha]
		array unset volumes
		array set volumes $newvolumes
		
		# Define mekal model from volume
		model 1:fakemek mekal & $t(1) & $rho(1) & $abundance & $redshift & & $volumes(7,7) &
		thaw fakemek:2 ; # Un-freeze density
		freeze fakemek:6 ; # Freeze volume
		fit ; # Do the fit
		
		# Read-out parameters
		# TEMPERATURE
		for {set count 0} {[string index [tcloutr param fakemek:1] $count]!=" "} {incr count} {
			append param_temp($i) [string index [tcloutr param fakemek:1] $count]
		}
		# DENSITY
		for {set count 0} {[string index [tcloutr param fakemek:2] $count]!=" "} {incr count} {
			append param_rho($i) [string index [tcloutr param fakemek:2] $count]
		}
		# VOLUME
		for {set count 0} {[string index [tcloutr param fakemek:6] $count]!=" "} {incr count} {
			append param_norm($i) [string index [tcloutr param fakemek:6] $count]
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
		
		incr i
	}
	
	# Find the alpha value corresponding to the best fit
	set minindex 1
	set chimin $chisquare(1)
	for {set n 2} {$n < $i} {incr n} {
		if {$chisquare($n) < $chimin} {
			set chimin $chisquare($n)
			set imin $n
		}
	}

	# Print the results
	puts "Best fit achieved for alpha = $nalpha($imin). Here the reduced chi-square was = $chisquare($imin)."
	puts "Parameters for best fit:"
	puts "Temperature kT = $param_temp($imin)"
	puts "Density nH = $param_temp($imin)"
	puts "Normalization (volume) HVAD ER ENHEDEN??? = $param_norm($imin)"
	
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



