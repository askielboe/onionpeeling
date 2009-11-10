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
	set alphamax 1
	set alphastep 0.1
	
	# Define procedure for reading out the reduced chi-square from a fit
	proc getchisquare { args } {
		# Get degrees of freedom
		set dof 1
		unset dof
		for {set count 0} {[string index [tcloutr dof] $count]!=" "} {incr count} {
			append dof [string index [tcloutr dof] $count]
		}
		# Get chi-square and compute reduced chi-square
		set chisquare [expr [tcloutr stat]/$dof]
		return $chisquare
	}
	
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
	
	# Then we loop through the fit for each value of alpha recording the reduced chi-square.
	set i 1
	for {set alpha $alphamin} {$alpha <= $alphamax} {set alpha [expr $alpha + $alphastep]} {
		
		# Calculate volumes for given alpha-value using sphvol
		set newvolumes [sphvol $alpha]
		array unset volumes
		array set volumes $newvolumes
		
		# Start by fitting the outermost shell (ishell = nshells)
		set ishell $nshells
		set meanchitemp 0
		data fm_shell$ishell.fak
		# Define mekal model from volume and name is after the shell (fakemekSHELLNUMBER)
		model mekal*constant & $t(1) & $nH & $abundance & $redshift & & $norm(1) & $volumes($ishell,$ishell) &
		freeze 7 ; # Freeze the volume
		fit 100 ; # Do the fit using atmost 100 iterations
		# Get the reduced chi-square af a function of iteration (i) and shell (ishell)
		set chisquare($i,$ishell) getchisquare

		# Then peel the rest of the shells while fitting with the given alpha-value
		for {set ishell [expr $nshells-1]} {$ishell >= 1} {set ishell [expr $ishell-1]} {
			# Copy the original spectra to the one we will subtract from
			exec cp fm_shell$ishell.fak sub_shell$ishell.fak
			data sub_shell$ishell.fak
			# Subtract the outer-lying spectra using mathpha
			for {set n $nshells} {$n > $ishell} {set n [expr $n-1]} {
				set volfrac [expr $volumes($n,$ishell)/$volumes($n,$n)]
				puts "Dividing volume ($n,$ishell) with volume ($n,$n)."
				puts "volfrac = $volumes($n,$ishell)/$volumes($n,$n) = $volfrac"
				puts "------------------------------------------------------------"
				puts "Subtracting $volfrac*sub_shell$n.fak from sub_shell$ishell.fak"
				puts "------------------------------------------------------------"
				if {$n == $nshells} {
					mathpha sub_shell$ishell.fak-$volfrac*fm_shell$n.fak R sub_shell.tmp sub_shell$ishell.fak NULL 0
				} else {
					mathpha sub_shell$ishell.fak-$volfrac*sub_shell$n.fak R sub_shell.tmp sub_shell$ishell.fak NULL 0	
				}
				exec rm sub_shell$ishell.fak
				exec mv sub_shell.tmp sub_shell$ishell.fak
				# # Plot the data to see whats going on
				# data sub_shell$ishell.fak
				# pl da
				# gets stdin dummyvar
				# if {$dummyvar == 1} {
				# 	break
				# }
			}
			# Do a fit with the volume corresponding to the shell of interest
			data sub_shell$ishell.fak
			response 1666_3.wrmf
			arf 1666_3.warf
			newpar 7 $volumes($ishell,$ishell)
			# freeze fakemek$ishell:7 ; # Freeze the volume
			puts "----------------------------------------------------------"
			puts "Fitting sub_shell$ishell.fak using volume $volumes($ishell,$ishell)"
			puts "----------------------------------------------------------"
			fit 100 ; # Do the fit using atmost 100 iterations
			# Get the reduced chi-square af a function of iteration (i) and shell (ishell)
			set meanchitemp [expr $meanchitemp + [getchisquare]]
		}
		# Save alpha value as a function of iteration (i) and calculate mean chi-square
		set nalpha($i) $alpha
		set meanchi($i) [expr $meanchitemp/$nshells]
		incr i
	}
	
	# Print mean chi's and alpha's
	# Make a file with meanchisquare as a function of alpha for external plotting
	set fout [open chi_vs_alpha.txt w]
	for {set n 1} {$n < $i} {incr n} {
		puts "Mean reduced chi-square = $meanchi($n) for Alpha = $nalpha($n)"
		puts $fout "$nalpha($n) $meanchi($n)"
	}
	close $fout
	
	# # Find the alpha value corresponding to the best fit
	# set minindex 1
	# set chimin $chisquare(1,1)
	# for {set n 2} {$n < $i} {incr n} {
	# 	puts "Chi-square = $chisquare($n), chimin = $chimin"
	# 	if {$chisquare($n) < $chimin} {
	# 		set chimin $chisquare($n)
	# 		set imin $n
	# 		puts "NEW MINIMUM FOUND!"
	# 	}
	# }
	# 
	# 
	# set fout [open vol_vs_alpha.txt w]
	# for {set n 1} {$n < $i} {incr n} {
	# 	puts $fout "$nalpha($n) $param_vol($n) $param_vol_error($n)"
	# }
	# close $fout
	# 
	# set fout [open temp_vs_alpha.txt w]
	# for {set n 1} {$n < $i} {incr n} {
	# 	puts $fout "$nalpha($n) $param_temp($n) $param_temp_error($n)"
	# }
	# close $fout
	# 
	# set fout [open norm_vs_alpha.txt w]
	# for {set n 1} {$n < $i} {incr n} {
	# 	puts $fout "$nalpha($n) $param_norm($n) $param_norm_error($n)"
	# }
	# close $fout
	# 
	# # Print the results
	# puts "Best fit achieved for alpha = $nalpha($imin). Here the reduced chi-square was = $chisquare($imin)."
	# puts "Parameters for best fit:"
	# puts "Temperature kT = $param_temp($imin)"
	# puts "Density nH = $param_rho($imin)"
	# puts "Normalization (volume) HVAD ER ENHEDEN??? = $param_norm($imin)"
	# 
	# # Wait for user to press any key
	# puts "Press any key to start the onionpeeling..."
	# gets stdin ans
	
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



