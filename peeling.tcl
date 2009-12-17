proc peeling { args } {
	
	# Define parameters
	set t(1) 1.
	set norm(1) 1.
	set nH 1.
	set abundance 0.3
	set redshift 0.18
	
	# First we select a step-size and fitting range for alpha
	set alpha 1.
	set alphamin 0.4
	set alphamax 0.9
	set alphastep 0.01
	
	# change the weighting to provide a better estimate of the variance in the small number limit.
	weight churazov
	
	# Calculate number of shells from lines in shells.dat
	set filenameradii "shells.dat"
	set fradii [open $filenameradii r]
	set i 1
	while {![eof $fradii]} {
		gets $fradii line
		for {set count 0} {[string index $line $count]!=" "} {incr count} {
			append radius($i) [string index $line $count]
		}
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
		data 1666_$ishell.pi
		# Define mekal model from volume and name is after the shell (fakemekSHELLNUMBER)
		model mekal*constant & $t(1) & $nH & $abundance & $redshift & & $norm(1) & $volumes($ishell,$ishell) &
		freeze 7 ; # Freeze the volume
		fit 100 ; # Do the fit using atmost 100 iterations
		# ignore 1:0-20
		# ignore 1:450-2000		
		# Get the reduced chi-square af a function of iteration (i) and shell (ishell)
		set chisquare($i,$ishell) [getchisquare]
		set fit_temp($i,$ishell) [getpar 1]
		set fit_normval($i,$ishell) [getpar 6]

		# Then peel the rest of the shells while fitting with the given alpha-value
		for {set ishell [expr $nshells-1]} {$ishell >= 1} {set ishell [expr $ishell-1]} {
			# Copy the original spectra to the one we will subtract from
			exec cp 1666_$ishell.pi sub_shell$ishell.fak
			data sub_shell$ishell.fak
			# Subtract the outer-lying spectra using mathpha
			for {set n $nshells} {$n > $ishell} {set n [expr $n-1]} {
				set volfrac [expr $volumes($n,$ishell)/$volumes($n,$n)]
				puts "Dividing volume ($n,$ishell) with volume ($n,$n)."
				puts "volfrac = $volumes($n,$ishell)/$volumes($n,$n) = $volfrac"
				puts "------------------------------------------------------------"
				puts "Subtracting $volfrac*sub_shell$n.fak from sub_shell$ishell.fak"
				puts "------------------------------------------------------------"
				# Exposuretime: 4.18E4, properr=’yes’ ?
				if {$n == $nshells} {
					mathpha expr=sub_shell$ishell.fak-$volfrac*1666_$n.pi outfil=!sub_shell.tmp units=R properr=yes areascal=sub_shell$ishell.fak exposure=sub_shell$ishell.fak ncomments=0
				} else {
					mathpha expr=sub_shell$ishell.fak-$volfrac*sub_shell$n.fak outfil=!sub_shell.tmp units=R properr=yes areascal=sub_shell$ishell.fak exposure=sub_shell$ishell.fak ncomments=0
				}
				exec mv sub_shell.tmp sub_shell$ishell.fak
			}
			# Do a fit with the volume corresponding to the shell of interest
			data sub_shell$ishell.fak
			response 1666_3.wrmf
			arf 1666_3.warf
			newpar 7 $volumes($ishell,$ishell)
			freeze 7 ; # Freeze the volume
			# ignore 1:0-20
			# ignore 1:450-2000
			puts "----------------------------------------------------------"
			puts "Fitting sub_shell$ishell.fak using volume $volumes($ishell,$ishell)"
			puts "----------------------------------------------------------"
			fit 100 ; # Do the fit using atmost 100 iterations
			# Get the reduced chi-square af a function of iteration (i) and shell (ishell)
			set meanchitemp [expr $meanchitemp + [getchisquare]]
			# Get parameters
			set chisquare($i,$ishell) [getchisquare]
			set fit_temp($i,$ishell) [getpar 1]
			set fit_norm($i,$ishell) [getpar 6]
		}
		# Save alpha value as a function of iteration (i) and calculate mean chi-square
		set nalpha($i) $alpha
		set meanchi($i) [expr $meanchitemp/$nshells]
		incr i
	}
	
	# Find the alpha value corresponding to the best fit
	set minindex 1
	set imin 1
	set chimin $meanchi(1)
	for {set n 2} {$n < $i} {incr n} {
		puts "Meanchi = $meanchi($n), chimin = $chimin"
		if {$meanchi($n) < $chimin} {
			set chimin $meanchi($n)
			set imin $n
		}
	}
	
	# Print the results
	puts "=================================== RESULTS ===================================="
	puts "Best fit achieved for alpha = $nalpha($imin). Here the reduced chi-square was = $meanchi($imin)."
	
	# Print mean chi's and alpha's
	# Make a file with meanchisquare as a function of alpha for external plotting
	set fout [open chi_vs_alpha.txt w]
	for {set n 1} {$n < $i} {incr n} {
		# puts "Mean reduced chi-square = $meanchi($n) for Alpha = $nalpha($n)"
		puts $fout "$nalpha($n) $meanchi($n)"
	}
	close $fout
	
	# Print all chi-squares as a function of alpha and shellno
	set fout [open chisquare.txt w]
	puts $fout "Colums are for different shells, starting with the innermost." 
	puts $fout "Rows are for different alpha-values."
	for {set n 1} {$n < $i} {incr n} {
		puts $fout "$nalpha($n) $chisquare($n,1) $chisquare($n,2) $chisquare($n,3) $chisquare($n,4) $chisquare($n,5) $chisquare($n,6)"
	}
	close $fout	

	# Write temperature profile for best fit to file
	set fout [open temp_vs_radius.txt w]
	for {set n 1} {$n < $nshells} {incr n} {
		#puts "$radius($n) $fit_temp($imin,$n)"
		puts $fout "$radius($n) $fit_temp($imin,$n)"
	}
	close $fout
	
	# Write density (norm) profile for best fit to file
	set fout [open norm_vs_radius.txt w]
	for {set n 1} {$n < $nshells} {incr n} {
		#puts "$radius($n) $fit_norm($imin,$n)"
		puts $fout "$radius($n) $fit_norm($imin,$n)"
	}
	close $fout
	
	# Write temperature profile for best fit to file
	set fout [open shellno_vs_chisquare.txt w]
	for {set n 1} {$n < $nshells} {incr n} {
		puts "$radius($n) $chisquare($imin,$n)"
		puts $fout "$radius($n) $chisquare($imin,$n)"
	}
	
	# Plot reduced chi-square as a function of alpha using gnuplot as defined in plot.gp
	exec gnuplot -persist plot_chi.gp
	# exec gnuplot -persist plot_temp.gp
}