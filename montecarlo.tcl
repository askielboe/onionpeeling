# Define the temperature profile, parameters: t0, tmin, rt, rcool a, acool, b, c
proc func_temp { R t0 tmin rt rcool a acool b c } {
	set x [expr pow($R/$rcool,$acool)]
	set tcool [expr ($x+$tmin/$t0)/($x+1)]
	set tr [expr pow($R/$rt,-$a)/pow(1+pow($R/$rt,$b),$c/$b)]
	set temp_model [expr $t0 * $tcool * $tr] 
	return $temp_model
}

set ntrials 1000
set npoints 100.
set rmax 10.

set diff_best 1000

set r1 0.3
set r2 0.6

# Do the Monte Carlo
for {set i 1} {$i <= $ntrials} {incr i} {
	# Choose random parameters
	set t0($i) [expr rand()+1.] ; # Value between 1 and 2
	set tmin($i) [expr rand()*5.] ; # Value between 0 and 5
	set rt($i) [expr rand()] ; # Value between 0 and 1
	set rcool($i) [expr rand()] ; # Value between 0 and 1
	set a($i) [expr (rand()-0.5)*4] ; # Value between -2 and 2
	set acool($i) [expr (rand()-0.5)*4] ; # Value between -2 and 2
	set b($i) [expr (rand()-0.5)*4] ; # Value between -2 and 2
	set c($i) [expr (rand()-0.5)*4] ; # Value between -2 and 2
	if {$rt($i) == 0 || $b($i) == 0} continue
	# Calculate points
	for {set point 1} {$point <= $npoints} {incr point} {
		set R [expr $point/$npoints*$rmax]
		set temp_real [expr 1.35*(pow($R/0.045,1.9)+0.45)/(pow($R/0.045,1.9)+1.)/pow(1+pow($R/0.6,2.),0.45)]
		set temp_model [func_temp $R $t0($i) $tmin($i) $rt($i) $rcool($i) $a($i) $acool($i) $b($i) $c($i)]
		set diff($point) [expr pow($temp_real-$temp_model,2.)]
	}
	# Calculate the reduced chi-square
	set diff_mean 0
	for {set point 1} {$point <= $npoints} {incr point} {
		set diff_mean [expr $diff_mean + $diff($point)] ; # First we sum the differences
	}
	set diff_mean [expr $diff_mean/$npoints] ; # Then divide by the number of points
	# If its better than the current best value, set new best fit.
	if {$diff_mean < $diff_best} {
		puts "New best fit found! (meandiff = $diff_mean). Iteration: $i."
		set diff_best $diff_mean
		set n $i
	}
}

puts "Best fit found at the $n'th iteration with meandiff = $diff_best."
puts "Parameters for best fit: $t0($n) $tmin($n) $rt($n) $rcool($n) $a($n) $acool($n) $b($n) $c($n)."

set fout [open mc_temp.gp w]
	puts $fout "f(x) = 1.35*((x/0.045)**1.9+0.45)/((x/0.045)**1.9+1.)/(1+(x/0.6)**2.)**0.45"
	puts $fout "g(x) = $t0($n) * ((x/$rcool($n))**$acool($n)+$tmin($n)/$t0($n))/((x/$rcool($n))**$acool($n)+1) * (x/$rt($n))**(-$a($n))/((1+(x/$rt($n))**$b($n))**($c($n)/$b($n)))"
	puts $fout "plot f(x),g(x)"
close $fout

exec gnuplot -persist mc_temp.gp