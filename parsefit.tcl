proc parsefit { args } {
	# Read parameters and errors from xspec output (using the tcloutr command)

	# TEMPERATURE
	for {set count 0} {[string index [tcloutr param fakemek:1] $count]!=" "} {incr count} {
		append param_temp($i) [string index [tcloutr param fakemek:1] $count]
		set n $count
	}
	for {set count [expr $count+1]} {[string index [tcloutr param fakemek:1] $count]!=" "} {incr count} {
		append param_temp_error [string index [tcloutr param fakemek:1] $count]
	}
	# NORM (DENSITY?)
	for {set count 0} {[string index [tcloutr param fakemek:6] $count]!=" "} {incr count} {
		append param_norm [string index [tcloutr param fakemek:6] $count]
	}
	for {set count [expr $count+1]} {[string index [tcloutr param fakemek:6] $count]!=" "} {incr count} {
		append param_norm_error [string index [tcloutr param fakemek:6] $count]
	}
	# VOLUME
	for {set count 0} {[string index [tcloutr param fakemek:7] $count]!=" "} {incr count} {
		append param_vol [string index [tcloutr param fakemek:7] $count]
	}
	for {set count [expr $count+1]} {[string index [tcloutr param fakemek:7] $count]!=" "} {incr count} {
		append param_vol_error [string index [tcloutr param fakemek:7] $count]
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
	set chisquare [expr [tcloutr stat]/$dof]

	set nalpha $alpha

}