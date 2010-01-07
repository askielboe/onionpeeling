# Define common procedures

# for reading out the reduced chi-square from a fit
proc getchisquare { args } {
	# Get degrees of freedom
	set dof 1
	unset dof
	for {set count 0} {[string index [tcloutr dof] $count]!=" "} {incr count} {
		append dof [string index [tcloutr dof] $count]
	}
	# Add a . at the end to ensure floating point results
	append dof "."
	# Get chi-square and compute reduced chi-square
	set chisquare [expr [tcloutr stat]/$dof]
	return $chisquare
}

proc getpar { args } {
	set par 1
	unset par
	for {set count 0} {[string index [tcloutr par $args] $count]!=" "} {incr count} {
		append par [string index [tcloutr par $args] $count]
	}
	return $par
}

proc geterror { args } {
	set error 1
	unset error
	set error [tcloutr sigma $args]
	# for {set count 0} {[string index [tcloutr sigma $args] $count]!=EOL} {incr count} {
	# 	append par [string index [tcloutr sigma $args] $count]
	# }
	return $error
}