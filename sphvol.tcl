# Define constants
set filenamein "shells.dat"
set filenameout "volumes.out"
set Pi 3.1415926535897932385

# Define procedures
proc Ve {rc a b} {
	set Pi 3.1415926535897932385
	if {$rc > $a} {
		set rc $a
	}
	set k [expr 1.-pow($rc,2)/pow($a,2)]
	set ve [expr 1./3.*$Pi*pow($a,2)*$b*(2.-3.*sqrt($k) + pow(sqrt($k),3))]
	return $ve
}
proc Vc {rc a b} {
	set Pi 3.1415926535897932385
	if {$rc > $a} {
		set rc $a
	}
	set k [expr 1.-pow($rc,2)/pow($a,2)]
	set vc [expr 2.*$Pi*pow($rc,2)*$b*sqrt($k)]
	return $vc
}
proc Vcut {rc a b} {
	set vcut [expr 2.*[Ve $rc $a $b] + [Vc $rc $a $b]]
	return $vcut
}
proc Vshell {rc a b amin bmin} {
	set vshell [expr [Vcut $rc $a $b] - [Vcut $rc $amin $bmin]]
	return $vshell
}

# Read in list with n's and m's.
set fin [open $filenamein r]
set i 1
while {![eof $fin]} {
	gets $fin line ; # Read line until we encounter a space (" ") and put into n
	for {set count 0} {[string index $line $count]!=" "} {incr count} {
		append n($i) [string index $line $count]
	}
	incr count ; # Skip the space.
	while {[string index $line $count]>=-1} { ; # And read rest of line into m.
		append m($i) [string index $line $count]
		incr count
	}
	
	# Define a(n), right now a = n
	set a($i) [expr 1*$n($i)]
	
	# Define alpha(r) and beta(r)
	set alpha($i) [expr 0.5*$n($i)]
	set beta($i) [expr 1.*1.]
	
	# Define b(a) = a*(beta + alpha * a)
	set b($i) [expr $a($i)*($beta($i) + $alpha($i)*$a($i))]
	
	if [eof $fin] break ;# otherwise loops one time too many
	incr i
}
close $fin

# Calculate volume elements velement(n,m)
set in 1
while {$in <= $i} {
	set im 1
	while {$im <= $i} {
		# Don't try and calculate if any var is 0 or if rc > a.
		if {$m($im) == 0 || $a($in) == 0 || $b($in) == 0 || $m($im) > $a($in)} {
			set velement($in,$im) 0 
		} elseif {$in == 1 && $im == 1} { ; # Don't subtract in the innermost shell.
			set velement($in,$im) [Vcut $m(1) $a(1) $b(1)]
		} elseif {$in == 1} { ; # For n = 1 only subtract inner cuts.
			set velement($in,$im) [expr [Vcut $m($im) $a($in) $b($in)] - [Vcut $m([expr $im-1]) $a($in) $b($in)]]
		} elseif {$im == 1} { ; # For m = 1 only subtract inner spheroids.
			set velement($in,$im) [Vshell $m($im) $a($in) $b($in) $a([expr $in-1]) $b([expr $in-1])]
		} else { ; # For anything else, subtract both.
			set velement($in,$im) [expr [Vshell $m($im) $a($in) $b($in) $a([expr $in-1]) $b([expr $in-1])] - [Vshell $m([expr $im-1]) $a($in) $b($in) $a([expr $in-1]) $b([expr $in-1])]]
		}
		incr im
	}
	incr in
}

# Write output to file $filenameout
set fout [open $filenameout w]
#puts $fout "a rc V(n,m) \n"
set in 1
while {$in <= $i} {
	set im 1
	while {$im <= $i} {	
		if {$velement($in,$im) != 0} {
			puts $fout "$in $im $velement($in,$im) \n"
		}
		incr im
	}
	incr in
}
close $fout
