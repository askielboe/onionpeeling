set a 0.291982
set b 0.145991
set k 2
set Pi 3.1415

set ve [expr 1/3.*$Pi*pow($a,2)*$b*(2.-3.*sqrt($k) + pow(sqrt($k),3))]

puts $ve