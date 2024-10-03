set inp [lindex $argv 0]
set top [lindex $argv 1]
set out [lindex $argv 2]
mol load gro $inp
set all_no_wats_inside [atomselect top "all and not same residue as waters and within 3 of protein" ]
set wats_inside [atomselect top "same residue as waters and within 3 of protein"]
set f [open wats.deleted a]
foreach i [lsort -unique [$wats_inside get resid]] {
	puts $f $i 
}
close $f
$all_no_wats_inside writegro $out

set length [llength [lsort -unique [$wats_inside get resid]]]
# Reduce SOL by $length
exec tail -1 $top > temp.length

set h [open temp.length r]
while { [gets $h data] >= 0 } {
	set sols $data
	set ns [lindex [split $sols " "] 0]
	set so [lindex [lreverse [split $sols " "]] 0]
}
close $h
if {$ns eq "SOL"} {
	set new_l [expr $so - $length]
} else {
	set new_l "Error"
}

file delete -force temp.length

if {$new_l ne "Error"} {
	exec head -n -1 $top > temp.length ; mv temp.length $top
	exec echo "SOL\t$new_l" > temp.length
	exec cat $top temp.length > temp.length2
	exec mv temp.length2 $top
	file delete temp.length
} else {
	exec head -n -1 $top > temp.length ; mv temp.length $top
	exec echo "$new_l\tRun at incorrect step -- run after solvation" > temp.length
	exec cat $top temp.length > temp.length2
	exec mv temp.length2 $top
	file delete temp.length
}

