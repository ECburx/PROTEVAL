set pdb [lindex $argv 0]
mol load pdb $pdb
set sel [atomselect top all]
$sel moveby [vectinvert [measure center $sel]]
$sel writepdb [lindex $argv 1]
$sel delete
