#!/Users/natha/OneDrive/Documents/Senior_Project/FinalCode

proc longvowels samples {
    set vowel 0
    set out [list]
    for {set i 0} {$i < $samples} {incr i 882} {
        set p [pitchy $i]
        if {$p > 0.0 && $vowel == 0} {
            set vowel 1
            set start $i
        }
        if {$p == 0.0 && $vowel > 0} {
            set vowel 0
            lappend out $start $i
        }
    }
    return $out
}

set ::minRMS
set ::maxRMS
set ::RMSMult [list]
set boosted 0
set vowelIndex 0
foreach {s e} [longvowels $samples] {
    if {$e > $samples} {set e $samples}
    set diff [expr {($e - $s) / 44100.0}]
    if {$diff > 0.05} {

        lappend ::RMSMult [expr {rand() * ($maxRMS - $minRMS) + $minRMS}]
        set RMSBoost [lindex $::RMSMult $vowelIndex]
	incr vowelIndex

	puts "BoostRMS: boost=$RMSBoost  ceiling=$::rmsCeiling  floor=$::rmsFloor"

    	if {$RMSBoost > $::rmsCeiling} {set dynamicBoost $::rmsCeiling}
	for {set i $s} {$i < $e} {incr i 882} {
            set frame [expr {$i / 882}]
            set p [pitchy $i]
            set r [RMS $i]
            if {$p == 0.0 || $r < $::rmsFloor} {
                continue
            }
	set RMSBoost [expr {1.0 + ($RMSBoost - 1.0)}]

        setrms $frame $RMSBoost
        }
    incr boosted
    }
}

puts "BoostRMS: applied to $boosted frames."
