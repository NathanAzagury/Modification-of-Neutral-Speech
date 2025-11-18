#!/Users/natha/OneDrive/Documents/Senior_Project/FinalCode

#set samples [readwav original.wav]	;#  get the wave (specify your own wav)
#set tags [tag]				;#  set input tags for TD-PSOLA

#puts "$samples samples, $tags tags"

proc silences samples {		;#  procedure to identify silences
	set silence 0
	set out [list]

	for {set i 0} {$i<$samples} {incr i 882} {
		set p [pitchy $i]
		set r [RMS $i]
		if {$p == 0.0 && $r < 700 && $silence==0} {
			set silence 1
			set start $i
#			puts "1. $r rms, $silence silence"

		}
		if {($p > 0.0 || $r >= 700) && $silence>0} {
			set silence 0
			lappend out $start $i
#			puts "2. $r rms, $silence silence"

		}
	}
	return $out
}

set ::silenceMult

foreach {s e} [silences $samples] {	;#  foreach start, end of vowel segment
	set diff [expr ($e-$s)/44100.0] ;#  compute duration in seconds
	if {$diff>0.1} {		;#  If at least a tenth of a second...
		puts "$s $e $diff"	;#  code below to uptalk
		set mult [expr {$::silenceMult * (1.0 + ((rand() - 0.5) * 0.15))}]
		for {set i $s} {$i<$e} {incr i 882} {
			set frame [expr {$i / 882}]
			setlong $i ::silenceMult
		}
	}
}

#puts "[trimangle]"	;#  compute output tags, triangle window synthesis
#writewav new.wav		;#  write output file
