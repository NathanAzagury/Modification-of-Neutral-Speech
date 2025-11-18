#!/Users/natha/OneDrive/Documents/Senior_Project/FinalCode

#set samples [readwav original.wav]	;#  get the wave (specify your own wav)
#set tags [tag]				;#  set input tags for TD-PSOLA

#puts "$samples samples, $tags tags"

proc longvowels samples {		;#  procedure to identify vowel segments
	set vowel 0
	set out [list]

	for {set i 0} {$i<$samples} {incr i 882} {
		set p [pitchy $i]
		if {$p>0.0 && $vowel==0} {
			set vowel 1
			set start $i
		}
		if {$p==0.0 && $vowel>0} {
			set vowel 0
			lappend out $start $i
		}
	}
	return $out
}

set ::minVowelStretch
set ::maxVowelStretch

foreach {s e} [longvowels $samples] {	;#  foreach start, end of vowel segment
	set diff [expr ($e-$s)/44100.0] ;#  compute duration in seconds
	if {$diff>0.1} {		;#  If at least a tenth of a second...

		set mult [expr {rand() * ($maxVowelStretch - $minVowelStretch) + $minVowelStretch}]
	        puts "Segment $s–$e duration=$diff  randomStretch=$mult"

		#set mult [expr $randomStretch/$diff]   ;# make each long vowel exactly 1s long
		 #                           ;# by stretching by a factor of 1/len

		for {set i $s} {$i<$e} {incr i 882} {
			setlong $i $mult             
		}
	}
}

#puts "[trimangle]"	;#  compute output tags, triangle window synthesis
#writewav new.wav		;#  write output file
