#!/Users/natha/OneDrive/Documents/Senior_Project/FinalCode

set samples [readwav original.wav]         ;# get the audio file
set tags [tag]           ;#set input tags for TD-PSOLA

puts "$samples samples, $tags tags"

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

set multiplier 1.8
set pitchMult 0.9

foreach {s e} [longvowels $samples] {	;#  foreach start, end of vowel segment
	set diff [expr ($e-$s)/44100.0] ;#  compute duration in seconds
	if {$diff>0.1} {		;#  If at least a tenth of a second...
		puts "$s $e $diff"	;#  code below to uptalk

		for {set i $s} {$i<$e} {incr i 882} {
			set frac [expr ($i-$s)*1.0/($e-$s)]  ;# position in vowel
			set mult [expr $pitchMult+$multiplier*$frac]        ;# range 1.0...1.5
			setpitch $i $mult                    ;# update pitch
		}
	}
}


puts "[outtag] [triangle]"         ;#compute output tags and triangle windo synthesis
writewav uptalk.wav                   ;#write output file

