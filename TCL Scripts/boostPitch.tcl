#!/Users/natha/OneDrive/Documents/Senior_Project/FinalCode

#set samples [readwav original.wav]         ;# get the audio file
#set tags [tag]           ;#set input tags for TD-PSOLA

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

set ::minPitch
set ::maxPitch

set vowelIndex 0
foreach {s e} [longvowels $samples] {
    if {$e > $samples} {set e $samples}
    set diff [expr {($e - $s)/44100.0}]
    if {$diff > 0.1} {
        # one random boost per vowel (same logic as RMS)
        set pitchBoost [expr {rand() * ($::maxPitch - $::minPitch) + $::minPitch}]
        puts "Vowel $vowelIndex uses pitchBoost=$pitchBoost"
        incr vowelIndex

        for {set i $s} {$i < $e && $i < $samples} {incr i 882} {
            setpitch $i $pitchBoost
        }
    }
}

#puts "[outtag] [triangle]"         ;#compute output tags and triangle windo synthesis
#writewav new.wav                   ;#write output file
