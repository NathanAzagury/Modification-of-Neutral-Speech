#!/Users/natha/OneDrive/Documents/Senior_Project/FinalCode

set samples [readwav original.wav]	;#  get the wave (specify your own wav)
set tags [tag]				;#  set input tags for TD-PSOLA

puts "$samples samples, $tags tags"

set mainpitch 0		;# will set to pitch of first pitch in file

for {set i 0} {$i<$samples-882} {incr i 882} {
	set time [expr int((2.0*$i)/44100.0)] ;#  time in half seconds

	set p [pitchy $i] 
	if {$p==0} continue
	if {$mainpitch==0} {set mainpitch $p; puts "mainpitch $mainpitch"}

	set factor ($p*1.0)/$mainpitch	 ;# makes everything monotone
	
	set note [expr $time%3]
	set note [lindex {0 4 7} $note]
	
	set ftwo [expr 2**(($note)/12.)]	;# convert note to frequency

	setpitch $i [expr $factor*$ftwo]
}

puts "[outtag] [triangle]"
#puts "doing pitch to [dopitch] samples"
writewav autotune.wav		;#  write output file

#exec open -a "Audacity" new.wav
