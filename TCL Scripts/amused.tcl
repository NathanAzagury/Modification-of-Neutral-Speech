
set samples [readwav original.wav]         ;# get the audio file
set tags [tag]           ;#set input tags for TD-PSOLA

puts "$samples samples, $tags tags"

#set all global variables
set ::minVowelStretch 3.0
set ::maxVowelStretch 4.5
set ::silenceMult 0.85
set ::minPitch 1.20
set ::maxPitch 1.40
set ::rmsCeiling 1.35
set ::rmsFloor   700
set ::minRMS 1.3
set ::maxRMS 1.6

if 0 {
	source boostRMS.tcl
	puts "BoostRMS Processed"
	puts "[outtag] [triangle]"
#	writewav boostRMS.wav
#	puts "Audio written"
	set samples [commit]
	set tags [tag]
	puts "Commit Executed Successfully"
}

if 0 {
	source boostPitch.tcl
	puts "Boost Pitch Code Processed"
	puts "[outtag] [triangle]"
#	writewav varyingUptalk.wav
#	puts "Audio written"
	set samples [commit]
	set tags [tag]
	puts "Commit Executed Successfully"
}

if 1 {
	source long.tcl
	puts "Lengthening Code Processed"
	puts "[outtag] [trimangle]"
	set samples [commit]
#	writewav longVowel.wav
#	puts "Audio written"
	puts "Commit Executed Successfully"
	set tags [tag]
#	puts "Commit Executed Successfully"
}

if 0 {
	source longSilence.tcl
	puts "Lengthening Code Processed"
	puts "[outtag] [trimangle]"
#	writewav longSilence.wav
#	puts "Audio written"
	set samples [commit]
	puts "Commit Executed Successfully"
	set tags [tag]
#	puts "Commit Executed Successfully"
}

puts "Writing Output Audio"
puts "[outtag] [triangle]"         ;#compute output tags and triangle windo synthesis
writewav longVowel.wav                   ;#write output file
