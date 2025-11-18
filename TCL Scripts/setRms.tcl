#!/Users/natha/OneDrive/Documents/Senior_Project/FinalCode

set samples [readwav "original.wav"]
set tags [tag]
set WIDDY 882
set totalFrames [expr {int(ceil(double($samples) / $WIDDY))}]

# Loop only within valid frame range (from 0 to totalFrames-1)
for {set i 0} {$i < $totalFrames} {incr i} {
    if {$i < 50} {
        if {[catch {setrms $i 0.8}]} break
    } elseif {$i < 100} {
        if {[catch {setrms $i 1.0}]} break
    } elseif {$i < 150} {
        if {[catch {setrms $i 1.4}]} break
    } else {
        if {[catch {setrms $i 1.8}]} break
    }
}

puts "[outtag] [triangle]"
writewav new.wav
