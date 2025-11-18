# Expressive Speech Modification Engine
*A C + TCL software pipeline for pitch shifting, RMS manipulation, vowel lengthening, and expressive speech synthesis.*

This project implements a standalone C program that embeds a TCL interpreter and exposes custom audio-processing commands for modifying speech. The engine performs pitch detection, RMS analysis, vowel/silence stretching, and time-domain resynthesis, while higher-level control logic is written entirely in TCL.

This architecture separates:

- C (processing): low-level DSP, WAV I/O, pitch/RMS calculations, time-warping  
- TCL (control): expressive transformations (uptalk, amused speech, RMS boost, vowel stretching)

The result is a flexible processing engine where DSP algorithms stay fast and consistent in C, while expressive behavior can be scripted easily in TCL.

---

# Features

## C Engine
- Reads and writes WAV files (44.1 kHz mono/stereo)
- Performs:
  - Pitch detection (`pitchy`)
  - RMS energy measurement (`RMSAt`)
  - Harmonic removal (`defund`)
  - Fundamental-only resynthesis (`fund`)
  - Noise injection
  - Sample-level edits (`sampleAt`, `samplePut`, etc.)
  - Marker/segment tagging (`tag`, `outtag`)
  - Multiple time-warp synthesis modes (`triangle`, `trimangle`)
- Allows dynamic per-frame adjustments:
  - Pitch multipliers  
  - RMS multipliers  
  - Length multipliers (vowel/silence stretching)

## TCL Scripting Layer
- Loads/controls processing pipeline
- Defines global parameter arrays:
  - `pitches(frame)`
  - `rmses(frame)`
  - `speeds(frame)`
- Implements high-level effects such as:
  - Uptalk
  - RMS boost
  - Pitch boost
  - Vowel lengthening
  - Silence stretching
  - Amused speech
  - Sine-based contours
- Supports multi-stage pipelines via `masterTcl.tcl`

---

# Repository Structure

```
C_Engine_With_TCL_Interface/
    SpeechModification.c     # Full DSP engine and TCL command registration

TCL Scripts/
    amused.tcl
    boostPitch.tcl
    boostRMS.tcl
    long.tcl
    longSilence.tcl
    masterTcl.tcl
    music.tcl
    setRms.tcl
    sinUptalk.tcl
    uptalk.tcl
    uptalkSin.tcl
    varyingUptalk.tcl

Examples/
    original.wav
    uptalk.wav
    longVowel.wav
    spongbob.wav
    autotune.wav
```

---

# Compiling the Engine

Compile with:

```bash
gcc C_Engine_With_TCL_Interface/SpeechModification.c -o speechengine -lm -ltcl
```

You may need to link explicitly to your system’s TCL version:

```
-ltcl8.6
```

---

# Running the Engine

### Run a standalone TCL script
```bash
./speechengine TCL\ Scripts/uptalk.tcl
```

### Run a master pipeline script
```bash
./speechengine TCL\ Scripts/masterTcl.tcl
```

---

# Standalone Scripts vs Master Scripts

Your TCL scripts fall into two categories:

---

## 1. Standalone Scripts (run individually)

Standalone scripts include:

```tcl
set samples [readwav original.wav]
set tags [tag]
```

…and end with:

```tcl
puts "[outtag] [triangle]"
writewav output.wav
```

Standalone scripts:
- Load the WAV file  
- Perform their effect  
- Run synthesis  
- Write the result  

**Example: `uptalk.tcl`**

```tcl
set samples [readwav original.wav]
set tags [tag]

# ... effect logic ...

puts "[outtag] [triangle]"
writewav uptalk.wav
```

---

## 2. Sub-Scripts Used Inside a Master Pipeline

When scripts are used inside `amused.tcl`, the master script handles:

- Reading input  
- Calling effects  
- Synthesis  
- Writing output  

Therefore, inside pipeline-compatible sub-scripts, you must **comment out**:

```tcl
set samples [readwav ...]
set tags [tag]
puts "[outtag] [triangle]"
writewav something.wav
```

Otherwise the effect will run twice, reload the WAV mid-pipeline, or overwrite processed DSP arrays.

Example from the master script:

```tcl
if 1 {
    source long.tcl
    set samples [commit]
    set tags [tag]
}
```

The master script performs final synthesis:

```tcl
puts "[outtag] [triangle]"
writewav final_output.wav
```

---

# Synthesis Modes: `triangle` vs `trimangle`

Based on the C engine:

## `triangle`
- Applies pitch modifications  
- Applies RMS modifications  
- Computes output tags  
- **Does NOT stretch time**  

Used for:
- Uptalk  
- Pitch boost  
- RMS expressive changes  
- Amused speech contours  

## `trimangle`
- Everything `triangle` does  
- **PLUS time-domain stretching**  

Used for:
- Vowel lengthening  
- Silence stretching  
- Slow speech effects  

**Quick rule:**

| Effect Type | Use |
|-------------|-----|
| Pitch / RMS change | `triangle` |
| Time stretching | `trimangle` |

---

# Example: Master Amused Speech Pipeline

```tcl
set samples [readwav original.wav]
set tags [tag]

# Set expressive parameters
set ::minVowelStretch 3.0
set ::maxVowelStretch 4.5
set ::silenceMult 0.85
set ::minPitch 1.20
set ::maxPitch 1.40
set ::rmsCeiling 1.35
set ::rmsFloor   700
set ::minRMS 1.3
set ::maxRMS 1.6

# Stage 1: Lengthening
if 1 {
    source long.tcl
    puts "Lengthening applied"
    set samples [commit]
    set tags [tag]
}

# Final synthesis
puts "[outtag] [triangle]"
writewav longVowel.wav
```

---

# Example: Standalone Uptalk Script

```tcl
set samples [readwav original.wav]
set tags [tag]

proc longvowels samples {
    # identify vowel segments by pitch
}

foreach {s e} [longvowels $samples] {
    for {set i $s} {$i<$e} {incr i 882} {
        set frac [expr ($i-$s)*1.0/($e-$s)]
        set mult [expr 0.9 + 1.8 * $frac]
        setpitch $i $mult
    }
}

puts "[outtag] [triangle]"
writewav uptalk.wav
```

---

# License

This project is intended for academic and research use.

---

# Author
Nathan Azagury  
Electrical & Computer Engineering  
Binghamton University
