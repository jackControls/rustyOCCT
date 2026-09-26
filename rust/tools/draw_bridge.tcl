# Test-only, headless DRAW compatibility host; works with Tcl 8.5 and newer.
# The original test files and assertion library are evaluated without rewriting.
# Run through run_upstream_tests.py, which enforces a process-group timeout.
set backend $::env(RUSTY_DRAW_BACKEND)
set unsupported {}
set missing {}
set messages ""
set queries 0
set commands {}
set caught ""

proc hex {value} {
    binary scan [encoding convertto utf-8 $value] H* encoded
    return $encoded
}
proc unhex {value} {
    return [encoding convertfrom utf-8 [binary format H* $value]]
}
proc unsupportedCommand {args} {
    lappend ::unsupported [join $args " "]
    error "unsupported command or signature: [join $args { }]"
}
proc logPuts {args} {
    set newline 1
    if {[lindex $args 0] eq "-nonewline"} {
        set newline 0
        set args [lrange $args 1 end]
    }
    if {[llength $args] == 2 && [lindex $args 0] in {stdout stderr}} {
        set args [lrange $args 1 end]
    }
    if {[llength $args] != 1} {return [unsupportedCommand puts {*}$args]}
    set text [lindex $args 0]
    append ::messages $text
    if {$newline} {append ::messages "\n"; puts $text} else {puts -nonewline $text}
}
proc metadata {args} {return ""}
proc cpuLimit {seconds} {
    # The outer process supervisor enforces a stricter wall-time limit.
    if {![string is integer -strict $seconds] || $seconds <= 0} {
        return [unsupportedCommand cpulimit $seconds]
    }
}
proc locateData {name} {
    # Data files must already be present. Never download or invent a fixture.
    # As upstream's locate_data_file (TestCommands.tcl): the case's own data
    # folder first, then each data directory and all of its subdirectories
    # breadth-first, skipping those whose names start with a dot (sorted
    # here, so the choice between equal names is deterministic).
    set own [file join $::env(RUSTY_DRAW_CASE_DIR) data $name]
    if {[file isfile $own]} {return [file normalize $own]}
    for {set i 0} {$i < $::env(RUSTY_DRAW_DATA_COUNT)} {incr i} {
        set queue [list $::env(RUSTY_DRAW_DATA_$i)]
        while {[llength $queue] > 0} {
            set dir [lindex $queue 0]
            set queue [lrange $queue 1 end]
            if {[string match .* [file tail $dir]]} {continue}
            set path [file join $dir $name]
            if {[file isfile $path]} {return [file normalize $path]}
            lappend queue {*}[lsort [glob -nocomplain -directory $dir -type d *]]
        }
    }
    lappend ::missing $name
    error "fixture file could not be found: $name"
}
# Viewer commands (REVIEW_NOTES.md R9) are recorded and not run, on both
# backends: the case then reports viewer_skipped instead of pass.
set viewer {}
proc viewerCommand {command args} {
    lappend ::viewer [join [linsert $args 0 $command] " "]
    return ""
}
proc ploadCommand {args} {
    # Modelling is always loaded; the viewer's plugin is never needed, since
    # viewer commands are recorded. Any other module is a capability.
    foreach module $args {
        if {$module ni {MODELING VISUALIZATION TOPTEST}} {
            if {$::backend eq "rust"} {return [unsupportedCommand pload {*}$args]}
            return [uplevel #0 [linsert $args 0 pload]]
        }
    }
    return ""
}
# Native selector (run_upstream_tests.py): record every pick of a native
# explode as kind, measure and centre of gravity, so the Rust adapter can
# select its own entity by geometry instead of OCCT's exploration order.
set explodes 0
proc recordPicks {ordinal arguments names} {
    set kind other
    # DBRep's explode reads the type from its first letter.
    switch -nocase -glob -- [lindex $arguments 1] {
        f* {set kind face}
        e* {set kind edge}
    }
    set stream [open $::env(RUSTY_DRAW_SELECTOR_OUT) a]
    set index 0
    foreach name $names {
        incr index
        if {$kind eq "other"} {
            puts $stream "pick $ordinal $index other"
            continue
        }
        set measure [expr {$kind eq "face" ? "sprops" : "lprops"}]
        set text [uplevel #0 [list $measure $name _rusty_gx _rusty_gy _rusty_gz]]
        if {![regexp {Mass +: +([-0-9.+eE]+)} $text all mass]} {error "selector: no mass for $name"}
        set centre {}
        foreach axis {_rusty_gx _rusty_gy _rusty_gz} {lappend centre [uplevel #0 [list dval $axis]]}
        puts $stream "pick $ordinal $index $kind $mass [join $centre { }]"
    }
    close $stream
}
proc runCommand {command args} {
    lappend ::commands $command
    if {$command eq "explode"} {incr ::explodes}
    if {$::backend eq "rust"} {
        # DRAW permits numeric Tcl expressions, e.g. 2e-7+1e-14. Preserve that
        # behavior through the actual interpreter, not a Python expression parser.
        # Names precede the numbers; prism's numbers are followed by a mode.
        set leading {box {1 end} trotate {1 end} ttranslate {1 end} polyline {1 end} prism {2 4} pcylinder {1 end} pcone {1 end}}
        if {[dict exists $leading $command]} {
            lassign [dict get $leading $command] first last
            set converted [lrange $args 0 [expr {$first - 1}]]
            foreach value [lrange $args $first $last] {
                if {[catch {interp eval testcase [list expr $value]} number]} {
                    return [unsupportedCommand $command {*}$args]
                }
                lappend converted $number
            }
            if {$last ne "end"} {lappend converted {*}[lrange $args [expr {$last + 1}] end]}
            set args $converted
        }
        set tokens {}
        foreach token [linsert $args 0 $command] {lappend tokens [hex $token]}
        puts $::worker [join $tokens " "]
        flush $::worker
        if {[gets $::worker response] < 0} {error "Rust worker closed its output"}
        if {![regexp {^(OK|ERROR|UNSUPPORTED) ([0-9a-f]*)$} $response all status encoded]} {
            error "invalid Rust worker response"
        }
        set result [unhex $encoded]
        if {$status eq "UNSUPPORTED"} {
            lappend ::unsupported $result
            error $result
        }
        if {$status ne "OK"} {error $result}
    } else {
        set result [uplevel #0 [linsert $args 0 $command]]
        if {$command eq "explode" && [info exists ::env(RUSTY_DRAW_SELECTOR_OUT)]} {
            recordPicks $::explodes $args $result
        }
    }
    if {$result ne ""} {logPuts $result}
    if {$command in {checkshape nbshapes vprops sprops lprops isbbinterf isdeleted}} {incr ::queries}
    # DBRep::Set binds DRAW shape names as Tcl variables as well.
    if {$::backend eq "rust"} {
        if {$command in {box pcylinder pcone polyline mkplane prism generated modified}} {
            interp eval testcase [list set [lindex $args 0] [lindex $args 0]]
        }
        if {$command in {copy restore}} {interp eval testcase [list set [lindex $args 1] [lindex $args 1]]}
        if {$command eq "explode"} {
            foreach name $result {interp eval testcase [list set $name $name]}
        }
    }
    return $result
}
proc evaluateFile {path} {
    set stream [open $path r]
    fconfigure $stream -encoding utf-8
    set script [read $stream]
    close $stream
    interp eval testcase $script
}

set version ""
if {[catch {
    if {$backend eq "occt"} {
        pload MODELING
        set version [dversion]
    } else {
        set worker [open [list | $::env(RUSTY_DRAW_WORKER) 2>@ stderr] r+]
        fconfigure $worker -buffering line -encoding utf-8
        set version "rusty-occt DRAW adapter"
    }
    interp create -safe testcase
    # Safe Tcl supplies real loops, substitutions, procedures and assertions.
    # Filesystem/process/viewer commands are unavailable to the test interpreter.
    interp alias testcase unknown {} unsupportedCommand
    interp alias testcase puts {} logPuts
    interp alias testcase help {} metadata
    interp alias testcase cpulimit {} cpuLimit
    interp alias testcase locate_data_file {} locateData
    foreach command {box copy ttranslate trotate isdraw whatis checkshape nbshapes vprops sprops lprops isbbinterf explode compound bcommon bfuse restore prism polyline mkplane savehistory generated modified isdeleted pcylinder pcone plane mkface line mkedge mkvolume bclearobjects bcleartools baddobjects baddtools bfillds bsplit bbuild} {
        interp alias testcase $command {} runCommand $command
    }
    # The variables upstream's _run_test (TestCommands.tcl) sets for a case.
    # The image directory is this case's output directory; viewer commands
    # are recorded, so nothing is written there.
    foreach {variable value} [list casename $::env(RUSTY_DRAW_CASE) \
            groupname $::env(RUSTY_DRAW_GROUP) gridname $::env(RUSTY_DRAW_GRID) \
            dirname [file join $::env(RUSTY_DRAW_ROOT) tests] test_image $::env(RUSTY_DRAW_CASE) \
            imagedir [file dirname $::env(RUSTY_DRAW_RESULT)]] {
        interp eval testcase [list set $variable $value]
    }
    # bugs/begin would load VISUALIZATION only when topology checks are absent.
    interp eval testcase {set Draw_Groups(TOPOLOGY\ Check\ commands) {checkshape}}
    interp alias testcase pload {} ploadCommand
    evaluateFile [file join $::env(RUSTY_DRAW_ROOT) resources DrawResources CheckCommands.tcl]
    # After the check library: checkview and checkcolor are its procedures.
    set stream [open [file join $::env(RUSTY_DRAW_ROOT) rust fixtures draw-viewer-commands.txt] r]
    foreach line [split [read $stream] "\n"] {
        if {$line eq "" || [string index $line 0] eq "#"} {continue}
        interp alias testcase $line {} viewerCommand $line
    }
    close $stream
    for {set i 0} {$i < $::env(RUSTY_DRAW_SOURCE_COUNT)} {incr i} {
        evaluateFile $::env(RUSTY_DRAW_SOURCE_$i)
    }
} caught options]} {
    logPuts "Tcl Exception: $caught"
} else {
    set caught ""
}
if {[info exists worker] && [catch {close $worker} closeError]} {
    set caught "Rust worker failed: $closeError"
    logPuts $caught
}

# Use the original log adjudicator and group/grid parse.rules too. This handles
# print-only assertion failures, REQUIRED/TODO rules and missing completion marks.
proc help {args} {return ""}
source [file join $::env(RUSTY_DRAW_ROOT) resources DrawResources CheckCommands.tcl]
source [file join $::env(RUSTY_DRAW_ROOT) resources DrawResources TestCommands.tcl]
_check_log [file join $::env(RUSTY_DRAW_ROOT) tests] \
    $::env(RUSTY_DRAW_GROUP) $::env(RUSTY_DRAW_GRID) $::env(RUSTY_DRAW_CASE) 1 $messages summary
set status failed
set adjudication ""
if {[regexp {^CASE [^:]+: ([A-Z]+)} $summary all adjudication]} {
    switch -- $adjudication {
        OK {set status pass}
        BAD {set status known_failure}
        IMPROVEMENT {set status unexpected_improvement}
        SKIPPED {set status skipped}
    }
}
if {$caught ne ""} {set status failed}
if {$status eq "pass" && $queries == 0} {set status unverified}
# Every geometric assertion ran; the image commands were only recorded.
if {$status eq "pass" && [llength $viewer] > 0} {set status viewer_skipped}
# A missing capability is latched outside the test, even if its Tcl catch hides it.
if {[llength $unsupported] > 0} {set status unsupported}
if {[llength $missing] > 0} {set status missing_fixture}

set output [open $::env(RUSTY_DRAW_RESULT) w]
foreach {key value} [list status $status backend $backend version $version queries $queries \
    unsupported [join $unsupported "\n"] missing [join $missing "\n"] viewer [join $viewer "\n"] \
    error $caught commands [join $commands "\n"] adjudication $summary] {
    puts $output "$key [hex $value]"
}
close $output
exit 0
