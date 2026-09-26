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
    for {set i 0} {$i < $::env(RUSTY_DRAW_DATA_COUNT)} {incr i} {
        set path [file join $::env(RUSTY_DRAW_DATA_$i) $name]
        if {[file isfile $path]} {return $path}
    }
    lappend ::missing $name
    error "fixture file could not be found: $name"
}
proc runCommand {command args} {
    lappend ::commands $command
    if {$::backend eq "rust"} {
        # DRAW permits numeric Tcl expressions, e.g. 2e-7+1e-14. Preserve that
        # behavior through the actual interpreter, not a Python expression parser.
        # Names precede the numbers; prism's numbers are followed by a mode.
        set leading {box {1 end} trotate {1 end} ttranslate {1 end} polyline {1 end} prism {2 4} pcylinder {1 end}}
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
    }
    if {$result ne ""} {logPuts $result}
    if {$command in {checkshape nbshapes vprops sprops lprops isbbinterf isdeleted}} {incr ::queries}
    # DBRep::Set binds DRAW shape names as Tcl variables as well.
    if {$::backend eq "rust"} {
        if {$command in {box pcylinder polyline mkplane prism generated modified}} {
            interp eval testcase [list set [lindex $args 0] [lindex $args 0]]
        }
        if {$command eq "copy"} {interp eval testcase [list set [lindex $args 1] [lindex $args 1]]}
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
    foreach command {box copy ttranslate trotate isdraw whatis checkshape nbshapes vprops sprops lprops isbbinterf explode compound bcommon bfuse restore prism polyline mkplane savehistory generated modified isdeleted pcylinder plane mkface line mkedge mkvolume bclearobjects bcleartools baddobjects baddtools bfillds bsplit bbuild} {
        interp alias testcase $command {} runCommand $command
    }
    # bugs/begin would load VISUALIZATION only when topology checks are absent.
    interp eval testcase {set Draw_Groups(TOPOLOGY\ Check\ commands) {checkshape}}
    evaluateFile [file join $::env(RUSTY_DRAW_ROOT) resources DrawResources CheckCommands.tcl]
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
# A missing capability is latched outside the test, even if its Tcl catch hides it.
if {[llength $unsupported] > 0} {set status unsupported}
if {[llength $missing] > 0} {set status missing_fixture}

set output [open $::env(RUSTY_DRAW_RESULT) w]
foreach {key value} [list status $status backend $backend version $version queries $queries \
    unsupported [join $unsupported "\n"] missing [join $missing "\n"] \
    error $caught commands [join $commands "\n"] adjudication $summary] {
    puts $output "$key [hex $value]"
}
close $output
exit 0
