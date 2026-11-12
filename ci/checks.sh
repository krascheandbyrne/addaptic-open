#!/usr/bin/env bash
#
#	checks.sh BINDIR [--sanitized]
#
#	The checks every CI job runs, against the binaries in BINDIR, from
#	the repository root. One script for every platform, so that the jobs
#	cannot drift apart -- two copies of a check are two things to keep in
#	step. It runs under Linux bash and the Git Bash of a Windows runner.
#
#	  accuracy     AddapticTests, the suite
#	  magnitudes   AddapticSweep, against AddapticSweepExpect.h
#	  threads      AddapticThreads
#	  audit        auditrounding.py: every sweep case against its correctly
#	               rounded value, worked in Python's decimal
#	  calibration  AddapticSweep --calibrate, which must reproduce the
#	               committed AddapticSweepExpect.h exactly: an engine that
#	               answers better than its expectations record fails until
#	               the header is regenerated, so they never fall behind
#
#	--sanitized runs the first three and the audit: its driver takes every
#	sweep case through the sanitizer, and a report fails it, since
#	auditrounding.py reads the driver's exit status. The calibration asks
#	only about answers, which a sanitizer build does not change.
#
#	Every check runs whatever the one before did, each writes its whole
#	output into BINDIR, and the exit status is the number that failed.
#

bin="$1"
sanitized=0
[ "$2" = "--sanitized" ] && sanitized=1

if [ -z "$bin" ] || [ ! -d "$bin" ]; then
	echo "usage: ci/checks.sh BINDIR [--sanitized]" >&2
	exit 64
fi

# A Windows build's binaries carry .exe; Git Bash runs them either way,
# but the test for their presence does not...
x=""
[ -f "$bin/AddapticTests.exe" ] && x=".exe"

# The first Python that runs and has decimal: on a Windows runner python3
#	can be a stub that exists and does not...
py=""
for c in python3 python; do
	if command -v "$c" > /dev/null && "$c" -c 'import decimal' > /dev/null 2>&1; then
		py="$c"
		break
	fi
done

failed=0

check() {
	local name="$1"; shift
	echo
	echo "=== $name"
	if "$@"; then
		echo "--- $name: passed"
	else
		echo "--- $name: FAILED"
		failed=$((failed + 1))
	fi
}

# A run that died without naming a check -- a sanitizer's report, a crash
#	-- is shown in its own words: the report from its first line, or else the
#	run's last lines...
explain() {
	local out="$1"
	local mark='ERROR: AddressSanitizer\|ERROR: LeakSanitizer\|runtime error:'

	if grep -q "$mark" "$out"; then
		grep -m1 -A40 "$mark" "$out"
	elif ! grep -q '^  FAIL' "$out"; then
		tail -15 "$out"
	fi
}

accuracy() {
	"$bin/AddapticTests$x" > "$bin/accuracy.txt" 2>&1
	local rc=$?
	[ $rc -ne 0 ] && explain "$bin/accuracy.txt"
	grep '^  FAIL' "$bin/accuracy.txt" | head -40
	tail -1 "$bin/accuracy.txt"
	return $rc
}

magnitudes() {
	"$bin/AddapticSweep$x" > "$bin/magnitudes.txt" 2>&1
	local rc=$?
	[ $rc -ne 0 ] && explain "$bin/magnitudes.txt"
	grep '^  FAIL\|^  UNKNOWN' "$bin/magnitudes.txt" | head -40
	grep 'held' "$bin/magnitudes.txt"
	return $rc
}

threads() {
	"$bin/AddapticThreads$x" > "$bin/threads.txt" 2>&1
	local rc=$?
	[ $rc -ne 0 ] && explain "$bin/threads.txt"
	tail -3 "$bin/threads.txt"
	return $rc
}

audit() {
	if [ -z "$py" ]; then
		echo "no Python with decimal on this machine"
		return 1
	fi
	"$py" auditrounding.py "$bin/auditrounding$x" > "$bin/audit.txt" 2>&1
	local rc=$?
	[ $rc -ne 0 ] && explain "$bin/audit.txt"
	head -1 "$bin/audit.txt"
	sed -n '/UNEXPLAINED/,$p' "$bin/audit.txt" | head -45
	return $rc
}

calibration() {
	"$bin/AddapticSweep$x" --calibrate 2>/dev/null | tr -d '\r' > "$bin/calibration.h"
	if tr -d '\r' < AddapticSweepExpect.h | cmp -s - "$bin/calibration.h"; then
		echo "identical to the committed AddapticSweepExpect.h ($(wc -l < "$bin/calibration.h" | tr -d " ") lines)"
		return 0
	fi
	echo "AddapticSweepExpect.h no longer matches the engine; regenerate it with AddapticSweep --calibrate"
	tr -d '\r' < AddapticSweepExpect.h | diff - "$bin/calibration.h" | head -40
	return 1
}

check accuracy accuracy
check magnitudes magnitudes
check threads threads

check audit audit

if [ $sanitized -eq 0 ]; then
	check calibration calibration
fi

echo
if [ $failed -eq 0 ]; then
	echo "ALL CHECKS PASSED"
else
	echo "$failed CHECK(S) FAILED"
fi

exit $failed
