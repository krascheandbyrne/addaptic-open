#!/usr/bin/env bash
#
#	warnings.sh BUILDLOG
#
#	Lists every warning in a build's log once -- GCC's and clang's
#	"file:line:col: warning:", MSVC's "file(line,col): warning Cnnnn:" --
#	and fails on any at all. The engine and the test programs have none
#	under GCC, clang or MSVC's /W4 /permissive-, so one appearing is news;
#	and on Windows the warnings are where 32-bit faults show first --
#	C4293 and C4244 named every one CI's first run found, all of them in
#	test programs. The count names the engine's own files apart, LZMath,
#	LZNumber and LZAddaptic.h, since one there is the more serious.
#
#	A log must be from a clean build -- an incremental one reports only
#	what it recompiled, and reads as clean. Every CI build is clean.
#

log="$1"

if [ -z "$log" ] || [ ! -f "$log" ]; then
	echo "usage: ci/warnings.sh BUILDLOG" >&2
	exit 64
fi

# Each warning once: a source compiled for every target warns once per
#	target, and MSVC adds the project's name at the end of the line...
all=$(grep -E ': warning: |\): warning C[0-9]+' "$log" | tr -d '\r' | sed -E 's/ \[[^]]*\.vcxproj\]$//' | sort -u)
engine=$(printf '%s\n' "$all" | grep -E '(^|[/\\])(LZMath|LZNumber)\.(cpp|h)[:(]|(^|[/\\])LZAddaptic\.h[:(]')

n_all=$(printf '%s\n' "$all" | grep -c . || true)
n_engine=$(printf '%s\n' "$engine" | grep -c . || true)

echo "$n_all warning(s) in all, $n_engine in the engine's own files"

if [ "$n_all" -gt 0 ]; then
	printf '%s\n' "$all" | head -60
fi

if [ "$n_engine" -gt 0 ]; then
	echo
	echo "THE ENGINE HAS A WARNING:"
	printf '%s\n' "$engine"
fi

if [ "$n_all" -gt 0 ]; then
	echo
	echo "A BUILD WITH NO WARNINGS HAS ONE NOW -- fix it, or suppress it where it arises with a note"
	echo "saying why, as AddapticTests.cpp does for GCC's -Wmismatched-new-delete"
	exit 1
fi

exit 0
