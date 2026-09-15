#!/usr/bin/env bash
# Regression suite for the A390049 sieve. Every case here is a defect that was
# live at commit a486447; see REVIEW.md for what each one was. Run via
# `make test`, or directly. Documented in AITESTING.md.
set -u
cd "$(dirname "$0")/.."

T=${T:-8}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

PASS=0; FAIL=0
ok()  { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad() { FAIL=$((FAIL+1)); printf '  FAIL  %s\n' "$1"; printf '        %s\n' "$2"; }
eq()  { if [ "$2" = "$3" ]; then ok "$1"; else bad "$1" "want [$2] got [$3]"; fi; }

terms() { ./sieve "$@" 2>/dev/null | awk '$1=="TERM"{print $2}' | sort -n | tr '\n' ' '; }
rc_of() { "$@" >"$TMP/out" 2>"$TMP/err"; echo $?; }

KNOWN="40 208 928 3904 260608 1045504 16764928 268386304 "

echo "== correctness =="

eq "eight known terms below 3e8" "$KNOWN" "$(terms 2 300000000 "$T" 1048576 "" 8192)"

REF=$(python3 tests/reference.py 3000000 | sort -n | tr '\n' ' ')
eq "agrees with tests/reference.py over [0,3e6)" "$REF" "$(terms 0 3000000 "$T" 1048576 "" 8192)"

eq "lo=0 does not hang or miss (rest=0 guard)" "$KNOWN" "$(terms 0 300000000 "$T" 1048576 "" 8192)"

echo "== window/block invariance =="
SWEEP_OK=1; SWEEP_MSG=""
EXPECT=$(terms 0 3000000 "$T" 1048576 "" 8192)
for W in 4096 65536 1048576 4194304; do
  for B in 64 1024 8192 65536 1048576; do
    G=$(terms 0 3000000 "$T" "$W" "" "$B")
    [ "$G" = "$EXPECT" ] || { SWEEP_OK=0; SWEEP_MSG="$SWEEP_MSG W=$W B=$B->[$G];"; }
  done
done
if [ $SWEEP_OK = 1 ]; then ok "20 window x block combinations agree"
else bad "20 window x block combinations agree" "$SWEEP_MSG"; fi

echo "== C1/C2 state file identity =="

# positive control first: a resume of the SAME run must still work
rm -f "$TMP/a.state"*
./sieve 16000000 20000000 "$T" 1048576 "$TMP/a.state" 8192 >/dev/null 2>&1
eq "terms file recorded"      "16764928" "$(cat "$TMP/a.state.terms" 2>/dev/null)"
R=$(./sieve 16000000 20000000 "$T" 1048576 "$TMP/a.state" 8192 2>&1 >/dev/null | grep -o 'already complete: [0-9]*/[0-9]*')
eq "same run resumes"         "already complete: 4/4" "$R"

# the defect: a bitmap from a different range used to be applied positionally,
# skipping windows that had never been sieved and still exiting "complete"
rm -f "$TMP/b.state"*
./sieve 0 3000000 "$T" 1048576 "$TMP/b.state" 8192 >/dev/null 2>&1
eq "different range refused"  "3" "$(rc_of ./sieve 16000000 20000000 "$T" 1048576 "$TMP/b.state" 8192)"
eq "  ...and says why"        "1" "$(grep -c 'different run' "$TMP/err")"
eq "different window refused" "3" "$(rc_of ./sieve 0 3000000 "$T" 65536 "$TMP/b.state" 8192)"

# a different block size is only a speed knob: it must still resume
eq "different block accepted" "0" "$(rc_of ./sieve 0 3000000 "$T" 1048576 "$TMP/b.state" 4096)"

head -c 20 "$TMP/b.state" > "$TMP/trunc.state"
eq "truncated state refused"  "3" "$(rc_of ./sieve 0 3000000 "$T" 1048576 "$TMP/trunc.state" 8192)"
head -c 200 /dev/urandom > "$TMP/junk.state"
eq "foreign state refused"    "3" "$(rc_of ./sieve 0 3000000 "$T" 1048576 "$TMP/junk.state" 8192)"

echo "== C3/C4/C5 argument validation =="
eq "hi <= lo rejected"        "2" "$(rc_of ./sieve 1000 500 "$T" 1048576 "" 8192)"
eq "hi == lo rejected"        "2" "$(rc_of ./sieve 500 500 "$T" 1048576 "" 8192)"
eq "window 0 rejected"        "2" "$(rc_of ./sieve 2 100000 "$T" 0 "" 8192)"
eq "block 0 rejected"         "2" "$(rc_of ./sieve 2 100000 "$T" 1048576 "" 0)"
eq "block > 2^40 rejected"    "2" "$(rc_of ./sieve 2 100000 "$T" 1048576 "" 2199023255552)"
eq "threads 0 rejected"       "2" "$(rc_of ./sieve 2 100000 0 1048576 "" 8192)"
eq "threads -4 rejected"      "2" "$(rc_of ./sieve 2 100000 -4 1048576 "" 8192)"
eq "oddpart threads -4 rejected" "2" "$(rc_of ./oddpart 1000 -4)"
eq "oddpart M=0 is a no-op"      "0" "$(rc_of ./oddpart 0 "$T")"

echo "== oddpart =="
OP=$(./oddpart 2000 "$T" 2>/dev/null | awk '$1=="TERM"{print $2}' | sed 's/oddpart=//' | sort -n | tr '\n' ' ')
eq "odd parts of the terms with 2^m-3 < 2000" "5 13 29 61 509 1021 " "$OP"

echo
echo "$PASS passed, $FAIL failed"
[ "$FAIL" -eq 0 ]
