#!/usr/bin/env python3
"""Report what a sieve checkpoint actually says.

The log is opened append and survives restarts, so a completion marker in it
proves nothing. The state file does: it names the range it belongs to and
carries one bit per window.

    tools/state-status.py run.state
    tools/state-status.py --legacy old-run.state   # bare pre-header bitmap

Exit 0 if every window is complete, 1 if not, 2 if the file is unreadable.
"""
import struct
import sys

HDR = "<8sIIQQQQQ"          # magic, version, hdrsize, lo, hi, window, block, nwin
MAGIC = b"A390049S"

def legacy(path):
    """A pre-header state file: a bare bitmap that names no range."""
    try:
        blob = open(path, "rb").read()
    except OSError as e:
        print("%s: %s" % (path, e.strerror), file=sys.stderr)
        return 2
    done = bin(int.from_bytes(blob, "little")).count("1")
    print("file    %s (legacy, no header)" % path)
    print("windows %d complete out of at most %d" % (done, len(blob) * 8))
    print("status  compare against the window count the run printed at startup;")
    print("        this file does not record which range it belongs to")
    return 0


def main(path):
    try:
        blob = open(path, "rb").read()
    except OSError as e:
        print("%s: %s" % (path, e.strerror), file=sys.stderr)
        return 2
    if len(blob) < struct.calcsize(HDR):
        print("%s: too short to be a state file" % path, file=sys.stderr)
        return 2
    magic, version, hdrsize, lo, hi, window, block, nwin = struct.unpack_from(HDR, blob)
    if magic != MAGIC:
        print("%s: not a sieve state file" % path, file=sys.stderr)
        return 2
    bits = blob[hdrsize:]
    if len(bits) != (nwin + 7) // 8:
        print("%s: truncated: %d of %d bitmap bytes"
              % (path, len(bits), (nwin + 7) // 8), file=sys.stderr)
        return 2
    done = bin(int.from_bytes(bits, "little")).count("1")
    print("file    %s (v%d)" % (path, version))
    print("range   [%d, %d)  window %d  block %d" % (lo, hi, window, block))
    print("windows %d of %d complete" % (done, nwin))
    if done == nwin:
        print("status  COMPLETE")
        return 0
    print("status  INCOMPLETE, %d window(s) remaining" % (nwin - done))
    return 1

if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "--legacy":
        sys.exit(legacy(args[1]) if len(args) > 1 else 2)
    sys.exit(main(args[0] if args else "run.state"))
