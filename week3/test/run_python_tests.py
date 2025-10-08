import argparse, sys, time
import pytest

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--just-ms", dest="just_ms", action="store_true",
                    help="print only the runtime in milliseconds")
    ap.add_argument("--tests", default="targeted_tests.py",
                    help="path to test file(s)")
    args = ap.parse_args()

    start = time.perf_counter()
    # run pytest programmatically
    rc = pytest.main(["-q", args.tests])
    elapsed_ms = int((time.perf_counter() - start) * 1000)

    if args.just_ms:
        print(elapsed_ms)
    else:
        print("Language    Runtime")
        print("-------------------")
        print(f"python      {elapsed_ms}ms")

    # propagate pytest exit code
    sys.exit(rc)

if __name__ == "__main__":
    main()
