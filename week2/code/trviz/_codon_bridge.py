import atexit
import subprocess
from pathlib import Path
from functools import lru_cache
from threading import Lock

# Single persistent worker
_worker = None
_lock = Lock()

def _start_worker():
    global _worker
    if _worker and _worker.poll() is None:
        return _worker
    root = Path(__file__).resolve().parents[1]  # .../week2/code
    # Prefer compiled worker if present; else run Codon source
    compiled = root / "trviz_codon_worker"
    if compiled.exists():
        cmd = [str(compiled)]
    else:
        worker_src = root / "trviz_codon" / "worker.codon"
        cmd = ["codon", "run", str(worker_src)]
    _worker = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,  # line buffered
    )
    return _worker

def _stop_worker():
    global _worker
    if _worker and _worker.poll() is None:
        try:
            _worker.stdin.close()
        except Exception:
            pass
        try:
            _worker.terminate()
        except Exception:
            pass
    _worker = None

atexit.register(_stop_worker)

def codon_call_args(cmd: str, *args: object) -> str:
    """
    Send one request to the persistent Codon worker.
    Returns the payload as a single string (joined by '\n' if multiple lines).
    Raises RuntimeError on worker-side error.
    """
    with _lock:
        p = _start_worker()
        if not p or p.stdin is None or p.stdout is None:
            raise RuntimeError("Cannot start Codon worker")

        # Request: header "CMD N", then N arg lines
        header = f"{cmd} {len(args)}\n"
        try:
            p.stdin.write(header)
            for a in args:
                p.stdin.write(f"{a}\n")
            p.stdin.flush()
        except Exception as e:
            raise RuntimeError(f"Failed to write to Codon worker: {e}")

        # Response: "OK" or "ERR: ..."; if OK, then "<L>" and then L lines
        status = p.stdout.readline().rstrip("\n")
        if not status:
            err = (p.stderr.read() or "").strip()
            raise RuntimeError(f"Codon worker died. stderr: {err}")

        if status.startswith("ERR:"):
            raise RuntimeError(status)

        line_count_s = p.stdout.readline().rstrip("\n")
        if not line_count_s:
            raise RuntimeError("Codon worker protocol error: missing line count")
        try:
            line_count = int(line_count_s)
        except Exception:
            raise RuntimeError(f"Codon worker protocol error: bad line count '{line_count_s}'")

        lines = []
        for _ in range(line_count):
            lines.append(p.stdout.readline().rstrip("\n"))
        return "\n".join(lines)

# Optional tiny cache for expensive, repeated commands:
@lru_cache(maxsize=128)
def codon_call_args_cached(cmd: str, args_tuple: tuple[str, ...]) -> str:
    return codon_call_args(cmd, *args_tuple)
