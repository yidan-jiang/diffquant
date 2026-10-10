"""Process-parallel map with a serial fallback (single item, jobs = 1, or platforms without process pools)."""
from __future__ import annotations


def pmap(fn, args, jobs: int = 1, chunksize: int = 4) -> list:
    args = list(args)
    if jobs <= 1 or len(args) <= 1:
        return [fn(a) for a in args]
    try:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(jobs) as ex:
            return list(ex.map(fn, args, chunksize=chunksize))
    except (PermissionError, OSError, NotImplementedError) as e:      # e.g. sandboxed / restricted environments
        print(f"[diffquant] process pool unavailable ({e.__class__.__name__}); running serially")
        return [fn(a) for a in args]
