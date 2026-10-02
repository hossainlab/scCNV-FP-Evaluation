
import os, time, requests

def fetch(url, dst, expected=None, tries=8, chunk=1<<22):
    """Resumable streaming download with Range retries. Returns final size."""
    for attempt in range(tries):
        have = os.path.getsize(dst) if os.path.exists(dst) else 0
        if expected and have == expected:
            return have
        headers = {"Range": f"bytes={have}-"} if have else {}
        try:
            with requests.get(url, headers=headers, stream=True, timeout=(30, 180)) as r:
                if r.status_code not in (200, 206):
                    raise IOError(f"HTTP {r.status_code}")
                total = expected
                if total is None:
                    cr = r.headers.get("Content-Range")
                    total = int(cr.split("/")[-1]) if cr else (
                        int(r.headers["Content-Length"]) + have if "Content-Length" in r.headers else None)
                mode = "ab" if have and r.status_code == 206 else "wb"
                if mode == "wb":
                    have = 0
                with open(dst, mode) as fh:
                    for blk in r.iter_content(chunk_size=chunk):
                        if blk:
                            fh.write(blk); have += len(blk)
            if total is None or have >= total:
                return have
            expected = total
        except Exception as e:
            time.sleep(min(2 ** attempt, 20))
    raise IOError(f"download failed after {tries} tries: {dst} ({have} bytes)")
