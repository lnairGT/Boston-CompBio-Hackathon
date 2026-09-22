"""Disk-cached, rate-limited HTTP with retries.

Caching is what makes the pipeline reproducible and polite: a stage re-run
costs no network traffic, and a reviewer can replay a whole run offline from
the cache directory.
"""

from __future__ import annotations

import hashlib
import json
import logging
import time
from pathlib import Path
from typing import Any

import requests

from .config import CACHE_DIR, USER_AGENT

log = logging.getLogger("ind2b.http")

# Per-host minimum seconds between requests.
_RATE_LIMITS: dict[str, float] = {
    "api.platform.opentargets.org": 0.20,
    "search.rcsb.org": 0.20,
    "data.rcsb.org": 0.20,
    "files.rcsb.org": 0.10,
    "rest.uniprot.org": 0.20,
    "www.ebi.ac.uk": 0.34,
    "clinicaltrials.gov": 0.34,
}
_DEFAULT_RATE = 0.25
_last_call: dict[str, float] = {}

_session: requests.Session | None = None

RETRY_STATUS = {429, 500, 502, 503, 504}


def session() -> requests.Session:
    global _session
    if _session is None:
        s = requests.Session()
        s.headers.update({"User-Agent": USER_AGENT, "Accept": "application/json"})
        _session = s
    return _session


def _host_of(url: str) -> str:
    return url.split("://", 1)[-1].split("/", 1)[0]


def _throttle(url: str) -> None:
    host = _host_of(url)
    gap = _RATE_LIMITS.get(host, _DEFAULT_RATE)
    now = time.monotonic()
    prev = _last_call.get(host)
    if prev is not None:
        wait = gap - (now - prev)
        if wait > 0:
            time.sleep(wait)
    _last_call[host] = time.monotonic()


def _cache_key(method: str, url: str, params: Any, body: Any) -> str:
    blob = json.dumps(
        {"m": method.upper(), "u": url, "p": params, "b": body},
        sort_keys=True,
        default=str,
    )
    return hashlib.sha256(blob.encode()).hexdigest()[:40]


def _cache_path(key: str, suffix: str = ".json") -> Path:
    return CACHE_DIR / key[:2] / f"{key}{suffix}"


def request_json(
    url: str,
    *,
    method: str = "GET",
    params: dict[str, Any] | None = None,
    json_body: dict[str, Any] | None = None,
    use_cache: bool = True,
    max_attempts: int = 4,
    timeout: int = 60,
) -> Any:
    """Fetch JSON with disk cache, throttling and bounded retries.

    Raises ``requests.HTTPError`` on a non-retryable status, and
    ``RuntimeError`` when every attempt at a retryable status is exhausted.
    """
    key = _cache_key(method, url, params, json_body)
    path = _cache_path(key)
    if use_cache and path.exists():
        try:
            return json.loads(path.read_text())
        except json.JSONDecodeError:
            log.warning("corrupt cache entry %s, refetching", path)

    last_status: int | None = None
    last_error: Exception | None = None
    for attempt in range(1, max_attempts + 1):
        _throttle(url)
        try:
            resp = session().request(
                method, url, params=params, json=json_body, timeout=timeout
            )
        except (
            requests.ConnectionError,   # includes ProxyError
            requests.Timeout,
            requests.exceptions.ChunkedEncodingError,
        ) as exc:
            # Transport-level failures are transient far more often than not:
            # a dropped proxy connection mid-run should cost a retry, not the
            # whole stage. Status-code retries alone do not cover these.
            last_error = exc
            sleep_s = min(2 ** attempt, 30)
            log.warning(
                "%s -> %s (attempt %d/%d), sleeping %.1fs",
                url, type(exc).__name__, attempt, max_attempts, sleep_s,
            )
            time.sleep(sleep_s)
            continue

        if resp.status_code in RETRY_STATUS:
            last_status = resp.status_code
            sleep_s = min(2 ** attempt, 30)
            retry_after = resp.headers.get("Retry-After")
            if retry_after:
                try:
                    sleep_s = max(sleep_s, min(float(retry_after), 60))
                except ValueError:
                    pass
            log.warning(
                "%s -> %s (attempt %d/%d), sleeping %.1fs",
                url, resp.status_code, attempt, max_attempts, sleep_s,
            )
            time.sleep(sleep_s)
            continue
        resp.raise_for_status()
        data = resp.json()
        if use_cache:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(data))
        return data

    if last_error is not None:
        raise RuntimeError(
            f"{url} unreachable after {max_attempts} attempts: "
            f"{type(last_error).__name__}: {last_error}"
        ) from last_error
    raise RuntimeError(
        f"{url} kept returning {last_status} after {max_attempts} attempts"
    )


def download(
    url: str,
    dest: Path,
    *,
    use_cache: bool = True,
    timeout: int = 120,
    max_attempts: int = 4,
) -> Path:
    """Download a file to ``dest``, skipping an existing non-empty copy.

    Writes to a ``.part`` file and renames on completion, so an interrupted
    download never leaves a truncated structure that a later run would treat
    as cached.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    if use_cache and dest.exists() and dest.stat().st_size > 0:
        return dest

    last_error: Exception | None = None
    for attempt in range(1, max_attempts + 1):
        _throttle(url)
        tmp = dest.with_suffix(dest.suffix + ".part")
        try:
            with session().get(url, stream=True, timeout=timeout) as resp:
                resp.raise_for_status()
                with tmp.open("wb") as fh:
                    for chunk in resp.iter_content(chunk_size=1 << 16):
                        if chunk:
                            fh.write(chunk)
            tmp.replace(dest)
            return dest
        except (
            requests.ConnectionError,
            requests.Timeout,
            requests.exceptions.ChunkedEncodingError,
        ) as exc:
            last_error = exc
            tmp.unlink(missing_ok=True)
            sleep_s = min(2 ** attempt, 30)
            log.warning(
                "download %s -> %s (attempt %d/%d), sleeping %.1fs",
                url, type(exc).__name__, attempt, max_attempts, sleep_s,
            )
            time.sleep(sleep_s)
        except requests.HTTPError:
            tmp.unlink(missing_ok=True)
            raise

    raise RuntimeError(
        f"{url} download failed after {max_attempts} attempts: "
        f"{type(last_error).__name__}: {last_error}"
    ) from last_error
