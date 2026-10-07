"""Audited anonymous public acquisition; does not interpret scientific payloads."""
import datetime
import hashlib
import json
from pathlib import Path
import time
import requests


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def fetch(url, path, timeout=(30, 180)):
    path = Path(path)
    receipt_path = Path(str(path) + '.access.json')
    if path.exists() and receipt_path.exists():
        receipt = json.loads(receipt_path.read_text())
        if receipt.get('complete') and receipt['url'] == url:
            return receipt
        raise RuntimeError('Existing incomplete acquisition must be preserved: ' + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = Path(str(path) + '.partial')
    if partial.exists():
        raise RuntimeError('Preserve/reconcile previous partial before retry: ' + str(partial))
    rec = {'url': url, 'started_utc': utc(), 'request_headers': {'Accept-Encoding': 'identity'},
           'path': str(path), 'complete': False}
    sha = hashlib.sha256(); md5 = hashlib.md5(); size = 0
    start = time.monotonic()
    try:
        with requests.get(url, stream=True, headers={'Accept-Encoding': 'identity'}, timeout=timeout) as response:
            rec.update(http_status=response.status_code, final_url=response.url,
                       response_headers=dict(response.headers))
            with partial.open('wb') as out:
                for block in response.iter_content(1024 * 1024):
                    if block:
                        out.write(block); sha.update(block); md5.update(block); size += len(block)
            rec['content_length_matches'] = (not response.headers.get('Content-Length') or
                                              int(response.headers['Content-Length']) == size)
            rec['complete'] = response.status_code == 200 and rec['content_length_matches']
            partial.rename(path)
    except Exception as exc:
        rec['error'] = repr(exc)
    rec.update(completed_utc=utc(), bytes=size, sha256=sha.hexdigest(), md5=md5.hexdigest(),
               seconds=time.monotonic() - start)
    receipt_path.write_text(json.dumps(rec, indent=2) + '\n')
    if not rec['complete']:
        raise RuntimeError('Acquisition did not complete: ' + str(receipt_path))
    return rec


if __name__ == '__main__':
    import sys
    print(json.dumps(fetch(sys.argv[1], sys.argv[2]), indent=2))
