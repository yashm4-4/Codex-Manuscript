#!/usr/bin/env python3
"""Bounded public-resource retrieval with exact acquired-byte provenance.

This is access/schema inspection only. It neither downloads complete large GWAS
files nor computes association statistics, LD, fine mapping or model scores.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.error
import urllib.request
import zlib

STAGE = Path(__file__).resolve().parents[1]


def now():
    return datetime.now(timezone.utc).isoformat()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--url', required=True)
    parser.add_argument('--label', required=True)
    parser.add_argument('--limit', type=int, default=1048576)
    parser.add_argument('--range', action='store_true')
    parser.add_argument('--gzip-prefix', action='store_true')
    parser.add_argument('--timeout', type=int, default=35)
    args = parser.parse_args()
    assert args.label.replace('_', '').replace('-', '').isalnum()
    assert 1 <= args.limit <= 4 * 1024 * 1024
    output = STAGE / 'audits/direct/acquired' / args.label
    output.mkdir(parents=True, exist_ok=False)
    headers = {'User-Agent': 'COPD-scientific-preflight/1.0 (bounded public metadata access)',
               'Accept-Encoding': 'identity'}
    if args.range:
        headers['Range'] = f'bytes=0-{args.limit - 1}'
    record = {'url': args.url, 'requested_utc': now(), 'request_headers': headers,
              'read_limit_bytes': args.limit, 'purpose': 'Public metadata/header/schema access audit; no statistical analysis'}
    request = urllib.request.Request(args.url, headers=headers)
    body = b''
    try:
        with urllib.request.urlopen(request, timeout=args.timeout) as response:
            body = response.read(args.limit + 1)
            record.update(http_status=response.status, effective_url=response.geturl(),
                          response_headers=dict(response.headers),
                          read_truncated=len(body) > args.limit, transport_error=None)
    except urllib.error.HTTPError as error:
        body = error.read(min(args.limit, 32768) + 1)
        record.update(http_status=error.code, effective_url=error.geturl(),
                      response_headers=dict(error.headers), read_truncated=len(body) > min(args.limit, 32768),
                      transport_error=str(error))
    except Exception as error:
        record.update(http_status=None, effective_url=args.url, response_headers={},
                      read_truncated=False, transport_error=f'{type(error).__name__}: {error}')
    body = body[:args.limit]
    with (output / 'response.bin').open('xb') as handle:
        handle.write(body)
    record.update(completed_utc=now(), acquired_bytes=len(body),
                  acquired_sha256=hashlib.sha256(body).hexdigest(),
                  body_path=str((output / 'response.bin').relative_to(STAGE)),
                  complete_remote_object_hash=False)
    content_range = record['response_headers'].get('Content-Range', record['response_headers'].get('content-range'))
    if content_range:
        record['content_range'] = content_range
    if args.gzip_prefix and body:
        try:
            decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
            decoded = decoder.decompress(body, 262144)
            # Save bounded decompressed bytes, not a reconstructed full file.
            with (output / 'decompressed_prefix.txt').open('xb') as handle:
                handle.write(decoded)
            record['gzip_prefix'] = {'bytes': len(decoded), 'sha256': hashlib.sha256(decoded).hexdigest(),
                                     'path': str((output / 'decompressed_prefix.txt').relative_to(STAGE)),
                                     'gzip_stream_eof_seen': decoder.eof,
                                     'decompressed_limit_bytes': 262144}
        except Exception as error:
            record['gzip_prefix_error'] = f'{type(error).__name__}: {error}'
    with (output / 'access.json').open('x') as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
