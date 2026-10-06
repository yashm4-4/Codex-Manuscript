#!/usr/bin/env python3
"""Bounded public metadata/header acquisition only; no association analysis."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import sys
import urllib.error
import urllib.request
import zlib

ROOT = Path(__file__).resolve().parents[1] / 'audits/alternatives/downloads'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--url', required=True)
    p.add_argument('--label', required=True)
    p.add_argument('--max-bytes', type=int, default=262144)
    p.add_argument('--range', action='store_true')
    p.add_argument('--range-start', type=int, default=0)
    args = p.parse_args()
    if not args.label.replace('_', '').replace('-', '').isalnum():
        raise ValueError('Unsafe output label')
    ROOT.mkdir(parents=True, exist_ok=True)
    body_path = ROOT / (args.label + '.bin')
    receipt_path = ROOT / (args.label + '.request.json')
    if body_path.exists() or receipt_path.exists():
        raise RuntimeError('Refuse overwrite; use a new attempt label')
    headers = {'User-Agent': 'COPD-scientific-preflight/1.0', 'Accept-Encoding': 'identity'}
    if args.range:
        if args.range_start < 0:
            raise ValueError('Range start must be nonnegative')
        headers['Range'] = f'bytes={args.range_start}-{args.range_start + args.max_bytes - 1}'
    receipt = {'requested_url': args.url, 'accessed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'request_headers': headers, 'maximum_bytes_requested_or_retained': args.max_bytes,
               'purpose': 'public source metadata or bounded schema sample; no fine-mapping', 'command': sys.argv}
    body = b''
    try:
        request = urllib.request.Request(args.url, headers=headers)
        try:
            response = urllib.request.urlopen(request, timeout=40)
        except urllib.error.HTTPError as exc:
            response = exc
        with response:
            receipt.update(http_status=response.code, response_url=response.geturl(),
                           response_headers=dict(response.headers.items()))
            body = response.read(args.max_bytes)
            length = response.headers.get('Content-Length')
            receipt['stored_body_reaches_declared_response_length'] = length is not None and len(body) == int(length)
            receipt['response_EOF_before_byte_limit'] = len(body) < args.max_bytes
            receipt['retained_prefix_only'] = response.code == 206 or not (receipt['stored_body_reaches_declared_response_length'] or receipt['response_EOF_before_byte_limit'])
    except Exception as exc:
        receipt.update(error=repr(exc), http_status=None, retained_prefix_only=True)
    with body_path.open('xb') as stream:
        stream.write(body)
    receipt.update(stored_path=str(body_path.relative_to(ROOT.parents[2])), stored_bytes=len(body),
                   stored_sha256=hashlib.sha256(body).hexdigest())
    if body.startswith(b'\x1f\x8b'):
        remaining, decoded = body, b''
        try:
            while remaining and len(decoded) < 1048576:
                inflater = zlib.decompressobj(31)
                chunk = inflater.decompress(remaining, 1048576 - len(decoded))
                decoded += chunk
                remaining = inflater.unused_data
                if not remaining:
                    break
            lines = decoded.decode('utf-8', errors='replace').splitlines()
            receipt['decompressed_prefix_bytes'] = len(decoded)
            receipt['header_and_first_records'] = lines[:6]
            receipt['prefix_decompression_does_not_verify_whole_file_CRC'] = True
        except Exception as exc:
            receipt['prefix_decompression_error'] = repr(exc)
    with receipt_path.open('x') as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write('\n')
    print(json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()
