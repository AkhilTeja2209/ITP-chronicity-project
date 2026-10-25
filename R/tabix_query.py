#!/usr/bin/env python3
"""Print the rows of a bgzipped, tabix-indexed file (local path or http(s) URL) within a region.

Usage: tabix_query.py <file-or-url> <chrom:start-end>

Ubuntu's apt `tabix` is built without libcurl, so remote eQTL Catalogue files cannot be queried
with it; pysam wheels bundle an htslib with network support. The index is cached in the working dir.
The bundled libcurl does not know where the system CA store lives, so it is pointed there, and
plain HTTP (also served by ftp.ebi.ac.uk) is tried if HTTPS still fails.

Remote reads can fail part-way through a region. Rows are therefore buffered and written only
after the whole region has been read; a failed read is retried, and if every attempt fails the
script exits non-zero so the caller never mistakes a failed query for a region with no data.
"""
import os
import sys
import time

for ca in ("/etc/ssl/certs/ca-certificates.crt", "/etc/pki/tls/certs/ca-bundle.crt"):
    if os.path.exists(ca):
        os.environ.setdefault("CURL_CA_BUNDLE", ca)
        break

import pysam  # noqa: E402

ATTEMPTS = 4

src, region = sys.argv[1], sys.argv[2]
chrom, rng = region.split(":")
start, end = (int(x) for x in rng.split("-"))

candidates = [src]
if src.startswith("https://"):
    candidates.append("http://" + src[len("https://"):])


def read_region(path: str) -> list[str]:
    with pysam.TabixFile(path) as tf:
        if chrom not in tf.contigs:
            return []  # the file has no such chromosome: a genuinely empty region
        return [row.rstrip("\r") for row in tf.fetch(chrom, start - 1, end)]


last_error = None
for attempt in range(1, ATTEMPTS + 1):
    for cand in candidates:
        try:
            rows = read_region(cand)
        except (OSError, ValueError) as e:
            last_error = e
            continue
        sys.stdout.write("".join(r + "\n" for r in rows))
        sys.exit(0)
    print(f"attempt {attempt} failed: {last_error}", file=sys.stderr)
    if attempt < ATTEMPTS:
        time.sleep(5 * attempt)
sys.exit(f"could not read {region} from {src} after {ATTEMPTS} attempts: {last_error}")
