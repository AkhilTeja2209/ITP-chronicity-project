#!/usr/bin/env python3
"""Print the rows of a bgzipped, tabix-indexed file (local path or http(s) URL) within a region.

Usage: tabix_query.py <file-or-url> <chrom:start-end>

Ubuntu's apt `tabix` is built without libcurl, so remote eQTL Catalogue files cannot be queried
with it; pysam wheels bundle an htslib with network support. The index is cached in the working dir.
The bundled libcurl does not know where the system CA store lives, so it is pointed there, and
plain HTTP (also served by ftp.ebi.ac.uk) is tried if HTTPS still fails.
"""
import os
import sys

for ca in ("/etc/ssl/certs/ca-certificates.crt", "/etc/pki/tls/certs/ca-bundle.crt"):
    if os.path.exists(ca):
        os.environ.setdefault("CURL_CA_BUNDLE", ca)
        break

import pysam  # noqa: E402

src, region = sys.argv[1], sys.argv[2]
chrom, rng = region.split(":")
start, end = (int(x) for x in rng.split("-"))

candidates = [src]
if src.startswith("https://"):
    candidates.append("http://" + src[len("https://"):])

last_error = None
for cand in candidates:
    try:
        tf = pysam.TabixFile(cand)
    except OSError as e:
        last_error = e
        continue
    with tf:
        for row in tf.fetch(chrom, start - 1, end):
            sys.stdout.write(row.rstrip("\r") + "\n")
    sys.exit(0)
raise SystemExit(f"could not open {src}: {last_error}")
