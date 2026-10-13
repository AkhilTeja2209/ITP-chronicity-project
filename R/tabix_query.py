#!/usr/bin/env python3
"""Print the rows of a bgzipped, tabix-indexed file (local path or https URL) within a region.

Usage: tabix_query.py <file-or-url> <chrom:start-end>

Ubuntu's apt `tabix` is built without libcurl, so remote eQTL Catalogue files cannot be queried
with it; pysam wheels bundle an htslib with HTTPS support. The index is cached in the working dir.
"""
import sys

import pysam

src, region = sys.argv[1], sys.argv[2]
with pysam.TabixFile(src) as tf:
    chrom, rng = region.split(":")
    start, end = (int(x) for x in rng.split("-"))
    out = sys.stdout
    for row in tf.fetch(chrom, start - 1, end):
        out.write(row.rstrip("\r") + "\n")
