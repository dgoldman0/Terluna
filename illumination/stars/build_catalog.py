#!/usr/bin/env python3
"""Build the compact bright-star catalogue from the Yale Bright Star Catalogue (BSC5).

Source: Hoffleit D., Warren Jr W.H. (1991), Bright Star Catalogue, 5th Revised Ed.,
NSSDC/ADC, distributed by CDS as VizieR catalogue V/50. The source file's SHA-256 is
pinned; a mismatch stops the build rather than substituting data.

Usage: python illumination/stars/build_catalog.py [--archive catalog.gz]
Writes illumination/stars/bright_stars.json: J2000 right ascension and declination
(radians), V magnitude and B-V colour for every entry with a position and V.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
URL = "https://cdsarc.cds.unistra.fr/ftp/V/50/catalog.gz"
SHA256 = "3dc44b1e90be8fbe5bcc7656032560f51275f985c7e3f783c9028e1838ec7bed"


def field(line: str, a: int, b: int) -> str:
    """Byte columns a..b (1-based, inclusive) of a fixed-width BSC5 record."""
    return line[a - 1:b].strip()


def parse(text: str) -> list[list[float]]:
    stars = []
    for line in text.splitlines():
        line = line.ljust(197)
        ra_h, dec_d, vmag = field(line, 76, 77), field(line, 85, 86), field(line, 103, 107)
        if not ra_h or not dec_d or not vmag:
            continue  # entries without a J2000 position or V (e.g. retracted objects)
        ra = (int(ra_h) + int(field(line, 78, 79)) / 60 + float(field(line, 80, 83)) / 3600) * math.pi / 12
        dec = int(dec_d) + int(field(line, 87, 88)) / 60 + int(field(line, 89, 90)) / 3600
        dec = math.radians(-dec if field(line, 84, 84) == "-" else dec)
        bv = field(line, 110, 114)
        stars.append([round(ra, 6), round(dec, 6), float(vmag), float(bv) if bv else 0.6])
    return stars


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, help="Local copy of catalog.gz")
    args = parser.parse_args()
    raw = args.archive.read_bytes() if args.archive else urllib.request.urlopen(URL, timeout=60).read()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != SHA256:
        raise SystemExit(f"BSC5 hash mismatch ({digest}); the catalogue was not rebuilt.")
    stars = parse(gzip.decompress(raw).decode("ascii"))
    out = {
        "schema": "terluna.illumination.bright-stars/1",
        "source": {"catalogue": "Yale Bright Star Catalogue, 5th Revised Ed. (Hoffleit & Warren 1991)",
                   "distributor": "CDS VizieR V/50", "url": URL, "sha256": SHA256},
        "columns": ["ra_j2000_rad", "dec_j2000_rad", "v_mag", "b_minus_v"],
        "notes": "B-V is 0.6 where the catalogue gives none. Proper motion is not applied.",
        "stars": stars,
    }
    (HERE / "bright_stars.json").write_text(json.dumps(out, separators=(",", ":")) + "\n")
    print(f"{len(stars)} stars, V {min(s[2] for s in stars):.2f} to {max(s[2] for s in stars):.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
