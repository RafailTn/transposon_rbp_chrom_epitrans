"""Search and download processed contact matrices from the 4D Nucleome portal.

    # what HiChIP mcools exist for K562?
    python src/hic/fetch_4dn.py search --experiment-type HiChIP --biosource K562

    # in situ Hi-C, any cell type, list only
    python src/hic/fetch_4dn.py search --experiment-type "in situ Hi-C"

    # grab one
    python src/hic/fetch_4dn.py get 4DNFI1ABCDEF --out data/hic

Released 4DN files need no credentials. Files still under embargo return 403 and
require a portal account plus an access key; pass --key KEY:SECRET for those.

WHICH FILE TO TAKE
The .mcool, essentially always. It is multi-resolution and already ICE-balanced,
so the weights this pipeline reads are present. The .pairs are 5-50 GB and are
only worth downloading if you intend to re-call loops or reprocess alignments for
multi-mapping reads -- which is the one thing that would let young subfamilies
(AluY, L1HS, SVA-E/F) into the analysis at all, and is a much bigger project than
this one.
"""
import argparse
import os
import sys

BASE = "https://data.4dnucleome.org"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)

    s = sp.add_parser("search", help="list matching processed files")
    s.add_argument("--experiment-type", default=None,
                   help='exact portal value: "HiChIP", "in situ Hi-C", "Micro-C", '
                        '"PLAC-seq", "in situ ChIA-PET", ... use --show-vocab to list')
    s.add_argument("--biosource", default=None,
                   help='SUBSTRING, case-insensitive: "k562" matches the portal\'s '
                        '"K562 (Tier 2)". Exact matching here is a trap -- '
                        'biosource_name is free text, not a controlled vocabulary')
    s.add_argument("--dataset", default=None, help="substring, case-insensitive")
    s.add_argument("--file-format", default="mcool")
    s.add_argument("--assembly", default="GRCh38")
    s.add_argument("--limit", default="all")
    s.add_argument("--show-vocab", action="store_true",
                   help="tally experiment_type and biosource_name over the matching "
                        "files instead of listing them; start here")
    s.add_argument("--filter", action="append", default=[], metavar="KEY=VALUE",
                   help="raw portal facet, repeatable")
    s.add_argument("--key", default=None, help="ACCESS_KEY:SECRET for embargoed data")

    g = sp.add_parser("get", help="download processed files by accession")
    g.add_argument("accession", nargs="+")
    g.add_argument("--out", default="data/hic")
    g.add_argument("--key", default=None)
    return p.parse_args()


def auth(key):
    if not key:
        return None
    if ":" not in key:
        raise SystemExit("--key must be ACCESS_KEY:SECRET")
    return tuple(key.split(":", 1))


def search(a):
    """Server-side facets for the controlled fields, client-side for the free-text ones.

    biosource_name and dataset are free text on the portal ("K562 (Tier 2)", not
    "K562"), so passing them as server facets silently returns zero hits rather
    than erroring. Filter those here on substring instead. experiment_type IS a
    controlled value and is safe to send server-side.
    """
    import requests
    from collections import Counter

    params = [("type", "FileProcessed"), ("format", "json"), ("limit", a.limit)]
    if a.file_format:
        params.append(("file_format.file_format", a.file_format))
    if a.assembly:
        params.append(("genome_assembly", a.assembly))
    if a.experiment_type:
        params.append(("track_and_facet_info.experiment_type", a.experiment_type))
    for f in a.filter:
        if "=" not in f:
            raise SystemExit(f"--filter needs KEY=VALUE, got {f!r}")
        params.append(tuple(f.split("=", 1)))

    r = requests.get(f"{BASE}/search/", params=params, auth=auth(a.key),
                     headers={"Accept": "application/json"}, timeout=180)
    if r.status_code == 404:
        print("no results for that query", file=sys.stderr)
        return
    r.raise_for_status()
    hits = r.json().get("@graph", [])

    def info(h):
        return h.get("track_and_facet_info") or {}

    if a.biosource:
        n = a.biosource.lower()
        hits = [h for h in hits if n in (info(h).get("biosource_name") or "").lower()]
    if a.dataset:
        n = a.dataset.lower()
        hits = [h for h in hits if n in (info(h).get("dataset") or "").lower()]

    if not hits:
        print("no results for that query", file=sys.stderr)
        print("try --show-vocab to see the exact experiment_type / biosource values",
              file=sys.stderr)
        return

    if a.show_vocab:
        for field in ("experiment_type", "biosource_name"):
            print(f"--- {field} ({len(hits)} files) ---")
            for k, v in Counter(info(h).get(field, "?") for h in hits).most_common():
                print(f"  {v:5d}  {k}")
            print()
        return

    print(f"{len(hits)} file(s)\n")
    print(f"{'accession':<14} {'size':>9}  {'assay':<15} {'biosource':<26} dataset")
    print("-" * 110)
    for h in sorted(hits, key=lambda x: -(x.get("file_size") or 0)):
        i = info(h)
        size = h.get("file_size")
        print(f"{h.get('accession', '?'):<14} "
              f"{(f'{size / 1e9:.2f} GB' if size else '?'):>9}  "
              f"{(i.get('experiment_type') or '?')[:15]:<15} "
              f"{(i.get('biosource_name') or '?')[:26]:<26} "
              f"{(i.get('dataset') or '')[:38]}")
    print(f"\ndownload with:  python src/hic/fetch_4dn.py get <accession> --out data/hic")


def get(a):
    import requests

    os.makedirs(a.out, exist_ok=True)
    for acc in a.accession:
        meta = requests.get(f"{BASE}/files-processed/{acc}/", params={"format": "json"},
                            auth=auth(a.key), headers={"Accept": "application/json"},
                            timeout=120)
        if meta.status_code == 403:
            raise SystemExit(f"{acc}: 403 -- embargoed, pass --key ACCESS_KEY:SECRET")
        meta.raise_for_status()
        m = meta.json()
        # open_data_url is the AWS Open Data S3 copy and needs no credentials;
        # every released GRCh38 mcool carries one. Fall back to the portal href,
        # which does go through auth, only if it is missing.
        url = m.get("open_data_url") or (f"{BASE}{m['href']}" if m.get("href") else None)
        if not url:
            raise SystemExit(f"{acc}: no open_data_url or href in metadata")
        total = m.get("file_size")
        dest = os.path.join(a.out, os.path.basename(url.split("?")[0]))

        # Resume rather than restart: these are multi-GB and a dropped connection
        # partway through an mcool is otherwise a full re-download.
        done = os.path.getsize(dest) if os.path.exists(dest) else 0
        if total and done == total:
            print(f"{acc}: already complete ({done / 1e9:.2f} GB) -> {dest}")
            continue
        headers = {"Range": f"bytes={done}-"} if done else {}
        print(f"{acc}: downloading {(total or 0) / 1e9:.2f} GB -> {dest}"
              + (f" (resuming at {done / 1e9:.2f} GB)" if done else ""))

        with requests.get(url, stream=True, auth=auth(a.key),
                          headers=headers, timeout=300) as r:
            r.raise_for_status()
            mode = "ab" if done and r.status_code == 206 else "wb"
            if mode == "wb":
                done = 0
            with open(dest, mode) as fh:
                for chunk in r.iter_content(chunk_size=1 << 22):
                    fh.write(chunk)
                    done += len(chunk)
                    if total:
                        pct = 100 * done / total
                        print(f"\r  {done / 1e9:6.2f} / {total / 1e9:.2f} GB "
                              f"({pct:5.1f}%)", end="", flush=True)
        print(f"\n{acc}: done")


def main():
    a = parse_args()
    (search if a.cmd == "search" else get)(a)


if __name__ == "__main__":
    main()
