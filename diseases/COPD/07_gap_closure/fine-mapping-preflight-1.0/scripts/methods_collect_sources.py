#!/usr/bin/env python3
"""Acquire bounded public method/LD documentation; never fetch LD matrix parts.

This is a metadata-access preflight, not fine-mapping or LD calculation. Records
the exact bytes returned (including error bodies), response metadata and SHA256.
Files are exclusively created and prior stages are never touched.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
from urllib.request import Request, urlopen
from urllib.error import HTTPError

STAGE = Path(__file__).resolve().parents[1]
AUDIT = STAGE / "audits/methods_ld"
SOURCES = {
    "panukb_ld_docs": "https://pan.ukbb.broadinstitute.org/docs/ld/index.html",
    "panukb_overview_docs": "https://pan.ukbb.broadinstitute.org/docs/technical-overview/index.html",
    "panukb_hail_docs": "https://pan.ukbb.broadinstitute.org/docs/hail-format/index.html",
    "panukb_flat_docs": "https://pan.ukbb.broadinstitute.org/docs/per-phenotype-files/index.html",
    "panukb_qc_docs": "https://pan.ukbb.broadinstitute.org/docs/qc/index.html",
    "panukb_ld_code": "https://raw.githubusercontent.com/atgu/ukbb_pan_ancestry/master/compute_ld_matrix.py",
    "panukb_ld_resources_code": "https://raw.githubusercontent.com/atgu/ukbb_pan_ancestry/master/resources/ld.py",
    "susie_rss_reference": "https://stephenslab.github.io/susieR/reference/susie_rss.html",
    "susie_diagnostics": "https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html",
    "susie_rss_vignette": "https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html",
    "susie_get_cs_reference": "https://stephenslab.github.io/susieR/reference/susie_get_cs.html",
    "susie_rss_paper": "https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010299",
    "finemap_official": "https://christianbenner.com/",
    "finngen_finemapping_docs": "https://finngen.gitbook.io/documentation/methods/finemapping",
    "finngen_pipeline_readme": "https://raw.githubusercontent.com/FINNGEN/finemapping-pipeline/master/README.md",
    "gbmi_slalom_paper": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9839193/",
    "gbmi_slalom_europepmc_xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9839193/fullTextXML",
    "gbmi_slalom_github": "https://api.github.com/repos/mkanai/slalom",
    "gbmi_resources": "https://www.globalbiobankmeta.org/resources",
    "reference_size_paper": "https://www.ukbiobank.ac.uk/publications/prospects-of-fine-mapping-trait-associated-genomic-regions-by-using-summary-statistics-from-genome-wide-association-studies/",
    "thousand_genomes_igsr": "https://www.internationalgenome.org/category/phase-3/",
    "bbj_pheweb_downloads": "https://pheweb.jp/downloads",
    "bbj_ld_github": "https://api.github.com/search/repositories?q=BioBank+Japan+LD+matrix",
    **{f"panukb_{pop}_blockmatrix_metadata": f"https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.{pop}.ldadj.bm/metadata.json" for pop in ("AFR", "AMR", "CSA", "EAS", "EUR", "MID")},
    **{f"panukb_{pop}_variant_metadata": f"https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.{pop}.ldadj.variant.ht/metadata.json.gz" for pop in ("EAS", "EUR")},
}


def now():
    return datetime.now(timezone.utc).isoformat()


def acquire(item, maximum):
    name, url = item
    raw_path = AUDIT / "sources" / (name + ".body")
    metadata_path = AUDIT / "sources" / (name + ".json")
    if raw_path.exists() or metadata_path.exists():
        raise RuntimeError("Refuse overwrite of existing source " + name)
    rec = {"source_id": name, "requested_url": url, "started_utc": now(),
           "request_method": "GET", "request_headers": {"User-Agent": "COPD-scientific-preflight/1.0 (public metadata audit)", "Accept-Encoding": "identity"},
           "maximum_response_bytes": maximum, "LD_matrix_values_requested": False}
    begin = time.monotonic()
    body = b""
    try:
        request = Request(url, headers=rec["request_headers"])
        try:
            response = urlopen(request, timeout=45)
        except HTTPError as error:
            response = error
        with response:
            rec.update(http_status=response.code, final_url=response.geturl(), response_headers=dict(response.headers.items()))
            content = response.read(maximum + 1)
        body = content[:maximum]
        rec["truncated"] = len(content) > maximum
        rec["status"] = "HTTP_SUCCESS" if 200 <= rec["http_status"] < 300 else "HTTP_ERROR"
    except Exception as error:
        rec.update(status="TRANSPORT_ERROR", error_type=type(error).__name__, error=str(error), truncated=False)
    rec.update(completed_utc=now(), wall_seconds=time.monotonic() - begin,
               saved_bytes=len(body), saved_sha256=hashlib.sha256(body).hexdigest(),
               saved_path=str(raw_path.relative_to(STAGE)),
               interpretation="Exact returned bytes only; HTTP success does not establish semantic usability, schema correctness, locus coverage or fine-mapping eligibility.")
    with raw_path.open("xb") as handle:
        handle.write(body)
    with metadata_path.open("x") as handle:
        json.dump(rec, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return rec


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", nargs="+")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--max-bytes", type=int, default=4000000)
    args = parser.parse_args()
    names = args.only or list(SOURCES)
    if not args.run_id.replace("_", "").replace("-", "").isalnum():
        raise ValueError("Unsafe run ID")
    output = AUDIT / ("source_collection_" + args.run_id + ".json")
    if output.exists():
        raise RuntimeError("Refuse overwrite of collection receipt")
    (AUDIT / "sources").mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        records = list(pool.map(lambda name: acquire((name, SOURCES[name]), args.max_bytes), names))
    receipt = {"completed_utc": now(), "stage": STAGE.name, "scope": "public methods/LD documentation and matrix metadata only",
               "scientific_computation": False, "records": records,
               "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    with output.open("x") as handle:
        json.dump(receipt, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"sources": len(records), "statuses": {name: sum(r["status"] == name for r in records) for name in {r["status"] for r in records}},
                      "receipt": str(output.relative_to(STAGE))}, indent=2))


if __name__ == "__main__":
    main()
