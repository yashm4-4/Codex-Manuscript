#!/usr/bin/env python3
"""Final metadata-only documentation probes for signed allele and LD contracts."""
import methods_collect_sources as collector

collector.SOURCES = {
    "hail_bgen_import_contract": "https://hail.is/docs/0.2/methods/impex.html",
    "finngen_pairwise_ld_access": "https://docs.finngen.fi/faq/about-finngen-data/can-i-download-all-pairwise-ld-data-across-the-genome-at-once",
    "susie_rss_source": "https://raw.githubusercontent.com/stephenslab/susieR/master/R/susie_rss.R",
    "finngen_ld_server_source": "https://raw.githubusercontent.com/FINNGEN/ld_server/master/ld_server.py",
}

if __name__ == "__main__":
    collector.main()
