#!/usr/bin/env python3
"""Bounded alternative public primary-document endpoints after recorded failures."""
import methods_collect_sources as collector

collector.SOURCES = {
    "reference_size_pmc_html": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5630179/",
    "bbj_finemapping_methods_pdf": "https://www.medrxiv.org/content/10.1101/2021.09.03.21262975v1.full.pdf",
    "bbj_pipeline_readme": "https://raw.githubusercontent.com/mkanai/finemapping-pipeline/master/README.md",
}

if __name__ == "__main__":
    collector.main()
