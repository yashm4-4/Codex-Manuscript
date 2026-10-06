#!/usr/bin/env python3
"""Public documentation follow-up; no LD matrix values or GWAS rows requested."""
import methods_collect_sources as collector

collector.SOURCES = {
    "panukb_genotype_resource_code": "https://raw.githubusercontent.com/atgu/ukbb_pan_ancestry/master/resources/genotypes.py",
    "susie_summarize_reference": "https://stephenslab.github.io/susieR/reference/susie_summarize.html",
    "finemap_official_www": "https://www.christianbenner.com/",
    "finemap_official_http": "http://www.christianbenner.com/",
    "reference_size_europepmc_xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5630179/fullTextXML",
    "finngen_methods_markdown": "https://finngen.gitbook.io/documentation/methods/finemapping.md",
    "finngen_ld_browser": "https://docs.finngen.fi/working-outside-the-sandbox/ld-browser",
    "finngen_ld_faq": "https://docs.finngen.fi/faq/about-finngen-data",
    "finngen_ld_server_readme": "https://raw.githubusercontent.com/FINNGEN/ld_server/master/README.md",
    "bbj_finemapping_availability": "https://raw.githubusercontent.com/mkanai/finemapping-insights/master/README.md",
    "bbj_finemapping_methods": "https://www.medrxiv.org/content/10.1101/2021.09.03.21262975v1.full",
    "bbj_hum0197_availability": "https://humandbs.dbcls.jp/en/hum0197-latest",
    "gnomad_ld_docs": "https://gnomad.broadinstitute.org/news/2018-10-gnomad-v2-1/",
    "onekg_sample_panel": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel",
    "gbmi_resources_current": "https://www.globalbiobankmeta.org/copy-of-participating-groups",
}

if __name__ == "__main__":
    collector.main()
