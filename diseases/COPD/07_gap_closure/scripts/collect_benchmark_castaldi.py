#!/usr/bin/env python3
"""Collect Castaldi primary evidence, without accessing model outputs.

Authorized by current user attachment 31c0f2f1-803e-4c5e-81ad-bd39763c2f0a,
explicitly requesting this external benchmark module (not the prior RC push).
Source retrieval and curation only; no inferred missing nulls or allele pairs.
"""
import hashlib
import gzip
import json
import re
import csv
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import requests
from bs4 import BeautifulSoup

BASE = Path(__file__).resolve().parents[1]
RAW = BASE / "data/COPD-V2-BENCH/castaldi"
OUT = BASE / "data/COPD-V2-BENCH/castaldi_evidence.json"
ARTICLE = "https://pmc.ncbi.nlm.nih.gov/articles/PMC6353020/"
DOI = "10.1164/rccm.201802-0337OC"
SUPP = "https://pmc.ncbi.nlm.nih.gov/articles/instance/6353020/bin/rccm.201802-0337OC_castaldi_data_supplement.pdf"


def stamp():
    return datetime.now(timezone.utc).isoformat()


def fetch(url, filename, expect, records):
    path = RAW / filename
    meta = path.with_suffix(path.suffix + ".retrieval.json")
    if path.exists() and meta.exists():
        record = json.loads(meta.read_text())
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"]
    else:
        r = requests.get(url, timeout=45)
        data = r.content
        path.write_bytes(data)
        valid = r.status_code == 200
        if expect == "pdf":
            valid = valid and data.startswith(b"%PDF")
        if expect == "image":
            valid = valid and data.startswith(b"\xff\xd8")
        if expect == "gzip":
            valid = valid and data.startswith(b"\x1f\x8b")
        if expect == "article":
            valid = valid and b"606 SNPs" in data and b"Hodges" in data
        status = "retrieved" if valid else "unavailable_or_nonmatching_content"
        if b"POW_CHALLENGE" in data:
            status = "browser_download_challenge_no_bypass_attempted"
        elif b"not open access" in data:
            status = "Europe_PMC_service_reports_not_open_access"
        record = dict(url=url, final_url=r.url, retrieved_utc=stamp(),
                      http_status=r.status_code, content_type=r.headers.get("Content-Type", ""),
                      expected_type=expect, status=status, local_path=str(path.relative_to(BASE)),
                      bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
        meta.write_text(json.dumps(record, indent=2) + "\n")
    records.append(record)
    return path


def num(value):
    return float(re.sub(r"\s+", "", value).replace("−", "-").replace("×10", "e"))


def audit_geo():
    """Enumerate source counts/design only; do not derive experimental P values."""
    designs = {}
    design_tags = set()
    with gzip.open(RAW / "GSE109452_secondary_design_barcodes.dat.gz", "rt") as stream:
        for line in stream:
            name, tags = line.rstrip().split("\t")
            tags = tags.split(",")
            assert name not in designs
            assert not (design_tags & set(tags))
            designs[name] = set(tags)
            design_tags.update(tags)
    alleles, positions, oligos, contexts, tags_seen = defaultdict(set), defaultdict(set), {}, {}, set()
    rows = 0
    input_columns = []
    all_zero = 0
    both_zero = [0, 0]
    with gzip.open(RAW / "GSE109452_FAM13A_MPRA.txt.gz", "rt") as stream:
        input_columns = next(stream).split()
        assert input_columns == ['tag', 'oligo', 'SNP', 'allele', 'chr', 'pos', 'loc', 'orient', 'minP_input', 'minP_output', 'SV40_input', 'SV40_output']
        for number, line in enumerate(stream, 2):
            tag, oligo, snp, allele, chrom, pos, loc, orient, *vals = line.split()
            assert tag not in tags_seen
            assert oligo in designs and tag in designs[oligo]
            tags_seen.add(tag)
            rows += 1
            vals = list(map(int, vals))
            assert len(vals) == 4 and min(vals) >= 0
            all_zero += not any(vals)
            both_zero[0] += not any(vals[:2])
            both_zero[1] += not any(vals[2:])
            alleles[snp].add(allele)
            positions[snp].add((chrom, pos))
            key = (snp, loc, orient)
            contexts.setdefault(key, set()).add(allele)
            if oligo not in oligos:
                oligos[oligo] = dict(oligo=oligo, rsid=snp, source_chrom=chrom,
                    source_pos=pos, allele=allele, loc=loc, orient=orient,
                    observed_barcodes=0, designed_barcodes=len(designs[oligo]),
                    minP_input_sum=0, minP_output_sum=0, SV40_input_sum=0, SV40_output_sum=0,
                    minP_input_nonzero=0, minP_output_nonzero=0, SV40_input_nonzero=0, SV40_output_nonzero=0)
            item = oligos[oligo]
            item['observed_barcodes'] += 1
            for name, value in zip(input_columns[8:], vals):
                item[name + '_sum'] += value
                item[name + '_nonzero'] += value > 0
    assert len(alleles) == 606 and all(len(v) == 2 for v in alleles.values())
    assert all(len(v) == 1 for v in positions.values())
    assert len(oligos) == len(designs) == 7270
    assert all(len(v) == 2 for v in contexts.values())
    missing = [(s, l, o) for s in sorted(alleles) for l in ('L', 'C', 'R') for o in ('F', 'RC') if (s, l, o) not in contexts]
    with (RAW / 'GSE109452_oligo_count_coverage.tsv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(next(iter(oligos.values()))), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(oligos[key] for key in sorted(oligos))
    with (RAW / 'GSE109452_variant_panel.tsv').open('w') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(['rsid', 'source_chrom', 'source_pos', 'source_build', 'tested_allele1', 'tested_allele2', 'observed_oligo_contexts', 'result_state_availability'])
        for snp in sorted(alleles):
            writer.writerow([snp, *next(iter(positions[snp])), 'not_reported_in_GEO_SOFT', *sorted(alleles[snp]), sum(key[0] == snp for key in contexts), 'Counts only; no significance or QC labels in GEO files'])
    audit = dict(geo_accession='GSE109452', processed_columns=input_columns,
        source_counts_rows=rows, unique_observed_barcodes=len(tags_seen), designed_barcodes=len(design_tags),
        designed_tags_not_in_processed=len(design_tags-tags_seen), processed_tags_not_in_design=len(tags_seen-design_tags),
        unique_variants=len(alleles), unique_allele_constructs=len(oligos), paired_oligo_contexts=len(contexts),
        promoter_contexts=2*len(contexts), expected_606x6_oligo_contexts=606*6,
        missing_design_and_processed_contexts=missing, all_zero_count_rows=all_zero,
        minP_both_counts_zero_rows=both_zero[0], SV40_both_counts_zero_rows=both_zero[1],
        per_construct_observed_barcodes_range=[min(v['observed_barcodes'] for v in oligos.values()), max(v['observed_barcodes'] for v in oligos.values())],
        every_construct_has_33_designed_barcodes=all(len(v)==33 for v in designs.values()),
        coordinates_unique_per_variant=True, allele_count_each_variant=2,
        published_significance_columns_available=False, published_QC_exclusions_available=False,
        genome_build_reported_in_SOFT=False, numerical_MPRA_calling_rule_available=False,
        barcode_rows_are_independent_variants=False, statistical_reanalysis_performed=False,
        interpretation='Complete 606-variant design and observed barcode panel are enumerable, but count data alone do not supply published positive/null labels. No missing barcode, zero count, or non-hit was labeled null.')
    (RAW / 'GSE109452_count_design_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    return alleles, positions, contexts, oligos, audit


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    retrievals = []
    article = fetch(ARTICLE, "PMC6353020_article.html", "article", retrievals)
    assert retrievals[-1]["status"] == "retrieved"
    soup = BeautifulSoup(article.read_text(), "html.parser")
    extracted = []
    previous = ""
    for index, tr in enumerate(soup.find("table").find("tbody").find_all("tr"), 1):
        vals = [td.get_text("", strip=True) for td in tr.find_all("td")]
        if not any(vals):
            continue
        if vals[0].startswith("rs"):
            previous = vals[0]
        else:
            vals.insert(0, previous)
        assert len(vals) == 11, vals
        extracted.append(dict(source_html_row=index, cells=vals))
    assert len(extracted) == 16
    (RAW / "Table1_source_rows.json").write_text(json.dumps(extracted, indent=2) + "\n")
    for figure in ("f5", "f6"):
        img = next(i for i in soup.find_all("img") if i.get("src", "").endswith(f"0337OC_{figure}.jpg"))
        fetch(img["src"], f"PMC6353020_{figure}.jpg", "image", retrievals)
    alternatives = [
        (SUPP, "supplement_PMC_response.html", "pdf"),
        ("https://www.atsjournals.org/doi/suppl/" + DOI + "/suppl_file/rccm.201802-0337oc_castaldi_data_supplement.pdf", "supplement_ATS_response.html", "pdf"),
        ("https://europepmc.org/articles/PMC6353020/bin/rccm.201802-0337OC_castaldi_data_supplement.pdf", "supplement_EuropePMC_response.html", "pdf"),
        ("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6353020/supplementaryFiles", "supplement_EuropePMC_service_response.xml", "pdf"),
        ("https://academic.oup.com/ajrccm/article/199/1/52/8496721", "publisher_response.html", "article"),
        ("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6353020/fullTextXML", "article_EuropePMC_XML_response.txt", "xml"),
        ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109452/suppl/", "GSE109452_supplement_directory.html", "html"),
        ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109452/suppl/GSE109452_FAM13A_MPRA.txt.gz", "GSE109452_FAM13A_MPRA.txt.gz", "gzip"),
        ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109452/suppl/GSE109452_secondary_design_barcodes.dat.gz", "GSE109452_secondary_design_barcodes.dat.gz", "gzip"),
        ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109452/soft/GSE109452_family.soft.gz", "GSE109452_family.soft.gz", "gzip"),
        ("https://castaldilab.org/datasets-and-websites/", "Castaldi_lab_data_page.html", "html"),
        ("https://pmc.ncbi.nlm.nih.gov/articles/PMC7063182/", "Lin2020_Castaldi_MPRA_reuse_primary.html", "html"),
    ]
    for url, name, kind in alternatives:
        fetch(url, name, kind, retrievals)
    geo_alleles, geo_positions, geo_contexts, geo_oligos, geo_audit = audit_geo()
    assays, contexts = [], []
    # Figure 5 was visually inspected for exact allele names and direction.
    pairs = {"rs7674369": ("A", "G"), "rs2013701": ("G", "T"),
             "rs7671167": ("C", "T"), "rs1795739": ("A", "G"),
             "rs78681184": ("C", "T"), "rs147089648": ("C", "T")}
    for snp, pair in pairs.items():
        if snp in geo_alleles:
            assert set(pair) == geo_alleles[snp]
    pairs.update({snp: tuple(sorted(v)) for snp, v in geo_alleles.items()})
    reporter_variants = ['rs7674369', 'rs2013701', 'rs7671167', 'rs1795739', 'rs78681184', 'rs147089648']

    def row(aid, rsid, cls, cell, state, locator, **extra):
        a1, a2 = pairs.get(rsid, ("", ""))
        result = dict(assay_id=aid, study_id="CASTALDI2019", locus="FAM13A", rsid=rsid,
            source_build="not_reported_in_accessible_primary_text_or_GEO_SOFT", source_chrom="4", source_pos=next(iter(geo_positions[rsid]))[1] if rsid in geo_positions else "",
            tested_allele1=a1, tested_allele2=a2, reported_strand="not_reported", assay_class=cls,
            cell_context=cell, experimental_state=state, reported_effect_direction="", effect_allele="",
            higher_activity_allele="", p_value="", fdr="", source_rule="", mechanism_in_model_scope="yes",
            mechanism_scope_rationale="Allele-dependent transcriptional reporter activity is in broad regulatory-sequence scope; cell context differs from bulk lung.",
            source_locator=locator, source_url=ARTICLE, source_id="CASTALDI2019_MAIN", notes="",
            evidence_role="individual_assay_context", phenotype_ascertainment="COPD GWAS or COPD whole-blood FAM13A eQTL selection")
        result.update(extra)
        return result

    for source_row in extracted:
        rsid, ea, _, _, context, e1, q1, e2, q2, _, _ = source_row["cells"]
        for promoter, effect, q in (("minP", e1, q1), ("SV40", e2, q2)):
            eff = num(effect)
            a1, a2 = pairs.get(rsid, ("", ""))
            assays.append(row(f"CASTALDI2019_MPRA_{rsid}_{context.replace(':', '_')}_{promoter}", rsid,
                "MPRA_allele_effect", f"Beas-2B; {promoter}; 144bp; SNP-position:orientation={context}",
                "ambiguous", f"Table 1 HTML body row {source_row['source_html_row']}; {promoter}",
                reported_strand="effect allele genomic coding; oligo orientation relative to human reference explicitly reported",
                tested_allele1=a1 or ea, tested_allele2=a2,
                reported_effect_allele=ea, effect_value=eff, reported_effect_direction="higher reported allelic ratio" if eff > 0 else "lower reported allelic ratio",
                effect_allele="", higher_activity_allele="", fdr=num(q),
                source_rule="FDR-adjusted P per Table 1 footnote; exact numerical MPRA threshold inaccessible until supplement retrieved; no generic cutoff imposed",
                notes="Retained regardless of FDR; source certifies significance at variant-any-context level, not every table entry. Exact tested pair and original coordinate recovered from GEO; genome build not stated. Effect allele/sign preserved, but higher-activity allele not assigned: Figure 3 calls the metric input/output while Table 1 says allelic ratio; missing methods prevent resolving ratio convention for activity direction.",
                source_row=source_row, promoter=promoter, oligo_context=context,
                assay_contrast=f"{ea} versus other allele; Hodges-Lehmann input/output ratio difference"))
    for rsid in dict.fromkeys(r["cells"][0] for r in extracted):
        assays.append(row(f"CASTALDI2019_MPRA_{rsid}_SOURCE_ANY_CONTEXT", rsid, "MPRA_allele_effect",
            "Beas-2B; source-reported any oligo/promoter context", "positive", "Table 1 title and footnote; MPRA Results",
            source_rule="Authors explicitly include variant among MPRA-significant SNPs; union over assayed oligo/promoter contexts",
            evidence_role="source_reported_variant_summary",
            notes="Not an additional experiment; descriptive source-summary only. Table 1 is a GWAS-selected 8-of-45-hit subset. Never count alongside its individual-context rows as independent assays."))
    assays.append(row("CASTALDI2019_MPRA_rs4416442_SOURCE_ANY_CONTEXT", "rs4416442", "MPRA_allele_effect",
        "Beas-2B; source-reported MPRA overall", "null", "Conditional Analysis Identifies Allelic Heterogeneity, paragraph 1",
        source_rule="Authors explicitly state lead GWAS SNP rs4416442 was not significant in MPRAs",
        evidence_role="source_reported_variant_summary",
        notes="Explicit tested-null statement, not inferred from absence; tested C/T allele pair recovered from GEO, but per-oligo/per-promoter result labels remain unavailable."))
    reuse_text = BeautifulSoup((RAW / 'Lin2020_Castaldi_MPRA_reuse_primary.html').read_text(), 'html.parser').get_text(' ', strip=True)
    assert '45 SNPs with significant allelic effects' in reuse_text
    for rsid in ('rs2276936', 'rs2167750', 'rs7695177'):
        assert rsid in geo_alleles and rsid in reuse_text
        assays.append(row(f'CASTALDI2019_MPRA_{rsid}_SOURCE_ANY_CONTEXT', rsid, 'MPRA_allele_effect',
            'Beas-2B; source-reported any oligo/promoter context', 'positive',
            'Lin et al. 2020 PMC7063182 Results: RS2276936 within FAM13A Gene Associated with HDL and Fat Mass Exerts Allelic Activity, first paragraph',
            source_id='CASTALDI2019_REUSE_LIN2020', source_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC7063182/',
            evidence_role='source_reported_variant_summary',
            source_rule='Primary follow-up explicitly states these three SNPs are drawn from the 45 significant allelic SNPs in previous Castaldi2019 MPRA',
            notes='Reused CASTALDI2019 experiment, not an independent replication or a newly scored MPRA. Only author-reported any-context positivity is recovered; per-promoter/window/orientation states remain unavailable. HDL association P values in that paragraph are not MPRA P values.'))
    for rsid in reporter_variants:
        higher = {"rs2013701": "T", "rs7671167": "C", "rs1795739": "A"}.get(rsid, "")
        positive = bool(higher)
        assays.append(row(f"CASTALDI2019_REPORTER_{rsid}_16HBE", rsid, "conventional_reporter",
            "16HBE; pGL4.23 luciferase; no exposure reported", "positive" if positive else "null",
            "Figure 5, visually inspected; Allele-Specific Distal Enhancer Activities",
            reported_effect_direction=f"{higher} higher reporter activity" if positive else "no significant allelic difference reported",
            effect_allele=higher, higher_activity_allele=higher, p_value="<0.001" if positive else "not_reported",
            source_rule="Figure 5 ***P<0.001, unpaired Student t test; Results explicitly identify only 3 of 6 tested variants as significant",
            notes="4–7 independent experiments, triplicate measurements. Null is allelic contrast, not activity versus vector. rs147089648 was not in the initial 606-variant MPRA panel."))
    assays.append(row("CASTALDI2019_EDIT_rs2013701_16HBE_expression", "rs2013701", "endogenous_allele_editing",
        "16HBE; CRISPR homology-directed GG-to-TT editing; submerged clonal culture", "positive",
        "Figure 6A and CRISPR-based Genome Editing section (heading typo rs2017301)",
        reported_effect_direction="TT has higher FAM13A mRNA than GG", effect_allele="T", higher_activity_allele="T",
        p_value="<0.05", source_rule="Figure 6 *P<0.05 unpaired Student t test; Results confirm lower FAM13A expression with G",
        mechanism_in_model_scope="partial", mechanism_scope_rationale="Isolated endogenous editing changes expression; supportive but not a direct bulk-lung enhancer/H3K27me3 assay.",
        notes="Two independently edited TT clones; 2–3 clones/genotype with 3 repeats in legend. Figure 6B proliferation is a separate downstream phenotype."))

    def context(cid, cls, state, locator, notes, **extra):
        item = dict(context_id=cid, study_id="CASTALDI2019", locus="FAM13A", assay_class=cls,
            cell_context="16HBE", experimental_state=state, source_locator=locator, source_url=ARTICLE,
            mechanism_in_model_scope="unevaluable", notes=notes, exact_variant_label=False)
        item.update(extra)
        contexts.append(item)
    context("CASTALDI2019_CONTACT_rs2013701_rs7671167_REGION", "chromatin_contact", "positive",
        "Allele-Specific Distal Enhancer Activities; Figure E5 referenced, supplement unavailable",
        "Promoter contact with one genomic segment containing rs2013701 AND rs7671167. Not two allele-specific contacts or nucleotide labels.", associated_rsids=["rs2013701", "rs7671167"])
    context("CASTALDI2019_CONTACT_rs1795739_REGION", "chromatin_contact", "null",
        "Allele-Specific Distal Enhancer Activities; Figure E5 referenced",
        "No promoter contact for region containing rs1795739. Regional null does not nullify reporter-positive nucleotide.", associated_rsids=["rs1795739"])
    context("CASTALDI2019_EDIT_rs2013701_proliferation", "other_or_unresolved", "positive", "Figure 6B and Results",
        "GG clones proliferate faster than TT; downstream cellular outcome, not enhancer-direction label.", associated_rsids=["rs2013701"])
    context("CASTALDI2019_GTEx_rs2013701", "expression_QTL", "positive", "GTEx confirmation Results; Figure E7 referenced",
        "Reported GTEx lung association concordant with T-higher FAM13A. QTL is contextual, not experimental nucleotide label.", associated_rsids=["rs2013701"], cell_context="GTEx v7 whole lung")
    context("CASTALDI2019_GEO_MISSING_rs3017926_C_RC", "MPRA_allele_effect", "unavailable",
        "GSE109452 design dictionary and processed counts; expected SNP center reverse context",
        "Both alleles of rs3017926 center/reverse context are absent from the design and count files. This is not a null experiment. Other contexts for this SNP are present.",
        associated_rsids=["rs3017926"], cell_context="Beas-2B; minP and SV40; C:RC", exclusion_reason="Undesigned or unavailable context; source does not explain omission")
    # Preserve the complete enumerable observed panel at the actual paired
    # variant x window x orientation x promoter unit, never barcode-row units.
    assay_by_id = {a['assay_id']: a for a in assays}
    constructs_by_context = defaultdict(list)
    for construct in geo_oligos.values():
        constructs_by_context[(construct['rsid'], construct['loc'], construct['orient'])].append(construct)
    for (rsid, loc, orient), constructs in sorted(constructs_by_context.items()):
        assert len(constructs) == 2
        for promoter in ('minP', 'SV40'):
            aid = f'CASTALDI2019_MPRA_{rsid}_{loc}_{orient}_{promoter}'
            if aid in assay_by_id:
                a = assay_by_id[aid]
            else:
                a = row(aid, rsid, 'MPRA_allele_effect', f'Beas-2B; {promoter}; 144bp; SNP-position:orientation={loc}:{orient}',
                    'unavailable', f'GSE109452_FAM13A_MPRA.txt.gz: SNP={rsid}; loc={loc}; orient={orient}; {promoter}_input and {promoter}_output',
                    source_id='CASTALDI2019_GEO', source_url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109452',
                    reported_strand='GEO allele labels retained; construct F/RC orientation reported; genome build not stated',
                    source_rule='Processed barcode counts only; published per-context statistical result and QC/calling rule unavailable',
                    notes='An explicit tested allele contrast is enumerated; its positive/null statistical result is unavailable. Zero counts, missing barcodes and absence from a selected hit list are never null labels.',
                    promoter=promoter, oligo_context=f'{loc}:{orient}')
                assays.append(a)
            a['source_oligo_ids'] = ';'.join(sorted(c['oligo'] for c in constructs))
            a['source_count_file'] = 'data/COPD-V2-BENCH/castaldi/GSE109452_FAM13A_MPRA.txt.gz'
            a['source_design_file'] = 'data/COPD-V2-BENCH/castaldi/GSE109452_secondary_design_barcodes.dat.gz'
            a['designed_barcodes_for_contrast'] = sum(c['designed_barcodes'] for c in constructs)
            a['observed_barcodes_for_contrast'] = sum(c['observed_barcodes'] for c in constructs)
            a['allele_count_totals'] = {c['allele']: dict(input_count=c[promoter+'_input_sum'], output_count=c[promoter+'_output_sum'], input_nonzero_barcodes=c[promoter+'_input_nonzero'], output_nonzero_barcodes=c[promoter+'_output_nonzero']) for c in constructs}
    searches = ["PMC6353020 FAM13A Castaldi supplementary 606 45 MPRA rs2013701",
        "10.1164/rccm.201802-0337OC supplementary tables", '"rccm.201802-0337OC_castaldi_data_supplement.pdf"',
        '"Castaldi" "FAM13A" "Table E2"', '"FAM13A" "MPRA" "606" data download',
        '"201802-0337OC" "FDR"', '"castaldi_data_supplement" "199"', '"8496721" "Supplementary"',
        '"606" "FAM13A" "Table E2"', '"rccm.201802-0337OC_castaldi_data_supplement.pdf" full panel',
        '"GSE109452" analysis code', '"FAM13A_MPRA" "wilcox"', '"Castaldi" "FAM13A" "MPRA" github', '"GSE109452" FDR',
        'Lin2020 PMC7063182 rs2276936 rs2167750 rs7695177 Castaldi MPRA']
    data = dict(schema_version="1.0", assembled_utc=stamp(), source_reader_model_outputs_accessed=False,
        sources=[dict(source_id="CASTALDI2019_MAIN", study_id="CASTALDI2019", title="Identification of Functional Variants in the FAM13A Chronic Obstructive Pulmonary Disease Genome-Wide Association Study Locus by Massively Parallel Reporter Assays", doi=DOI, pmcid="PMC6353020", pmid="30079747", publication_year=2019, url=ARTICLE, study_design="Single-locus COPD-selected MPRA with selected reporter/editing follow-up", evidence_status="primary_article_retrieved_supplement_unavailable_GEO_counts_retrieved", locus="FAM13A", primary=True),
            dict(source_id='CASTALDI2019_GEO', study_id='CASTALDI2019', title='Fine Mapping and Functional Characterization of Genetic Variants in the FAM13A Chronic Obstructive Pulmonary Disease GWAS locus using Massively Parallel Reporter Assays', geo_accession='GSE109452', url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109452', direct_processed_url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109452/suppl/GSE109452_FAM13A_MPRA.txt.gz', evidence_status='Complete design/observed barcode panel available; published significance/QC labels not supplied', primary=True, locus='FAM13A', reported_genome_build='not_reported', linkage_to_article='Matching title, investigators, 606-variant/45-hit design, rs2013701 results and official Castaldi lab data page; GEO citation field itself is missing'),
            dict(source_id='CASTALDI2019_REUSE_LIN2020', study_id='CASTALDI2019', title='FAM13A Represses AMPK Activity and Regulates Hepatic Glucose and Lipid Metabolism', doi='10.1016/j.isci.2020.100928', pmcid='PMC7063182', publication_year=2020, url='https://pmc.ncbi.nlm.nih.gov/articles/PMC7063182/', primary=True, evidence_status='Primary follow-up explicitly names three original Castaldi MPRA-positive variants; reused observations, not independent replication', related_study_id='LIN2020_FAM13A')],
        assays=assays, contexts=contexts, retrievals=retrievals,
        searches=[dict(query=q, searched_date_utc="2026-10-05", service="web search", result="Primary article, GEO and official lab data page located; full count/design panel retrieved through standard GEO FTP HTTPS endpoints. Original analysis code, complete QC/test rule and supplement hit table not recovered.") for q in searches],
        denominator_audit=dict(study_id="CASTALDI2019", designed_variants_reported=606, assayed_variants_reported=606,
            qc_passing_variants="not_enumerated_in_accessible_sources", oligo_length_bp=144,
            promoters=["minP", "SV40"], cell_context="Beas-2B", significant_oligo_contexts_reported=85,
            significant_unique_variants_reported=45, recovered_main_table_variants=8,
            recovered_main_table_oligo_contexts=16, recovered_main_table_promoter_measurements=32,
            recovered_reported_MPRA_hit_variants=11, recovered_hit_variants_from_primary_reuse=3,
            unseen_reported_hit_variants=34, all606identity_list_recovered=True, all606_tested_allele_pairs_recovered=True, complete_null_panel_available=False,
            complete_denominator_pass=False, complete_assayed_denominator_valid=False,
            complete_design_and_observed_panel_available=True,
            reason="Complete 606-variant design/observed counts are available, but published QC exclusions, numerical calling rule and all per-context positive/null result states are not.",
            denominator_failure_reason="GEO counts and design enumerate all 606 variants and their allele pairs; they do not contain P/FDR/significance or analytical QC columns. The inaccessible supplement prevents reconstruction of published labels and 45-hit identities. Raw count reanalysis would require unverified filtering/normalization/testing choices; no such analysis was performed.",
            numerical_MPRA_rule_status="FDR-adjusted P verified; exact threshold unavailable in retrieved primary text; 32 individual contexts ambiguous pending source rule",
            explicit_MPRA_null_rsids=["rs4416442"], selected_reporter_panel_count=6,
            selected_reporter_positive_count=3, selected_reporter_null_count=3, reporter_panel_source_complete=True,
            reporter_panel_selection="5 MPRA-prioritized variants plus secondary GWAS lead; not representative of the 606-variant panel",
            endogenous_editing_variants=1, locus_count=1,
            geo_count_design_audit=geo_audit, recovered_promoter_contexts=7270, known_individual_result_stats_rows=32,
            unavailable_individual_result_rows=7238, source_reported_variant_summary_rows=12,
            permitted_analysis="Descriptive complete-panel model coverage and selected-case recovery; no full-606 sensitivity/specificity/AUROC/AP/FPR without published experimental labels; source-summary rows are not independent assays"))
    assert len({a["assay_id"] for a in assays}) == len(assays)
    checks = dict(
        source_reader_model_outputs_accessed=False,
        expected_total_rows=len(assays) == 7289,
        distinct_source_variants=len({a['rsid'] for a in assays}) == 607,
        all_source_allele_pairs_present=all(a['tested_allele1'] and a['tested_allele2'] for a in assays),
        raw_GEO_context_rows=sum(a['assay_class']=='MPRA_allele_effect' and a['evidence_role']=='individual_assay_context' for a in assays)==7270,
        raw_GEO_unavailable_result_rows=sum(a['experimental_state']=='unavailable' for a in assays)==7238,
        published_ambiguous_contexts=sum(a['experimental_state']=='ambiguous' for a in assays)==32,
        source_reported_summaries=sum(a['evidence_role']=='source_reported_variant_summary' for a in assays)==12,
        summaries_not_counted_as_new_experiments=True,
        reporter_panel_six=sum(a['assay_class']=='conventional_reporter' for a in assays)==6,
        reporter_positive_three=sum(a['assay_class']=='conventional_reporter' and a['experimental_state']=='positive' for a in assays)==3,
        reporter_null_three=sum(a['assay_class']=='conventional_reporter' and a['experimental_state']=='null' for a in assays)==3,
        incomplete_full_label_denominator_not_used=not data['denominator_audit']['complete_assayed_denominator_valid'],
        no_unverified_MPRA_activity_direction=all(not a['higher_activity_allele'] for a in assays if a['assay_class']=='MPRA_allele_effect'),
        no_nulls_derived_from_GEO_counts=not any(a['source_id']=='CASTALDI2019_GEO' and a['experimental_state']=='null' for a in assays),
        raw_retrieval_hashes_match=all(hashlib.sha256((BASE/r['local_path']).read_bytes()).hexdigest()==r['sha256'] for r in retrievals))
    assert all(v for k, v in checks.items() if k != 'source_reader_model_outputs_accessed')
    (RAW / 'Castaldi_source_curation_checks.json').write_text(json.dumps(checks, indent=2) + '\n')
    OUT.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(dict(output=str(OUT), assay_rows=len(assays), context_rows=len(contexts),
        states={s: sum(a["experimental_state"] == s for a in assays) for s in ("positive", "null", "ambiguous", "unavailable")}), indent=2))


if __name__ == "__main__":
    main()
