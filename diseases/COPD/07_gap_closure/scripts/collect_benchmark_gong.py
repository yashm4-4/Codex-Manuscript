#!/usr/bin/env python3
"""Retrieve publicly linked Gong et al. source evidence; never reads model scores.

Run from anywhere. Downloads are cached unchanged; metadata records URLs and
SHA-256. The public supplement URL was discovered on the publisher abstract page.
No credentials, access-control bypass, participant data, or inference is used.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "data/COPD-V2-BENCH/gong"
SOURCES = {
    "pubmed.html": "https://pubmed.ncbi.nlm.nih.gov/42745457/",
    "publisher.html": "https://academic.oup.com/ajrcmb/advance-article/doi/10.1093/ajrcmb/aanag179/8798774",
    "supplement.zip": "https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/ajrcmb/PAP/10.1093_ajrcmb_aanag179/1/aanag179_supplementary_data.zip",
    "europepmc.json": "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:42745457%20AND%20SRC:MED&format=json",
    "europepmc_core.json": "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:42745457%20AND%20SRC:MED&format=json&resultType=core",
    "crossref.json": "https://api.crossref.org/works/10.1093/ajrcmb/aanag179",
}
PUBLISHER_URL = SOURCES["publisher.html"]
SUPPLEMENT_URL = SOURCES["supplement.zip"]


def curate():
    """Re-express source facts under the locked protocol; no model inputs."""
    import pandas as pd
    from numbers_parser import Document
    from pypdf import PdfReader
    supp = DEST / "supplement"
    design = pd.read_csv(supp / "Table_S1_mpraDesignOligos.csv", keep_default_na=False)
    experimental = sorted(design.loc[design.type == "EXP", "ID"].unique())
    (DEST / "experimental_rsids.txt").write_text("\n".join(experimental) + "\n")
    numbers_rows = Document(supp / "Table_S3_mpraSigGwasResults.numbers").sheets[0].tables[0].rows(values_only=True)
    pd.DataFrame(numbers_rows[1:], columns=numbers_rows[0]).to_csv(DEST / "Table_S3_extracted.tsv", sep="\t", index=False)
    numbers_map = {r[0]: r[1] for r in numbers_rows[1:]}
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    tree = ET.fromstring(zipfile.ZipFile(supp / "MPRA_Sul_File_Final_CLEAN_July2026.docx").read("word/document.xml"))
    paragraphs = ["".join(t.text or "" for t in p.findall(".//w:t", ns)) for p in tree.findall(".//w:p", ns)]
    (DEST / "supplement_methods_extracted.txt").write_text("\n".join(f"{i}: {p}" for i, p in enumerate(paragraphs)) + "\n")
    pages = PdfReader(supp / "Supplementary_figures_V3_July2026.pdf").pages
    (DEST / "supplement_figures_extracted.txt").write_text("\n".join(f"PAGE {i + 1}\n{p.extract_text()}" for i, p in enumerate(pages)))
    rule = "Positive only when QuASAR, mpralm_sum and DESeq2 all have BH-adjusted P<0.05 (Supplementary Methods, MPRA data analysis); complete reported nonsignificant contrasts are null under this consensus rule. Missing any method or absent row is unavailable, never null."
    assays, contexts, sequence_rows, cell_audits = [], [], [], []
    resultbook = pd.ExcelFile(supp / "Table_S2_MPRA_RESULT_three_cell_lines.xls")
    summary = {}
    construct_qc = []
    complement = str.maketrans("ACGT", "TGCA")
    prefix = "ACTGGCCGCTTCACTG"
    linker = "GGTACCTCTAGA"
    for rsid, group in design.groupby("ID", sort=True):
        typed = {}
        for orientation in ["FOR", "REV"]:
            for role in ["ref", "alt"]:
                selected = group[(group.orientation == orientation) & group.oligoID.str.contains(f"_{role}_")]
                if selected.allele.nunique() != 1:
                    assert rsid == "rs76419734"
                    group.to_csv(DEST / "rs76419734_design_discrepancy.tsv", sep="\t", index=False)
                first = selected.iloc[0]
                sequence = first.sequence
                assert sequence.startswith(prefix) and linker in sequence
                insert = sequence[len(prefix):sequence.rfind(linker)]
                typed[(orientation, role)] = dict(allele=first.allele, construct_sequence=sequence,
                                               insert_sequence=insert, oligo_id=first.oligoID,
                                               csv_row=int(first.name) + 2)
                sequence_rows.append(dict(rsid=rsid, panel_role=first.type, orientation=orientation,
                                         allele_role=role, **typed[(orientation, role)]))
        assert typed[("FOR", "ref")]["allele"] == typed[("REV", "ref")]["allele"]
        assert typed[("FOR", "alt")]["allele"] == typed[("REV", "alt")]["allele"]
        role_checks = {}
        for role in ["ref", "alt"]:
            forward = typed[("FOR", role)]["insert_sequence"]
            reverse = typed[("REV", role)]["insert_sequence"]
            expected = forward.translate(complement)[::-1]
            exact = reverse == expected
            end_shift = (len(reverse) == len(expected) and
                         (reverse[1:] == expected[:-1] or reverse[:-1] == expected[1:]))
            role_checks[role] = "exact_reverse_complement" if exact else "one_base_end_shift" if end_shift else "internal_insert_mismatch"
            construct_qc.append(dict(rsid=rsid, panel_role=group.type.iloc[0], allele_role=role,
                                     forward_insert=forward, reverse_insert=reverse, expected_reverse_insert=expected,
                                     status=role_checks[role]))
        summary[rsid] = dict(panel_role=group.type.iloc[0], oligo_count=len(group), alleles=typed,
                             construct_checks=role_checks)
    pd.DataFrame(sequence_rows).to_csv(DEST / "allele_constructs.tsv", sep="\t", index=False)
    pd.DataFrame(construct_qc).to_csv(DEST / "construct_sequence_qc.tsv", sep="\t", index=False)
    indel_rows = [r for r in sequence_rows if r["panel_role"] == "EXP" and (len(r["allele"]) != 1 or r["allele"] == "-")]
    pd.DataFrame(indel_rows).to_csv(DEST / "indel_allele_constructs.tsv", sep="\t", index=False)
    for sheet in resultbook.sheet_names:
        cell = sheet.replace(" CELLS", "")
        results = pd.read_excel(resultbook, sheet_name=sheet, header=1)
        results.to_csv(DEST / f"Table_S2_{cell}_extracted.tsv", sep="\t", index=False)
        indexed = {r.simN: (i + 3, r) for i, r in results.iterrows()}
        counts = {"positive": 0, "null": 0, "unavailable": 0}
        for rsid, info in summary.items():
            f = info["alleles"]
            for orientation in ["FOR", "REV"]:
                key = f"{rsid}:{orientation}"
                match = indexed.get(key)
                stats = {m: (float(match[1][m]) if match and pd.notna(match[1][m]) else None)
                         for m in ["QuASAR", "mpralm_sum", "DESeq2"]}
                missing = [m for m, p in stats.items() if p is None]
                state = "unavailable" if missing else ("positive" if all(p < .05 for p in stats.values()) else "null")
                qnote = ("Source result row absent; reason not documented. S1 duplicates this variant's design oligos (160 versus usual80) with swapped ref/alt roles and differing flanks; first source representation retained here and all160rows preserved separately. C/Tcontrast remains, no outcome inferred." if not match else
                         "Reported row has missing method statistic(s): " + ",".join(missing) if missing else
                         "All three reported adjusted P values present; barcode/replicate QC was performed, but per-variant exclusion ledger is not supplied.")
                marker = numbers_map.get(rsid, "")
                locus = "unresolved_source_locus"
                if marker.startswith("chr3:"): locus = "EEFSEC"
                elif marker.startswith("chr4:"): locus = "GSTCD"
                elif marker.startswith("chr6:"): locus = "ADGRG6"
                elif marker.startswith("chr5:157"): locus = "ADAM19"
                elif marker.startswith("chr5:148"): locus = "HTR4"
                if rsid == "rs57658727": locus = "unresolved_source_locus"
                row = dict(assay_id=f"GONG2026_MPRA_{cell}_{orientation}_{rsid}", study_id="GONG2026",
                           source_id="GONG2026_S2", locus=locus, rsid=rsid, source_build="GRCh38",
                           source_chrom="", source_pos="", tested_allele1=f[("FOR", "ref")]["allele"],
                           tested_allele2=f[("FOR", "alt")]["allele"], reported_strand="+",
                           genomic_strand_evidence="Source S1 FOR allele-resolved constructs; build38 design stated in supplementary fine-mapping Methods; FOR/REV construct orientation does not complement the source allele labels. Requires central reference validation.",
                           assay_class="MPRA_allele_effect", cell_context=cell, assay_orientation=orientation,
                           experimental_state=state, reported_effect_direction="unavailable", effect_allele="",
                           p_value="", fdr=stats, source_rule=rule,
                           mechanism_in_model_scope="yes", mechanism_scope_rationale="Allele-specific transcriptional reporter is broadly regulatory-sequence scope; cell line and short episomal construct differ from frozen bulk-lung sequence-model context.",
                           source_locator=f"Table S2; sheet={sheet}; Excel row={match[0] if match else 'ABSENT'}; simN={key}; S1 FOR ref/alt CSV rows={f[('FOR','ref')]['csv_row']},{f[('FOR','alt')]['csv_row']}",
                           source_url=PUBLISHER_URL, source_result_present=bool(match), source_stats_complete=not missing,
                           source_panel_role=info["panel_role"], source_gwas_marker=marker,
                           source_gwas_marker_build="not separately declared in Table S3; not used as assay coordinate",
                           allele1_forward_construct=f[("FOR", "ref")]["construct_sequence"],
                           allele2_forward_construct=f[("FOR", "alt")]["construct_sequence"],
                           allele1_forward_insert=f[("FOR", "ref")]["insert_sequence"],
                           allele2_forward_insert=f[("FOR", "alt")]["insert_sequence"],
                           notes=qnote + " Source reference/alternate design labels are not assumed canonical GRCh38 REF/ALT. S2 contains no allele-effect estimate or direction; Table S3 effectallele is GWAS, never reporter direction.")
                mismatch = orientation == "REV" and "internal_insert_mismatch" in info["construct_checks"].values()
                row["construct_sequence_checks"] = info["construct_checks"]
                row["construct_identity_status"] = "unresolved_reverse_insert_mismatch" if mismatch else "allele_resolved_constructs_require_reference_validation"
                row["exact_identity_eligible"] = not mismatch and rsid != "rs76419734"
                row["construct_exact_allele_identity_valid"] = row["exact_identity_eligible"]
                row["harmonization_exclusion_reason"] = ("Deposited reverse ALT insert is not reverse complement of forward ALT while REF is; do not join this assay to forward-label allele identity without independent contrast reconciliation." if mismatch else
                                                         "Duplicated source design has swapped ref/alt roles and different flanks; no assay results reported." if rsid == "rs76419734" else "")
                if info["panel_role"] == "EXP":
                    assays.append(row)
                    counts[state] += 1
                else:
                    row.update(context_id=row["assay_id"], exclusion_reason="S1 CON: preselected positive-control design, outside1,120variant COPD experimental denominator; preserve observed result, do not treat as unselected benchmark null/positive.")
                    contexts.append(row)
        cell_audits.append(dict(cell_context=cell, designed_variants=1120, designed_orientation_contrasts=2240,
                               reported_variants=1119, reported_orientation_contrasts=2238,
                               positive_orientation_contrasts=counts["positive"], null_orientation_contrasts=counts["null"],
                               unavailable_orientation_contrasts=counts["unavailable"],
                               missing_result_rows=2, missing_method_statistics=counts["unavailable"]-2,
                               complete_denominator_valid=False, reason="rs76419734 has no results and no exclusion reason; barcode/replicateQC documented but variant-levelQC ledger absent; retain missing statistics as unavailable."))
    # Region-targeting experiments are contextual, never allele-specific labels.
    regional = [("rs7684442", "GSTCD", "TBCK", .00022, "positive"),
                ("rs7684442", "GSTCD", "PPA2", .00339, "positive"),
                ("rs7684442", "GSTCD", "GSTCD", .08222, "null"),
                ("rs35421223", "EEFSEC", "RAB7A", .00484, "positive"),
                ("rs35421223", "EEFSEC", "RUVBL1", .09719, "null")]
    for rsid, locus, gene, p, state in regional:
        contexts.append(dict(context_id=f"GONG2026_CRISPRI_16HBE_{rsid}_{gene}", study_id="GONG2026", source_id="GONG2026_METHODS", rsid=rsid, locus=locus,
                             assay_class="CRISPRi_region", cell_context="16HBE", experimental_state=state,
                             target_gene=gene, p_value=p, source_rule="Table S7: Stouffer combined P<0.01; two gRNAs and three replicates each.",
                             mechanism_in_model_scope="unevaluable", mechanism_scope_rationale="Region-targeting knockdown is not an allele contrast.",
                             source_locator="Supplementary methods Table S7 and Figure S6", source_url=PUBLISHER_URL,
                             notes="Table S7 caption claims significant pairs listed, but GSTCD=.08222 and RUVBL1=.09719 exceed its stated .01 threshold. Abstract describes RUVBL1 regulation in16HBE; this conflicts with S7 combined test and is not resolved by changing the threshold."))
    for gene in ["RAB7A", "RUVBL1"]:
        contexts.append(dict(context_id=f"GONG2026_REGION_NHBE_rs35421223_{gene}", study_id="GONG2026", source_id="GONG2026_ART", rsid="rs35421223", locus="EEFSEC",
                             assay_class="CRISPRi_region", cell_context="NHBE", experimental_state="positive", target_gene=gene,
                             source_rule="Published abstract reports regional regulation; quantitative transcript statistics unavailable in accessible abstract.",
                             mechanism_in_model_scope="unevaluable", mechanism_scope_rationale="Regional perturbation, not an allele-isolating edit.",
                             source_locator="Abstract; supplementary Figure S7 protein follow-up", source_url=PUBLISHER_URL,
                             notes="Publication calls this CRISPRi, but supplementary NHBE delivery method describes HiFiCas9 nuclease RNP rather than dCas9repressor. Preserve authors' label with unresolved mechanistic inconsistency; never nucleotide validation."))
    for cell in ["HUVEC", "IMR90"]:
        for rsid, locus in [("rs35421223", "EEFSEC"), ("rs2955083", "EEFSEC"), ("rs6439124", "EEFSEC"), ("rs7684442", "GSTCD")]:
            contexts.append(dict(context_id=f"GONG2026_HIC_{cell}_{rsid}", study_id="GONG2026", source_id="GONG2026_FIGURES", rsid=rsid, locus=locus,
                                 assay_class="chromatin_contact", cell_context=cell, experimental_state="unevaluable",
                                 mechanism_in_model_scope="unevaluable", mechanism_scope_rationale="Public5kb Hi-C regional context does not isolate tested alleles.",
                                 source_locator="Figure S5; GSE63525; hg19", source_url=PUBLISHER_URL, notes="No allele-specific positive assigned from visual contact matrix."))
    contexts.append(dict(context_id="GONG2026_REPORTER_RESULTS_ACCESS", study_id="GONG2026", source_id="GONG2026_ART", rsid="", locus="multiple",
                         assay_class="conventional_reporter", cell_context="16HBE", experimental_state="unavailable",
                         mechanism_in_model_scope="yes", mechanism_scope_rationale="Methods describe150bp dual-luciferase ref/alt constructs in bothorientations; allele-specific outcome rows not publicly recoverable from accessed abstract/supplement.",
                         source_locator="Supplementary Methods Reporter assay; main article subscriber-only", source_url=PUBLISHER_URL,
                         notes="Methods state>=1.5-fold allele difference; no fabricated tested roster or positive labels from MPRA hit membership."))
    contexts.append(dict(context_id="GONG2026_QBIC_PREDICTIONS", study_id="GONG2026", source_id="GONG2026_S4", rsid="", locus="multiple",
                         assay_class="other_or_unresolved", cell_context="in_silico", experimental_state="unevaluable", mechanism_in_model_scope="unevaluable",
                         source_locator="Table S4", source_url=PUBLISHER_URL, notes="QBiC-Pred TF-binding predictions are not experimentally assayed allele-specific TFbinding outcomes."))
    source_files = [("ART", "Published accepted manuscript abstract preserved by Europe PMC", "europepmc_core.json"),
                    ("S1", "Complete design oligonucleotide table", "supplement/Table_S1_mpraDesignOligos.csv"),
                    ("S2", "All reported per-cell/orientation MPRA results", "supplement/Table_S2_MPRA_RESULT_three_cell_lines.xls"),
                    ("S3", "Selected significant variant GWAS results", "supplement/Table_S3_mpraSigGwasResults.numbers"),
                    ("S4", "Predicted TF-binding changes, not experimental labels", "supplement/Table_S4_TFs_analysis.csv"),
                    ("METHODS", "Supplementary Methods, Tables S5-S7 and figure legends", "supplement/MPRA_Sul_File_Final_CLEAN_July2026.docx"),
                    ("FIGURES", "Supplementary Figures S1-S7", "supplement/Supplementary_figures_V3_July2026.pdf")]
    sources = [dict(source_id=f"GONG2026_{key}", study_id="GONG2026", title=title, doi="10.1093/ajrcmb/aanag179", pmid="42745457", year=2026,
                    source_url=PUBLISHER_URL if key == "ART" else SUPPLEMENT_URL, local_path=str((DEST / path).relative_to(ROOT)),
                    sha256=hashlib.sha256((DEST / path).read_bytes()).hexdigest(),
                    access_status="abstract accessible via browser; terminal publisher retrieval403; main text subscriber-only" if key == "ART" else "public supplement retrieved",
                    ascertainment="COPD GWAS; five loci from Hobbs2017/Sakornsakolpat2019") for key, title, path in source_files]
    retrievals = [json.loads(p.read_text()) for p in sorted(DEST.glob("*.retrieval.json"))]
    queries = ['"42745457" COPD Gong', '"10.1093/ajrcmb/aanag179" supplement', '"Interrogation of functional variants" Gong preprint',
               '"rs35421223" "MPRA"', '"aanag179" pdf', '"Gong" "MPRA" "COPD" "rs35421223"',
               '"Interrogation of functional variants in COPD" data', '"rs35421223" "reporter"']
    searches = [dict(search_id=f"GONG2026_SEARCH_{i+1:02}", query=q, searched_utc_date="2026-10-05", engine="web search",
                     outcome="Primary publisher/PubMed article located; public supplement recovered via publisherlink. No public mainarticle or complete reporter outcomes located.") for i,q in enumerate(queries)]
    audit = dict(study_id="GONG2026", designed_unique_variants=1120, designed_experimental_variants=1120, designed_controls=34, designed_oligos=92400,
                 designed_experimental_oligos=89680, repeated_design_id="rs76419734", abstract_significant_variants=25,
                 observed_unique_positive_experimental_variants=len({r["rsid"] for r in assays if r["experimental_state"] == "positive"}),
                 complete_denominator_valid=False, complete_assayed_denominator_valid=False, per_cell=cell_audits, direction_available=False,
                 recoverable_complete_case_rows_by_context={c["cell_context"]:2240-c["unavailable_orientation_contrasts"] for c in cell_audits},
                 significance_rule=rule, source_build_evidence="Methods paragraph27: designs removed variants not liftable tobuild38; TableS1 allele constructs lackexplicitcoordinates.",
                 qc_details="Low-correlation DNA replicates excluded; mean DNA barcode counts <=10 in16HBE/MRC5 and <=3 inHUVEC excluded (Word Symbol F0A3 denotes <=). Per-variant QC/exclusion ledger not supplied.",
                 construct_internal_mismatch_ids=sorted({r["rsid"] for r in construct_qc if r["status"] == "internal_insert_mismatch"}),
                 construct_end_shift_ids=sorted({r["rsid"] for r in construct_qc if r["status"] == "one_base_end_shift"}),
                 raw_design_sequence_lengths=dict(design.sequence.str.len().value_counts().sort_index().astype(int).items()),
                 warning="Supplement calls144bp oligos but deposited fullconstruct lengths vary123-149; preserved exactly. TableS3 has24GWASrows/23unique IDs, repeatedrs6439124, and is not full MPRAdenominator or reporterdirection table.",
                 null_definition="Tested orientation/cellcontrast with all three reported methodstatistics present but consensusrule not met; not proof ofbiological inactivity.",
                 omitted_from_metrics_reason="IncompleteQC/resultdenominator; preserve recoverable nulls and unavailable states; no fullpanel sensitivity/specificity/AUROC claim authorized.")
    audit["locus_source_regions"] = [
        dict(locus=locus, source_build="GRCh37/hg19", chrom=chrom, start=start, end=end,
             designed_variant_count=n, source_locator="Supplementary Figure S1A counts and S1B hg19_MPRA_SNPs track display bounds; visually inspected", source_id="GONG2026_FIGURES")
        for locus,chrom,start,end,n in [
            ("GSTCD","4",106280877,106899567,492),
            ("EEFSEC","3",127795922,128133249,224),
            ("ADAM19","5",156731272,157143218,204),
            ("ADGRG6","6",142351024,142912277,108),
            ("HTR4","5",147700515,148074785,92),
        ]
    ]
    payload = dict(sources=sources, assays=assays, contexts=contexts, retrievals=retrievals, searches=searches, denominator_audit=audit)
    checks = {
        "experimental_design_ids": len(experimental) == 1120,
        "experimental_rows_complete_design_cross_product": len(assays) == 6720,
        "unique_assay_ids": len({r["assay_id"] for r in assays}) == len(assays),
        "unique_positive_variants_match_abstract": audit["observed_unique_positive_experimental_variants"] == 25,
        "positive_rows": sum(r["experimental_state"] == "positive" for r in assays) == 29,
        "explicit_complete_result_nulls": sum(r["experimental_state"] == "null" for r in assays) == 6671,
        "unavailable_not_null": sum(r["experimental_state"] == "unavailable" for r in assays) == 20,
        "missing_statistics_never_null": all(r["experimental_state"] == "unavailable" for r in assays if not r["source_stats_complete"]),
        "controls_separated": len([r for r in contexts if r.get("source_panel_role") == "CON"]) == 204,
        "all_experimental_rows_have_source_alleles": all(r["tested_allele1"] and r["tested_allele2"] for r in assays),
        "all_experimental_rows_have_source_locator": all(r["source_locator"] for r in assays),
        "no_inferred_reporter_direction": all(not r["effect_allele"] for r in assays),
        "internal_reverse_insert_mismatches_retained": len(audit["construct_internal_mismatch_ids"]) == 14,
        "unreconciled_construct_rows_blocked": sum(not r["construct_exact_allele_identity_valid"] for r in assays) == 48,
        "full_denominator_gate_closed": audit["complete_assayed_denominator_valid"] is False,
        "source_hashes_recorded": all(len(s["sha256"]) == 64 for s in sources),
    }
    assert all(checks.values()), checks
    (DEST / "curation_qc.json").write_text(json.dumps({"checks":checks, "all_pass":all(checks.values()), "model_inputs_read":False}, indent=2) + "\n")
    destination = ROOT / "data/COPD-V2-BENCH/gong_evidence.json"
    destination.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"assays":len(assays), "contexts":len(contexts), "denominator_audit":audit}, indent=2))


def fetch(name, url):
    path = DEST / name
    meta = path.with_name(path.name + ".retrieval.json")
    if meta.exists():
        return json.loads(meta.read_text())
    record = {"filename": name, "url": url,
              "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat()}
    req = urllib.request.Request(url, headers={"User-Agent": "COPD-V2-BENCH public-evidence-audit/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=60) as response:
            body = response.read()
            record.update(status=response.status, resolved_url=response.url,
                          content_type=response.headers.get("Content-Type", ""))
    except urllib.error.HTTPError as exc:
        body = exc.read()
        record.update(status=exc.code, resolved_url=exc.url,
                      content_type=exc.headers.get("Content-Type", ""), error=str(exc))
    except Exception as exc:
        body = b""
        record.update(status="error", error=str(exc))
    path.write_bytes(body)
    record.update(bytes=len(body), sha256=hashlib.sha256(body).hexdigest(),
                  local_path=str(path.relative_to(ROOT)))
    meta.write_text(json.dumps(record, indent=2) + "\n")
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--supplement-url", help="Publisher-provided public download URL if canonical URL requires its signed query")
    parser.add_argument("--roster-only", action="store_true")
    parser.add_argument("--curate", action="store_true", help="Parse cached sources only; requires pandas,xlrd,numbers-parser,pypdf")
    args = parser.parse_args()
    DEST.mkdir(parents=True, exist_ok=True)
    if args.curate:
        curate()
        return
    if args.roster_only:
        with (DEST / "supplement/Table_S1_mpraDesignOligos.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        ids = sorted({r["ID"] for r in rows if r["type"] == "EXP"})
        (DEST / "experimental_rsids.txt").write_text("\n".join(ids) + "\n")
        print(len(ids))
        return
    jobs = dict(SOURCES)
    if args.supplement_url:
        jobs["supplement_signed.zip"] = args.supplement_url
    rows = [fetch(name, url) for name, url in jobs.items()]
    for cached in sorted(DEST.glob("*.retrieval.json")):
        record = json.loads(cached.read_text())
        if record["filename"] not in jobs:
            rows.append(record)
    for row in rows:
        path = ROOT / row["local_path"]
        if path.suffix == ".zip" and zipfile.is_zipfile(path):
            with zipfile.ZipFile(path) as archive:
                members = archive.infolist()
                for member in members:
                    target = (DEST / "supplement" / member.filename).resolve()
                    if not target.is_relative_to((DEST / "supplement").resolve()):
                        raise ValueError("Unsafe archive member")
                    archive.extract(member, DEST / "supplement")
                row["archive_members"] = [m.filename for m in members]
    (DEST / "retrieval_manifest.json").write_text(json.dumps(rows, indent=2) + "\n")
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
