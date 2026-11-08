# AFLA CFTR Nanopore

Analyses targeted CFTR sequencing from Nanopore R10.4.1 reads and writes a **one-page report per sample**: small variants (bcftools by
default, or Clair3), optional read-based phasing (WhatsHap), CFTR2 and ClinVar classes, HGVS and legacy names, poly-T/TG, exon deletions
and duplications, and coverage gaps. Runs in [EPI2ME Desktop](https://epi2me.nanoporetech.com/). **Education and research use only. Not for clinical use.**

## What you need
| Item | Where it comes from |
|---|---|
| EPI2ME Desktop with Docker | [epi2me.nanoporetech.com](https://epi2me.nanoporetech.com/) |
| Basecalled FASTQ, one folder or file per sample | MinKNOW barcode folders, or the AFLA Demux workflow |
| GRCh38 reference FASTA (with `chr` names) | For example the GATK/Broad hg38 analysis set |
| CFTR2 variant list (`.xlsx`) | Download from [cftr2.org](https://cftr2.org) after accepting its usage agreement. Not included here. |
| Ensembl VEP cache, version 115, GRCh38 (optional) | Ensembl; needed for annotation |
| ClinVar GRCh38 VCF with `.tbi` (optional) | NCBI |

## Good to know
- **Use the default bcftools caller.** Clair3 is available but experimental (see Testing status).
- Nanopore miscounts repeats, so **poly-T and poly-TG results are labelled lower confidence**; confirm 5T and TG counts another way.
- Include **at least 4 samples** from the run for exon deletions and duplications.
- Read-based phasing only joins variants carried by the same reads. Variants not joined are reported as cis/trans unknown.
- Primer bases are not removed in this workflow, and the report says so.
- To save time, put a minimap2 index named `<reference name>.map-ont.mmi` next to the reference; otherwise each run builds one.

Install in EPI2ME with **Import workflow** and this repository's address.

## Testing status
Tested on **one real Nanopore amplicon sample** with the default bcftools caller. Small variants, F508del and poly-T agreed with the
laboratory's results; the TG repeat count differed by one on one allele.

**Clair3 is experimental here: on that sample it missed F508del in a full run**, although 46% of the reads carry it. Run over that
amplicon alone, it called it only with low quality (5 to 22), and not at all with the older model. Clair3 v1.0.8 from ONT's EPI2ME image
behaved the same way across five R10.4.1 models (called at quality 5 to 12 by some, not called by others). Use bcftools unless you have validated
Clair3 on your own known samples. Exon copy number was only checked for false calls, by splitting one real sample into four; real
deletions on real Nanopore runs have not been tested. Not validated for any clinical purpose: validate it on your own known samples first.
