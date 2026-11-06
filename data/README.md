# Data

| File | Content |
|---|---|
| `accessions.tsv` | Table 2 of the Supplementary Material: one row per mtDNA control region haplotype (H1-H64) |
| `fetch_genbank.R` | Downloads all haplotype sequences with a GenBank accession to `sequences.fasta` (FASTA names are the haplotype IDs) |

## Columns of `accessions.tsv`

`haplotype`, `marker`, `genbank_accession`, `reference` (study that first reported the haplotype), `new_in_this_study` (H61-H64, accessions OK188832-OK188835), `total_n` (individuals, 400 in total), then one column per sampling site with the number of individuals carrying the haplotype.

## Sampling site codes

ES: Espírito Santo; RJN: northern Rio de Janeiro; RJS: southern Rio de Janeiro; SPN: northern São Paulo; SPC: central São Paulo; SPS: southern São Paulo; PR: Paraná; SCN: northern Santa Catarina; BAB: Babitonga bay; RS: Rio Grande do Sul; URAO: Atlantic Ocean of Uruguay; RP: Río de la Plata; SCL: San Clemente; SB: San Bernardo; NC: Necochea; CL: Claromecó; MH: Monte Hermoso; BB: Bahía Blanca; RN: Río Negro; UNK: unknown location.

## Notes

- Haplotypes H12 and H14-H17 are listed as "unsubmitted" in the paper, so `fetch_genbank.R` skips them.
- The BioGeoBEARS tip ranges (`biogeography/input/geography.data`) are built from this table once the biogeographic areas of the paper are defined, grouping the sampling sites into areas.
