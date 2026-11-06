# Downloads every sequence listed in accessions.tsv from GenBank into sequences.fasta
# Run from the data/ folder. Requires the ape package.
library(ape)

acc <- read.delim("accessions.tsv", stringsAsFactors = FALSE)
acc <- acc[acc$genbank_accession != "unsubmitted", ]
seqs <- read.GenBank(acc$genbank_accession)
names(seqs) <- acc$haplotype[match(names(seqs), acc$genbank_accession)]
write.dna(seqs, "sequences.fasta", format = "fasta", nbcol = -1, colsep = "")
