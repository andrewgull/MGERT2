# Project requirements

This document outlines the core functional and non-functional requirements for the MGERT2 project.

---

## Functional Requirements

| ID | Type | Description | Status | Acceptance criteria | Tool |
| :--- | :--- | :--- | :--- | :--- | :--- |
| R1.1 | Functional | BuildDabase | DONE | run BuildDatabase script to build a database from the input genome | RepeatModeler |
|R1.2 | Functional | de novo TE search and classification | DONE | Given a set of genomic sequences, this pipeline rule shall identify and classify mobile genetic elements (MGEs) present in the sequences. | [RepeatModeler2](https://github.com/Dfam-consortium/RepeatModeler) |
|R2 | Functional | Collect TE seqs of interest from R1 output | DONE | Given the output from R1 (a fasta file), this rule shall extract those consenus sequences specified by user. | Custom script |
|R3 | Functional | Annotate TE seqs from R2 output | DONE | Given the output from R2 (a fasta file), this rule shall annotate those consensus sequences specified by user. | [RepeatMasker](http://www.repeatmasker.org/) |
|R4 | Functional | Cut out TE regions from genome | DONE | Given the output from R3 (a GFF or TSV file), this rule shall extract those TE regions specified by user from the genome. | [BEDTools](https://bedtools.readthedocs.io/en/latest/) |
|R5 | Functional | ORF search | DONE | Given the output from R4 (a fasta file), this rule hall identify open reading frames within the sequences. | [ORFfinder](https://www.ncbi.nlm.nih.gov/orffinder/) |
|R6 | Functional | Domain search | DONE | Given the output from R5 (a fasta file), this rule shall identify conserved protein domains within the sequences. | [RPS-BLAST](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/rpsblast/) |
| R7 | Functional | Filtering | DONE | Given the output from R4, R5, R6, this rule shall filter out the TE seqs without ORFs, or CDs inside them | Custom script |
| R8 | Functional | Flanking regions | WONTDO | Given the output from R7, this rule shall add flanking regions to the input sequences | Custom script |
| R9 | Functional | Build a phylogenetic tree | DONE | Given the ORF with CD (RT for example), this rule shall build a phylogenetic tree | [IQ-Tree](https://iqtree.github.io/) |
---

This functionality has been implemented. The next steps can be found under the `Issues` tab.
