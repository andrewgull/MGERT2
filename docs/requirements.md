# Project requirements

This document outlines the core functional and non-functional requirements for the MGERT2 project.

---

## Functional Requirements

| ID | Type | Description | Status | Acceptance criteria | Tool |
| :--- | :--- | :--- | :--- | :--- | :--- |
| R1.1 | Functional | BuildDabase | ... | run BuildDatabase script to build a database from the input genome | RepeatModeler |
|R1.2 | Functional | de novo TE search and classification | ... | Given a set of genomic sequences, this pipeline rule shall identify and classify mobile genetic elements (MGEs) present in the sequences. | [RepeatModeler2](https://github.com/Dfam-consortium/RepeatModeler) |
|R2 | Funcitonal | Collect TE seqs of interest from R1 output | ... | Given the output from R1 (a fasta file), this rule shall exatract those consenus sequences specified by user. | Custom script |
|R3 | Functional | Annotate TE seqs from R2 output | ... | Given the output from R2 (a fasta file), this rule shall annotate those consenus sequences specified by user. | [RepeatMasker](http://www.repeatmasker.org/) |
|R4 | Functional | Cut out TE regions from genome | ... | Given the output from R3 (a GFF or TSV file), this rule shall extract those TE regions specified by user from the genome. | [BEDTools](https://bedtools.readthedocs.io/en/latest/) |
|R5 | Functional | ORF search | ... | Given the output from R4 (a fasta file), this rule hall identify open reading frames within the sequences. | [ORFfinder](https://www.ncbi.nlm.nih.gov/orffinder/) |
|R6 | Functional | Domain search | ... | Given the output from R5 (a fasta file), this rule shall identify conserved protein domains within the sequences. | [RPS-BLAST](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/rpsblast/) |
| R7 | Functional | Filtering | ... | Given the output from R4, R5, R6, this rule shall filter out the TE seqs without ORFs, or CDs inside them | Custom script |
| R8 | Fuctional | Flanking regions | ... | Given the output from R7, this rule shall add flanking regions to the input sequences | Custom script |
| R9 | Functional | ... | Build a phylogenetic tree | Given the ORF iwth CD (RT for example), this rule shall build a phylogenetic tree | [IQ-Tree](https://iqtree.github.io/) |
---

## More detailed description

**R1 de novo TE search and classification**:
  - tool: RepeatModeller2
  - input: genome assembly, .fasta.gz
  - output: consesnus sequences, .fasta
  - tasks:
    - [ ] make & test conda env
    - [ ] rule
    - [ ] test
