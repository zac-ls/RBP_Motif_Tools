# RNA Binding Protien Motif Tools (for rMATS)

RNA-binding proteins (RBPs) are proteins that bind to RNA molecules and are involved in many aspects of RNA processing, including alternative splicing. RBPs bind to regulatory elements on RNA via their RNA recognition motif (RRM). Many RRMs have been experimentally characterized, and RBP-specific consensus motifs have been probabilistically deterined (see [ESEfinder](https://esefinder.ahc.umn.edu/tools/ESE2/) and [RBPmap](http://rbpmap.technion.ac.il)).

The tools presented here are functions in R designed to provide motif scores based on nucleotide sequences and known RBP position probability matrices (PPM). **These tools are intended for downstream analysis of outputs from [rMATS](https://rmats.sourceforge.io)** (a computational tool for detecting differential alternative splicing events in RNA-Seq data), using a data frame input structured as shown with the necessary columns below:

| gene  | strand | chr | exonStart_0base | exonEnd   |
|-------|--------|-----|-----------------|-----------|
| gene1 | +      | 17  | 81110914        | 81110955  |
| gene2 | +      | 3   | 193626092       | 193626202 |
| gene3 | –      | 12  | 21657751        | 21657835  |

The tools are listed below:

#### **`extract_sequence`** - *Get Sequences from Genomic Coordinates*

This tool leverages the Bioconductor packages [BSgenome](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) and its associated genome (e.g. BSgenome.Hsapiens.UCSC.hg38) to extract nucleotide sequences based on genomic coordinates — exon start (exonStart_0base) to exon end (exonEnd). For entries where the strand is negative (-), representing a reverse orientation compared to the standard convention, the extracted sequence is reverse complemented. A new data frame will look like this:
  
| gene  | strand | chr | exonStart_0base | exonEnd   | nucleotide_sequence                                 |
|-------|--------|-----|-----------------|-----------|-----------------------------------------------------|
| gene1 | +      | 17  | 81110914        | 81110955  | CCGCCGATGTGGAAGTGGCCAGATTCTGAGCCGCCTGACTAGA         |
| gene2 | +      | 3   | 193626092       | 193626202 | GGTCTGCTTGGTGAGCTCATTCTCTTACAACAACAAATTCAAGAGCATGAAGAGGAAGCGCGCAGAGCCGCTGGCCAATATAGCACGAGCTATGCCCAACAGAAGCGCAAG           |
| gene3 | –      | 12  | 21657751        | 21657835  | AGCCGGCGCCGGAGGAGACGCACGCAGCTGACTTTGTCTTCTCCGCACGACTGTTACAGAGGTCTCCAGAGCCTTCTCTCTCCTG             |

#### **`extract_sequence_with_flanks`** - *Get Sequences and thier Flanking Regions from Genomic Coordinates*

This function is similar to `extract_sequence` except that it extracts both sequence defined by exonStart_0base to exonEnd and a desired number of flanking seqeunces (such as 50 +/- the exon region). For easy visualizaiton, the flanking (intronic) regions are in lowercase.

| gene  | strand | chr | exonStart_0base | exonEnd   | nucleotide_sequence                                 
|-------|--------|-----|-----------------|-----------|--------------------------------------------------------------------|
| gene1 | +      | 17  | 81110914        | 81110955  | agtcctcagaccccatgctgcctccaactgagccttgtgtttccttgcagCGCCGATGTGGAAGTGGCCAGATTCTGAGCCGCCTGACTAGAgttagtaagttgcctggcgttctcgtgcagtcactggcctctccagtggt   |
| gene2 | +      | 3   | 193626092       | 193626202 | attattctcctccccaatttcctcttctcctcattgtgaactcgtggcagGGTCTGCTTGGTGAGCTCATTCTCTTACAACAACAAATTCAAGAGCATGAAGAGGAAGCGCGCAGAGCCGCTGGCCAATATAGCACGAGCTATGCCCAACAGAAGCGCAAGgtgatggatggtttaagggggctaccgatacattcacactaatcagccat                                                                                                                  |
| gene3 | –      | 12  | 21657751        | 21657835  | taagaggctgcggtggttgtggggccccgccccctcctccctccttgcagAGCCGGCGCCGGAGGAGACGCACGCAGCTGACTTTGTCTTCTCCGCACGACTGTTACAGAGGTCTCCAGAGCCTTCTCTCTCCTGgtaggtttcggctcaggaccctgaatcctggcccacaggcaagcctgatg                                                                                                                                            |

#### **`score_motif`** - *Get a Sequence Motif Score for an RBP Based on its PPM*


| Position |   1  |   2   |   3   |   4  |
|----------|-------|-------|-------|-------|
| **A**   | 0.25  | 0.30  | 0.35  | 0.10  |
| **C**    | 0.20  | 0.25  | 0.40  | 0.15  |
| **G**    | 0.15  | 0.10  | 0.65  | 0.10  |
| **T**    | 0.40  | 0.20  | 0.30  | 0.10  |
| **A**    | 0.30  | 0.20  | 0.25  | 0.25  |


