# RNA Binding Protien Motif Tools (for rMATS)

RNA-binding proteins (RBPs) are proteins that bind to RNA molecules and are involved in many aspects of RNA processing, including alternative splicing. RBPs bind to regulatory elements on RNA via their RNA recognition motif (RRM). Many RRMs have been experimentally characterized, and RBP-specific consensus motifs have been probabilistically deterined (see [ESEfinder](https://esefinder.ahc.umn.edu/tools/ESE2/) and [RBPmap](http://rbpmap.technion.ac.il)).

The tools presented here are functions in R designed to provide motif scores based on nucleotide sequences and known RBP position probability matrices (PPM). **These tools are intended for downstream analysis of outputs from [rMATS](https://rmats.sourceforge.io)** (a computational tool for detecting differential alternative splicing events in RNA-Seq data), using a data frame input structured as shown with the necessary columns below:

| gene  | chr | strand | exonStart_0base | exonEnd   |
|-------|--------|-----|-----------------|-----------|
| gene1 | 17     | +  | 81110914        | 81110955  |
| gene2 | 3      | +  | 193626092       | 193626202 |
| gene3 | 12     | -  | 21657751        | 21657835  |

The tools are listed below:

#### **`extract_sequence`** - *Get Sequences from Genomic Coordinates*

This tool leverages the Bioconductor packages [BSgenome](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) and its associated genome (e.g. BSgenome.Hsapiens.UCSC.hg38) to extract nucleotide sequences based on genomic coordinates — exon start (exonStart_0base) to exon end (exonEnd). For entries where the strand is negative (-), representing a reverse orientation compared to the standard convention, the extracted sequence is reverse complemented. A new data frame will look like this:
  
| gene  | chr | strand | exonStart_0base | exonEnd   | nucleotide_sequence                                 |
|-------|--------|-----|-----------------|-----------|-----------------------------------------------------|
| gene1 | 17     | +  | 81110914        | 81110955  | CCGCCGATGTGGAAGTGGCCAGATTCTGAGCCGCCTGACTAGA         |
| gene2 | 3      | +  | 193626092       | 193626202 | GGTCTGCTTGGTGAGCTCATTCTCTTACAACAACAAATTCAAGAGCATGAAGAGGAAGCGCGCAGAGCCGCTGGCCAATATAGCACGAGCTATGCCCAACAGAAGCGCAAG           |
| gene3 | 12     | -  | 21657751        | 21657835  | AGCCGGCGCCGGAGGAGACGCACGCAGCTGACTTTGTCTTCTCCGCACGACTGTTACAGAGGTCTCCAGAGCCTTCTCTCTCCTG             |

#### **`extract_sequence_with_flanks`** - *Get Sequences and thier Flanking Regions from Genomic Coordinates*

This function is similar to `extract_sequence` except that it extracts both sequence defined by exonStart_0base to exonEnd and a desired number of flanking seqeunces (such as 50 +/- the exon region). For easy visualizaiton, the flanking (intronic) regions are in lowercase.

| gene  | chr | strand | exonStart_0base | exonEnd   | nucleotide_sequence                                 
|-------|--------|-----|-----------------|-----------|--------------------------------------------------------------------|
| gene1 | 17      | +  | 81110914        | 81110955  | agtcctcagaccccatgctgcctccaactgagccttgtgtttccttgcagCGCCGATGTGGAAGTGGCCAGATTCTGAGCCGCCTGACTAGAgttagtaagttgcctggcgttctcgtgcagtcactggcctctccagtggt   |
| gene2 | 3      | +  | 193626092       | 193626202 | attattctcctccccaatttcctcttctcctcattgtgaactcgtggcagGGTCTGCTTGGTGAGCTCATTCTCTTACAACAACAAATTCAAGAGCATGAAGAGGAAGCGCGCAGAGCCGCTGGCCAATATAGCACGAGCTATGCCCAACAGAAGCGCAAGgtgatggatggtttaagggggctaccgatacattcacactaatcagccat                                                                                                                  |
| gene3 | 12      | -  | 21657751        | 21657835  | taagaggctgcggtggttgtggggccccgccccctcctccctccttgcagAGCCGGCGCCGGAGGAGACGCACGCAGCTGACTTTGTCTTCTCCGCACGACTGTTACAGAGGTCTCCAGAGCCTTCTCTCTCCTGgtaggtttcggctcaggaccctgaatcctggcccacaggcaagcctgatg                                                                                                                                            |

#### **`score_motif`** - *Get a Sequence Motif Score for an RBP Based on its PPM*

Here, we can generate sequence motif scores based on either the position probability matrix (PPM) or the position weighted matrix (PWM) of a characterized RBP. Such matrix information can extracted from databases, such as [RBPmap](http://rbpmap.technion.ac.il/download.html). The below table represents a PPM of the *TIA1* gene. 

| Position |   1    |   2    |   3    |   4    |   5    |   6    |   7    |
|----------|--------|--------|--------|--------|--------|--------|--------|
| **A**    | 0.1502 | 0.0148 | 0.0148 | 0.0148 | 0.0148 | 0.0148 | 0.0843 |
| **C**    | 0.0620 | 0.0148 | 0.0148 | 0.0148 | 0.0148 | 0.2068 | 0.0843 |
| **G**    | 0.0620 | 0.0148 | 0.0148 | 0.0148 | 0.0148 | 0.2048 | 0.5665 |
| **T**    | 0.7258 | 0.9556 | 0.9556 | 0.9556 | 0.9556 | 0.5736 | 0.2649 |

PPMs represent the probabilities for each nucleotide at each position. These probabilities are derived from a matrix of frequencies (or counts) of each nucleotide at each position, known as a position count matrix (PCM). Note that in the PPM, all positions sum to 1.



S(N) = log2(P(C_N) / B_N)

 
