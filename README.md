# RNA Binding Protien Motif Tools (for rMATS)

RNA-binding proteins (RBPs) are proteins that bind to RNA molecules and are involved in many aspects of RNA processing, including alternative splicing. RBPs bind to RNA sequence motifs (about 3–5 bases) via their RNA recognition motif (RRM). Many RRMs and their binding sites have been experimentally characterized, and RBP-specific consensus motif sequences have been probabilistically determined (see [ESEfinder](https://esefinder.ahc.umn.edu/tools/ESE2/) and [RBPmap](http://rbpmap.technion.ac.il)).

The tools presented here are functions in R designed to provide motif scores based on nucleotide sequences and known RBP position probability matrices (PPM). **These tools are intended for downstream analysis of outputs from [rMATS](https://rmats.sourceforge.io)** (a computational tool for detecting differential alternative splicing events in RNA-Seq data), using a data frame input structured as shown with the necessary columns below:

| geneSymbol  | chr | strand | exonStart_0base | exonEnd   |
|-------|--------|-----|-----------------|-----------|
| SLK | 10     | +  | 81110914        | 81110955  |
| OPA1 | 3      | +  | 193626092       | 193626202 |
| LDHB | 12     | -  | 21657751        | 21657835  |

The tools are listed below:

#### **`extract_sequence`** - *Get Sequences from Genomic Coordinates*

This tool leverages the Bioconductor packages [BSgenome](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) and its associated genome (e.g. BSgenome.Hsapiens.UCSC.hg38) to extract nucleotide sequences based on genomic coordinates — exon start (exonStart_0base) to exon end (exonEnd). For entries where the strand is negative (-), representing a reverse orientation compared to the standard convention, the extracted sequence is reverse complemented. A new data frame will look like this:
  
| geneSymbol  | chr | strand | exonStart_0base | exonEnd   | nucleotide_sequence                                 |
|-------|--------|-----|-----------------|-----------|-----------------------------------------------------|
| SLK | 10      | +  | 81110914        | 81110955  | CCGCCGATGTGGAAGTGGCCAGATTCTGAGCCGCCTGACTAGA         |
| OPA1 | 3      | +  | 193626092       | 193626202 | GGTCTGCTTGGTGAGCTCATTCTCTTACAACAACAAATTCAAGAGCATGAAGAGGAAGCGCGCAGAGCCGCTGGCCAATATAGCACGAGCTATGCCCAACAGAAGCGCAAG           |
| LDHB | 12     | -  | 21657751        | 21657835  | AGCCGGCGCCGGAGGAGACGCACGCAGCTGACTTTGTCTTCTCCGCACGACTGTTACAGAGGTCTCCAGAGCCTTCTCTCTCCTG             |

#### **`extract_sequence_with_flanks`** - *Get Sequences and thier Flanking Regions from Genomic Coordinates*

This function is similar to `extract_sequence` except that it extracts both sequence defined by exonStart_0base to exonEnd and a desired number of flanking seqeunces (such as 50 +/- the exon region). For easy visualizaiton, the flanking (intronic) regions are in lowercase.

| geneSymbol  | chr | strand | exonStart_0base | exonEnd   | nucleotide_sequence                                 
|-------|--------|-----|-----------------|-----------|--------------------------------------------------------------------|
| SLK | 10      | +  | 81110914        | 81110955  | agtcctcagaccccatgctgcctccaactgagccttgtgtttccttgcagCGCCGATGTGGAAGTGGCCAGATTCTGAGCCGCCTGACTAGAgttagtaagttgcctggcgttctcgtgcagtcactggcctctccagtggt   |
| OPA1 | 3      | +  | 193626092       | 193626202 | attattctcctccccaatttcctcttctcctcattgtgaactcgtggcagGGTCTGCTTGGTGAGCTCATTCTCTTACAACAACAAATTCAAGAGCATGAAGAGGAAGCGCGCAGAGCCGCTGGCCAATATAGCACGAGCTATGCCCAACAGAAGCGCAAGgtgatggatggtttaagggggctaccgatacattcacactaatcagccat                                                                                                                  |
| LDHB | 12      | -  | 21657751        | 21657835  | taagaggctgcggtggttgtggggccccgccccctcctccctccttgcagAGCCGGCGCCGGAGGAGACGCACGCAGCTGACTTTGTCTTCTCCGCACGACTGTTACAGAGGTCTCCAGAGCCTTCTCTCTCCTGgtaggtttcggctcaggaccctgaatcctggcccacaggcaagcctgatg                                                                                                                                            |

### **Motif Scoring & Plots** 
#### **`find_rbp_binding_and_generate_table`** - *finds RBP binding sites, plots PWM motifs, and exports plots as a PDF to subfolder*

This is a wrapper function that contains the following:

- **`generate_pwm_from_ppm`** - *calculates an RBP PWM from input PPM*
- **`calculate_max_score`** - finds the maximum possible score for a PWM*
- **`score_motif_with_positions_pwm`** - *scores motifs using PWM or PPM input and returns the positions of high scores and their ranges*

Here, we can generate RBP sequence motif scores based on either the position probability matrix (PPM) or the position weighted matrix (PWM) of a characterized RBP. Such matrix information can be found and downloaded from certain online databases, such as [RBPmap](http://rbpmap.technion.ac.il/download.html). The below table represents a PPM of the *TIA1* gene. 

| Position |   1    |   2    |   3    |   4    |   5    |   6    |   7    |
|----------|--------|--------|--------|--------|--------|--------|--------|
| **A**    | 0.1502 | 0.0148 | 0.0148 | 0.0148 | 0.0148 | 0.0148 | 0.0843 |
| **C**    | 0.0620 | 0.0148 | 0.0148 | 0.0148 | 0.0148 | 0.2068 | 0.0843 |
| **G**    | 0.0620 | 0.0148 | 0.0148 | 0.0148 | 0.0148 | 0.2048 | 0.5665 |
| **T**    | 0.7258 | 0.9556 | 0.9556 | 0.9556 | 0.9556 | 0.5736 | 0.2649 |

PPMs represent the probabilities for each nucleotide at each position. These probabilities are derived from a matrix of raw frequencies (or counts) of each nucleotide letter at each position, known as a position count matrix (PCM). The PPM normalizes such counts (note that all positions sum to 1 in the table). A pseudocount of 1 is applied.

To correct for background nucleotide frequencies, position weight matrices (PWMs) can be used. A PWM presents the "expected" motif sequence scores by incorporating background nucleotide probabilities, assuming uniform distribution (i.e., $1/4$ for each nucleotide). To compute a PWM from a PPM, the `generate_pwm_from_ppm` function applies the following log-odds scoring formula:

### $S(N) = \log_2 \left( \frac{P(C_N)}{B_N} \right)$

Where $S(N)$ is the score at position $N$, $P(C_N)$ is the probability of the observed counts $C_N$ (which corresponds to the PPM), and $B_N$ represents the background nucleotide frequency at position $N$.

For a more detailed explanation of the mathematical principles behind sequence motif matrices, refer to [this vignette](https://bioconductor.org/packages/devel/bioc/vignettes/universalmotif/inst/doc/IntroductionToSequenceMotifs.pdf).

The function `score_motif_with_positions_pwm` slides a window of (i.e. scans) the motif's length across the sequence one position at a time and generates a score based on the matching nucleotides of the PWM at each window. Of course, we shouldn't accept all the motifs that are identified, so a threshold is applied. The thresholding method used here discards scores less than or equal to a percentage of the maximal possible score for a PWM (the max score is obtained by the `calculate_max_score` function.

The `find_rbp_binding_and_generate_table` function is the wrapper function that outputs results. (1) including a plot for each sequence, displaying the motif scores of each RBP along the sequence with reference to the position of the exon (see below). These plots are exported as PDFs to subfolder in the working directory.

<img width="458" alt="image" src="https://github.com/user-attachments/assets/31bfa738-4f54-444e-8a63-aa7abaf25085">

The function also generates a table of the RBP motif analyses, including the gene symbol, RBP, motif sequence, the motif's score, its relative start and end positions within the sequence, and it's genomic start and end coordinates. 

 | ID | geneSymbol | RBP  | motif_sequence | score  | start | end | genomic_start | genomic_end |
|----|------|------|----------------|--------|-------|-----|---------------|-------------|
| 1  | SLK  | REM47| AAATGA         | 4.52175| 56    | 61  | 104010871     | 104010876   |
| 2  | OPA1 | ILF2 | GCACGGT        | 5.47395| 46    | 53  | 193626137     | 193626144   |
| 2  | OPA1 | REM47| CAATTT         | 5.81979| 15    | 20  | 193626106     | 193626111   |
| 2  | OPA1 | REM47| AAATTC         | 5.507349| 85   | 90  | 193626176     | 193626181   |
| 2  | OPA1 | REM47| TAATCA         | 8.626369| 201  | 206 | 193626292     | 193626297   |
| 3  | LDHB | ILF2 | GCTCGCG        | 5.158135| 6    | 13  | 21657756      | 21657763    |
| 3  | LDHB | ILF2 | GTTGTCG        | 5.154187| 15   | 22  | 21657765      | 21657772    |
| 3  | LDHB | ILF2 | TTGTCGGG       | 5.864684| 17   | 24  | 21657767      | 21657774    |
| 3  | LDHB | ILF2 | CGGACGAG       | 5.948750| 60   | 67  | 21657810      | 21657817    |
| 3  | LDHB | REM47| GAATCC         | 4.500737| 158  | 163 | 21657908      | 21657913    |



## What's Next?

For these tools, I plan to streamline them further and adapt them to accept other rMATS input formats (i.e. IR, A3SS, A5SS, and MXE tables). For a future RBP tool, I aim to develop functions that plot motif enrichment of alternatively spliced events against their ΔPSI values. This will enable investigations into how different RBPs contribute to the differential motif enrichment of alternatively spliced events across two conditions.
  
