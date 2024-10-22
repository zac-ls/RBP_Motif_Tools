# RNA Binding Protien Motif Tools

RNA-binding proteins (RBPs) are proteins that bind to RNA molecules and are involved in many aspects of RNA processing, including alternative splicing. RBPs bind to regulatory elements on RNA via their RNA recognition motif (RRM). Many RRMs have been experimentally characterized, and RBP-specific consensus motifs have been probabilistically deterined (see [ESEfinder](https://esefinder.ahc.umn.edu/tools/ESE2/) and [RBPmap](http://rbpmap.technion.ac.il)).

The tools presented here are functions in R designed to provide motif scores based on nucleotide sequences and known RBP position probability matrices (PPM). **These tools are intended for downstream analysis of outputs from [rMATS](https://rmats.sourceforge.io)** (a computational tool for detecting differential alternative splicing events in RNA-Seq data), using a data frame input structured as shown with the necessary columns below:

| gene  | strand | chr | exonStart_0base | exonEnd   |
|-------|--------|-----|-----------------|-----------|
| gene1 | +      | 17  | 81110914        | 81110955  |
| gene2 | +      | 3   | 193626092       | 193626202 |
| gene3 | –      | 12  | 21657751        | 21657835  |
