Fractionation-seq physically separates cellular components like the nucleus and cytoplasm. By sequencing RNA from these fractions, it quantitatively maps the subcellular location of transcripts. 
This reveals crucial insights into localized gene regulation and transport.

The Fractionation-seq used here is from the cytoplasmic and nuclear fractions of K562 cells. Data source: ENCODE Project Consortium. 
An integrated encyclopedia of DNA elements in the human genome. Nature 2012 Sep 6;489(7414):57-74. PMID: 22955616

The analysis script includes data processing steps such as downloading and preprocessing data, performing quality control (QC), mapping with Hisat, quantifying with HTSeq-count or featureCounts, 
and conducting Differential Gene Expression (DGE) analysis.
