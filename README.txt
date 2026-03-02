This pipeline was created to massively screen the upstream regions of genes to assess their energetic and structural properties (ESp) and identify with statistical significance the genes with ESp values ​​most deviating from the mean. The system was designed for the in silico prediction of putative multigene transcriptional units.

REQUIREMENTS:
A) A conda environment with all dependencies. Required packages can be found in the ESp.yml file.
B) A file with the genomic sequence (fasta-like) and its annotation file (e.g. .gff output of PROKKA, https://github.com/tseemann/prokka).
C) The three tab_n.csv files containing the ESp values ​​for different dinucleotide pairs
D) All five scripts that make up the pipeline to be launched in sequential order. 

PIPELINE OVERVIEW
step_01_parameter_extraction.py	-> Extracts upstream sequences and calculates average DNA ESp of each gene.
step_02_iterative_screening.py	-> Performs an iterative percentile-based screening to define general and biological outlier profiles.
step_03_anomaly_frequency.py		-> Aggregates the results of the screening and generates a trend analysis graph.
step_04_multi_scoring.py		-> Calculates anomaly scores (weighted, cluster-based, and Z-score) for a specific threshold.
step_05_combined_scoring.py		-> Computes a continuous combined score to rank gene candidates.

CONTRIBUTING
Please feel free to suggest improvements or report issues by opening a pull request.
