# ADVICER - Analysis Dashboard for Virus-Induced CEll Response based on RNA-seq data

ADVICER is an R Shiny application developed to interactively visualize the differential gene expression analysis results of the [SFB1021][sfb] joint virus project (https://www.sfb1021.de/). In this project HuH7 cells were infected with different human pathogenic viruses for defined time periods. Total RNA was analyzed to identify common and unique transcriptional changes.
The application uses DESeq2 result tables to: 
* generate MAplots and Volcano plots for each table
* compare the differentially expressed genes (DEGs) between different time points of a virus 
* compare DEGs between different viruses 
* plot the expression of a specific gene over time 
* show SNPs and Indels for each virus over time

## Access to ADVICER

### https://advicer.computational.bio/

## Citation/Paper

Hofmann N, Bartkuhn M, Becker S, Biedenkopf N, Böttcher-Friebertshäuser E, Brinkrolf K, Dietzel E, Fehling SK, Goesmann A, Heindl MR, Hoffmann S, Karl N, Maisner A, Mostafa A, Kornecki L, Müller-Kräuter H, Müller-Ruttloff C, Nist A, Pleschka S, Sauerhering L, Stiewe T, Strecker T, Wilhelm J, Wuerth JD, Ziebuhr J, Weber F, Schmitz ML.2024.Distinct negative-sense RNA viruses induce a common set of transcripts encoding proteins forming an extensive network. J Virol98:e00935-24.https://doi.org/10.1128/jvi.00935-24

## Contact

Nina Hofmann<br />
advicer@computational.bio<br />

Bioinformatics and Systems Biology<br />
Justus Liebig University Giessen<br />
35392 Giessen<br />
Germany

[sfb]: https://www.sfb1021.de/ ""
