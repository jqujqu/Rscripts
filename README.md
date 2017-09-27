# Rscripts

## sraInfo.R

### Dependencies:  
- R package [`optparse`](https://cran.rstudio.com/web/packages/optparse/index.html)
- Bioconductor package [`SRAdb`](https://www.bioconductor.org/packages/release/bioc/html/SRAdb.html)
- (Optional--can be downloaded during first execution of `sraInfo.R`) Pre-downloaded and uncompressed  `SRAmetadb.sqlite` (see [`SRAdb` manual](https://www.bioconductor.org/packages/devel/bioc/vignettes/SRAdb/inst/doc/SRAdb.pdf))

### Usage
Get usage information
``` bash
Rscript sraInfo.R -h
```
Example usage
```
Rscript sraInfo.R -d /Users/jenny/SRAmetadb.sqlite  -o ${PWD}/new_study -l ChIP-Seq SRP057141 
```
### Output:
An `<Accession>.csv` file containing the following columns: 

Column name | Example values
------------|---------------
study	| SRP057141
submission 	| SRA258544
sample	| SRS908811
experiment| SRX994293
run	| SRR1977542
ftp	| ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX994/SRX994293/SRR1977542/SRR1977542.sra
experiment_alias	| human1_pachytene_H3K4me3
sample_alias	| human1_pachytene
scientific_name	| Homo sapiens
library_name	| human1_PS_K4
library_strategy|ChIP-Seq	
library_layout	|SINGLE - 
sample_attribute| isolate: BL2 \|\| age: missing \|\| biomaterial_provider: Dr. Sherman Silber, Infertility Center of St. Louis, St. Luke's Hospital, St. Louis, MO, USA \|\| sex: male \|\| tissue: testis \|\| cell_type: pachytene spermatocyte \|\| phenotype: fertile \|\| BioSampleModel: Human
