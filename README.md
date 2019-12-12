# GoDb - a simple data store system for multiple SNP panels (assays)
## Background
For single institution bio-resources genotyping of subjects may have taken place over a period of some years and on differing SNP assay platforms. The resulting data sets would reside in separate files, possibly in different genptype formats (PLINK BED, Oxoford .gen or VCF for example


## Requirements
The over-arching requirement is the provision of a simple system architecture and the software to bring together the genomic data for a population, which originate from multiple SNP assay plaftorms.

### Functional
Some additional functional requirements were set out at the start of the project:
* Curation of genomic data from multiple SNP panels
* Provision of the means to query data using well-known SNP identifiers (dbSNP rsids). 
* Maximising samples size via the automatic combination of genotype records across assay platforms, resolving overlaps in sample sets
* Add genotype data from new assays as they become avaiable.   
 
### Non-Functional
Two non-functional requirements were identified:
* Operate in an environment where compute resources may be limited
* Allow for scaling up of the size of the data, in particular sample-size, with little performance degradation   

## Description
A hybrid data store was built using MongoDb to hold collection of variant, sample and file location data, with genotype data held in VCF files.

Software was developed to take advantage of the rapid access times offered by both MongoDb for storing variant id (rsid) vs genomic co-ordininates, and tabix indexing for random access to compressed VCF files via genomic co-ordinates. This includes a web application to allow querying of the data store by variant_id and lists of variant ids and command line tools for bulk data extract.

 
The following subdirectories can be found in the repository:

- *cfg/* Config files containing environment variables to locate data, software and database host 

- *load/py/* Data store load Python scripts

- *load/sh/* Data store load bash wrapper scripts

- *webapp/* All Python code, templates and image files related to the web application

- *extract/py/* Command-line Python code for genotype data extract 

- *extract/src/* Root directory for Go(lang) source code to build a multi-threaded command-line extract tool

- *extract/sh/* bash scripts, wrappers for command-line extract tools

## Dependencies
## Running

## Data Store Architecture 
![](images/godb_architecture.png)

## Combining genotype records 
![](images/combining_geno_data.png)
## Ackowledgments

