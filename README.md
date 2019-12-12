# GoDb - a simple, hybrid data store system for multiple SNP panels (assays)
## Background
For single institution bio-resources genotyping of subjects may have taken place over a period of some years and on differing SNP assay platforms. The resulting data sets would reside in separate files, possibly in different genptype formats (PLINK BED, Oxford .gen or VCF, for example).


## Requirements
The over-arching requirement is the provision of a simple system architecture and the software to bring together the genomic data for a population, which originate from multiple SNP assay plaftorms.

### Functional
Some additional functional requirements were set out at the start of the project:
* Curation of genomic data from multiple SNP panels.
* Provision of the means to query data using well-known SNP identifiers (dbSNP rsids). 
* Maximising samples size via the automatic combination of genotype records across assay platforms, resolving overlaps in sample sets
* Add genotype data from new assays as they become available.   
 
### Non-Functional
Two non-functional requirements were identified:
* Operate in an environment where compute resources may be limited.
* Allow for scaling up of the size of the data, in particular sample-size, with little performance degradatioa.n   

## Description
A hybrid data store was designed and built using MongoDb to hold collections of variant, sample and file location data, with genotype data held in VCF files.

Software was developed to take advantage of the rapid access times offered by both MongoDb for storing variant id (rsid) vs genomic co-ordininates, and tabix indexing for random access to compressed VCF files via genomic co-ordinates. This includes a web application to allow querying of the data store by variant_id and lists of variant ids and command line tools for bulk data extract.

### This repository 
The following subdirectories can be found in the repository:

- *cfg/* Config files containing environment variables to locate data, software and database host 

- *load/py/* Data store load Python scripts

- *load/sh/* Data store load bash wrapper scripts

- *webapp/* All Python code, templates and image files related to the web application

- *extract/py/* Command-line Python code for genotype data extract 

- *extract/src/* Root directory for Go(lang) source code to build a multi-threaded command-line extract tool

- *extract/sh/* bash scripts, wrappers for command-line extract tools

### MongoDb Collection Examples
variants (one document per variant per SNP panel (assaytype)):
```
{
	"_id" : ObjectId("5dee1245d5d298277178bcfb"),
	"assaytype" : "metabo",
	"ref_maf" : 0.6,
	"rsid" : "rs7294904",
	"info" : 1,
	"alleleB" : "C",
	"position" : 199532,
	"alleleA" : "T",
	"chromosome" : "12"
}
```

samples (one document per sample per SNP panel (assaytype)):
```
{
	"_id" : ObjectId("5decf26e64b5031da4b9c5cd"),
	"assaytype" : "affy",
	"list_posn" : 0,
	"sample_id" : "006561"
}
```

filepaths (one document per SNP panel (assaytype)):
```
{
	"_id" : ObjectId("5decf307339f134beefc2419"),
	"assaytype" : "affy",
	"files" : [
		{
			"CHROM" : "12",
			"filename" : "chr12.vcf.gz"
		},
		{
			"CHROM" : "02",
			"filename" : "chr02.vcf.gz"
		},
		{
			"CHROM" : "06",
			"filename" : "chr06.vcf.gz"
		},
		{
			"CHROM" : "22",
			"filename" : "chr22.vcf.gz"
		},
		{
			"CHROM" : "01",
			"filename" : "chr01.vcf.gz"
		},
		{
			"CHROM" : "19",
			"filename" : "chr19.vcf.gz"
		}
	],
	"fpath_prefix" : "<data root >/godb",
	"filepath" : "<data root >/affy/",
	"fpath_suffix" : "affy/"
}
```



## Dependencies
- MongoDb community edition (version 3 upwards, tested to 4.2.1)
- Tabix (0.2.5)
- pysam python library (0.15.3)
- pymongo python library (3.9.0)
- gopkg.in/mgo.v2 golang library
- gopkg.in/mgo.v2/bson golang library
- github.com/brentp/bix golang library for tabix index access
- github.com/brentp/irelate/intercases golang library for tabix index access

## Data Store Architecture 
High-level block diagram of the data store and software layers

![](images/godb_architecture.png)

## Combining genotype records 
![](images/combining_geno_data.png)
## Ackowledgments

