//------------------------------------------------------------------------------
// Access GoDb via the imported mongo libraries and (in 'vcfmerge')
// tabix libraries to produce combined vcf records
//
// Uses goroutines for file i/o, care need to ba taken over the 
// number of open file handles as a long list of SNPs can mean 
// 1000's of file reads
//
// Steps:
// 1) Read in a file of rs numbers
// 2) For each:
//    2.1 get 'variants' and 'filepaths' data from mongodb, access VCF records
//    2.2 save vcf records in maps of arrays (rsid -> array of VCF records
// 3) Organise assaytypes found in the data and build a combined column list and 
//    VCF header record for all present.
// 4) Build combined VCF records, applying genotype resolution rules
// 5) Output combined records.
//------------------------------------------------------------------------------
package main

import (
	"bufio"
	"flag"
	"fmt"
	"genometrics"
	"gopkg.in/mgo.v2"
	"log"
	"godb"
	"os"
	"sample"
	"strings"
	"sync"
	"variant"
	"vcfmerge"
)
//------------------------------------------------
// file-scope vars, accessed by multiple funcs
//------------------------------------------------
var rsFilePath string
var vcfPathPref string
var gdb string
var var_collection string
var fp_collection string
var samp_collection string
var dbhost string
var threshold float64
var assayTypes string
var validAssaytypes = map[string]bool{}
var fsem chan struct{}
var csem chan struct{}
//------------------------------------------------
// main package routines
//------------------------------------------------
//------------------------------------------------
// init() set up and parse cmd line flags
//------------------------------------------------
func init() {
	const (
		defaultRsFilePath  = "./data/rslist1.txt"
		rsusage            = "File containing list of rsnumbers"
		defaultGdb         = "genomicsdb"
		gusage             = "dbname for genomics data"
		defaultvars        = "variants"
		varusage           = "variant collection name"
		defaultfps         = "filepaths"
		fpusage            = "filepath collection name"
		defaultsamps       = "samples"
		sampusage          = "samples collection name"
		defaultDbhost      = "localhost"
		dusage             = "mongodb hostname"
		defaultvcfPathPref = ""
		vusage             = "default path prefix for vcf files"
		defaultThreshold   = 0.9
		thrusage           = "Prob threshold"
		defaultAssayTypes  = "affy,illumina,broad,metabo,exome"
		atusage            = "Assay types"
	)
	flag.StringVar(&rsFilePath, "rsfile", defaultRsFilePath, rsusage)
	flag.StringVar(&rsFilePath, "r", defaultRsFilePath, rsusage+" (shorthand)")
	flag.StringVar(&gdb, "gdb", defaultGdb, gusage)
	flag.StringVar(&gdb, "g", defaultGdb, gusage+" (shorthand)")
	flag.StringVar(&var_collection, "variants", defaultvars, varusage)
	flag.StringVar(&var_collection, "m", defaultvars, varusage+" (shorthand)")
	flag.StringVar(&fp_collection, "filepaths", defaultfps, fpusage)
	flag.StringVar(&fp_collection, "f", defaultfps, fpusage+" (shorthand)")
	flag.StringVar(&samp_collection, "samples", defaultsamps, sampusage)
	flag.StringVar(&samp_collection, "s", defaultsamps, sampusage+" (shorthand)")
	flag.StringVar(&dbhost, "dbhost", defaultDbhost, dusage)
	flag.StringVar(&dbhost, "d", defaultDbhost, dusage+" (shorthand)")
	flag.StringVar(&vcfPathPref, "vcfprfx", defaultvcfPathPref, dusage)
	flag.StringVar(&vcfPathPref, "v", defaultvcfPathPref, dusage+" (shorthand)")
	flag.Float64Var(&threshold, "threshold", defaultThreshold, thrusage)
	flag.Float64Var(&threshold, "t", defaultThreshold, thrusage+" (shorthand)")
	flag.StringVar(&assayTypes, "assaytypes", defaultAssayTypes, atusage)
	flag.StringVar(&assayTypes, "a", defaultAssayTypes, atusage+" (shorthand)")
	flag.Parse()
}
//------------------------------------------------
// check(error) crude general error handler
//------------------------------------------------
func check(e error) {
	//  fmt.Println("check")
	if e != nil {
		fmt.Println("err != nil")
		log.Fatal(e)
	}
}
//------------------------------------------------
// main() program entry point
//------------------------------------------------
func main() {
  var wg sync.WaitGroup
	// Control # of active goroutines (in this case also the # of open files)
	fsem = make(chan struct{}, 64)
	// 10,000 here is arbitrary, needs a re-think
	file_records := make(chan string, 10000)

	f, err := os.Open(rsFilePath)
	check(err)
	defer f.Close()
	session, err := mgo.Dial(dbhost)
	check(err)
	defer session.Close()

	atList := strings.Split(assayTypes, ",")
	for at := range atList {
		validAssaytypes[atList[at]] = true
	}

	rsid_list := make([]string, 0, 1000)
	rsid_count := 0
	scanner := bufio.NewScanner(f)
  // For each rsid, access godb and get the lists of variants vs filepaths
	for scanner.Scan() {
		rsid := scanner.Text()
		rsid_count++
		rsid_list = append(rsid_list, rsid)
    variants, filepaths := godb.Getvardbdata(session, gdb, var_collection, fp_collection, vcfPathPref, rsid)
    // For each variant, filepath combination get a file record
    for idx, variant := range variants {
		  if _, ok := validAssaytypes[variant.Assaytype]; ok {
        //log.Printf("##VAR FILEPATH %v, %s\n", variant, filepaths[idx])
        wg.Add(1)
		    go getvarfiledata(filepaths[idx], variant, file_records, &wg)
      }
    }
	}
	check(err)

	// Wait for the file-reading go routines (defined in the godb package) to complete
	wg.Wait()
	close(fsem)
	close(file_records)

	// Map rsid's to their retrieved vcf file records
	rsids := make(map[string][][]string, rsid_count)
	// And to their vcf data
	rsids_data := make(map[string][]vcfmerge.Vcfdata, rsid_count)

	// What assaytypes do we have?
	assaytypes := make(map[string]int, 10) // 10 is a guess
	assaytype_list := make([]string, 0)

	// Read the channel of file records
	for record := range file_records {
		var recdata vcfmerge.Vcfdata
		fields := strings.Split(record, "\t")
    if _, ok := validAssaytypes[fields[0]]; !ok {
      continue
    }
		if validAssaytypes[fields[0]] != true {
			continue
		}
		prfx, _ := variant.GetVCFPrfx_Sfx(fields[1:])
		recdata.Probidx = variant.GetProbidx(prfx)
		if _, ok := assaytypes[fields[0]]; !ok {
			assaytypes[fields[0]] = 1
			assaytype_list = append(assaytype_list, fields[0])
		}
		varid := variant.GetVarid(prfx)
		if _, ok := rsids[varid]; !ok {
			rsids[varid] = make([][]string, 0)
			rsids_data[varid] = make([]vcfmerge.Vcfdata, 0)
		}
		rsids[varid] = append(rsids[varid], fields)
		rsids_data[varid] = append(rsids_data[varid], recdata)
	}
  // get all sample data from godb and organise into maps of maps: 
  // assaytype -> sample name -> sample posn (sample_name_map)
  // assaytype -> sample posn -> sample name (sample_posn_map)
	sample_name_map, sample_posn_map := godb.GetSamplesByAssaytype(session, gdb, samp_collection)
  // Condense all sample_names into a combined map samplename -> record position
	combocols := sample.GetCombinedSampleMapByAssaytypes(sample_name_map, assaytype_list)
	// Get column headers as a single tab delimited string, with prefix in place, and as a list, both in postion order
	colhdr_str, combo_names := vcfmerge.GetCombinedColumnHeaders(combocols)
	//fmt.Printf("%s\n", "combined"+"\t"+colhdr_str)
	fmt.Printf("%s\n", colhdr_str)

	var genomet genometrics.AllMetrics

	// output the vcf records in input order, can also log the 'NOT FOUND's at this point
	for _, rsid := range rsid_list {
		if records, ok := rsids[rsid]; ok {
      rec_str := vcfmerge.Combine_one(records, rsids_data[rsid], rsid, sample_posn_map, combocols, combo_names, threshold, &genomet)
      //fmt.Printf("%s\n", "combined"+"\t"+rec_str)
      fmt.Printf("%s\n", rec_str)
		}
	}
	log.Printf("##METRICS (ALL),AllGenos=%d,UniqueGenos=%d,Alloverlap=%d,Two=%d,GTTwo=%d,Odiff=%d,OMiss=%d,OMissRes=%d,NoAssay=%d\n", 
      genomet.AllGenoCount, genomet.UniqueGenoCount, genomet.OverlapTestCount, 
      genomet.TwoOverlapCount, genomet.GtTwoOverlapCount, genomet.MismatchCount, 
      genomet.MissTestCount, genomet.MissingCount, genomet.NoAssayCount)
}
//------------------------------------------------------------------------------
// wrap the godb.Getvarfiledata func, for use as a goroutine
//------------------------------------------------------------------------------
func getvarfiledata(f string, dbv godb.DBVariant, recs chan string, wg *sync.WaitGroup) {
  fsem <- struct{}{}
  defer func() { <-fsem }()
  godb.Getvarfiledata(f, dbv, recs)
  wg.Done()
}
