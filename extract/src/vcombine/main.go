//------------------------------------------------------------------------------
// Access GoDb via the imported mongo libraries and (in 'vcfmerge')
// tabix libraries to produce combined vcf records
//
// Uses goroutines for file i/o, care need to ba taken over the 
// number of open file handles as a long list of SNPs can mean 
// 1000's of file reads
//
// Steps:
// 1) Read in a file of rs numbers or a single rsid
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
var logFilePath string
var rsFilePath string
var rsId string
var vcfPathPref string
var threshold float64
var errpctthr float64
var assayTypes string
var logLevel int
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
    defaultLogFilePath   = "./data/vcombine_output.log"
    lusage               = "Log file"
		defaultRsFilePath  = "./data/rslist1.txt"
		rsusage            = "File containing list of rsnumbers"
		defaultRsId        = ""
		rsidusage          = "Single rsid"
		defaultvcfPathPref = ""
		vusage             = "default path prefix for vcf files"
		defaultThreshold   = 0.9
		thrusage           = "Prob threshold"
		defaultAssayTypes  = "affy,illumina,broad,metabo,exome"
		atusage            = "Assay types"
		defaultLogLevel    = 0
		loglusage          = "0=Minimal 1=Sum 2=max"
	)
  flag.StringVar(&logFilePath, "logfile", defaultLogFilePath, lusage)
  flag.StringVar(&logFilePath, "l", defaultLogFilePath, lusage+" (shorthand)")
	flag.StringVar(&rsFilePath, "rsfile", defaultRsFilePath, rsusage)
	flag.StringVar(&rsFilePath, "r", defaultRsFilePath, rsusage+" (shorthand)")
	flag.StringVar(&rsId, "rsid", defaultRsId, rsidusage)
	flag.StringVar(&rsId, "i", defaultRsId, rsidusage+" (shorthand)")
	flag.StringVar(&vcfPathPref, "vcfprfx", defaultvcfPathPref, vusage)
	flag.StringVar(&vcfPathPref, "v", defaultvcfPathPref, vusage+" (shorthand)")
	flag.Float64Var(&threshold, "threshold", defaultThreshold, thrusage)
	flag.Float64Var(&threshold, "t", defaultThreshold, thrusage+" (shorthand)")
	flag.StringVar(&assayTypes, "assaytypes", defaultAssayTypes, atusage)
	flag.StringVar(&assayTypes, "a", defaultAssayTypes, atusage+" (shorthand)")
	flag.IntVar(&logLevel, "logopt", defaultLogLevel, loglusage)
	flag.IntVar(&logLevel, "o", defaultLogLevel, loglusage+" (shorthand)")
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
  // set up logging to a file, TODO move to init()
  lf, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0640)
  check(err)
  defer lf.Close()

  log.SetOutput(lf)
  var wg sync.WaitGroup
	// Control # of active goroutines (in this case also the # of open files)
	fsem = make(chan struct{}, 64)
	// 10,000 here is arbitrary, needs a re-think
	file_records := make(chan string, 10000)

	rsid_list := make([]string, 0, 1000)
	rsid_count := 0

  if rsId == "" {
	  f, err := os.Open(rsFilePath)
	  check(err)
	  defer f.Close()
	  scanner := bufio.NewScanner(f)
	  for scanner.Scan() {
		  rsid := scanner.Text()
		  rsid_count++
		  rsid_list = append(rsid_list, rsid)
    }
  } else {
		rsid_list = append(rsid_list, rsId)
  }

	atList := strings.Split(assayTypes, ",")
	for at := range atList {
		validAssaytypes[atList[at]] = true
	}

  for _, rsid := range rsid_list {
    // For each rsid, access godb and get the lists of variants vs filepaths
    variants, filepaths := godb.Getvardbdata(vcfPathPref, rsid)
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
	sample_name_map, sample_posn_map := godb.GetSamplesByAssaytype()
  // Condense all sample_names into a combined map samplename -> record position
	combocols := sample.GetCombinedSampleMapByAssaytypes(sample_name_map, assaytype_list)
	// Get column headers as a single tab delimited string, with prefix in place, and as a list, both in postion order
	colhdr_str, combo_names := vcfmerge.GetCombinedColumnHeaders(combocols)
	//fmt.Printf("%s\n", "combined"+"\t"+colhdr_str)
	fmt.Printf("%s\n", colhdr_str)

	var genomet genometrics.AllMetrics

	// output the vcf records in input order, can also log the 'NOT FOUND's at this point
  var snpcount int
	for _, rsid := range rsid_list {
		if records, ok := rsids[rsid]; ok {
	    var rsid_genomet genometrics.AllMetrics
      rec_str := vcfmerge.Combine_one(records, rsids_data[rsid], rsid, sample_posn_map, combocols, combo_names, threshold, &rsid_genomet)
      fmt.Printf("%s\n", rec_str)
      snpcount  += 1
      genometrics.Increment(&genomet, &rsid_genomet)
      genometrics.Log_metrics(logLevel, rsid, 1, "##VARIANT", &rsid_genomet)
		}
	}
  genometrics.Log_metrics(3, "all", snpcount, "##TOTAL", &genomet)
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
