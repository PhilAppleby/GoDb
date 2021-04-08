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
//    2.2 save vcf records in maps of arrays
// 3) Organise assaytypes found in the data and build a combined column list
//    for all present.
// 4) Build combined VCF records, applying genotype resolution rules
// 5) Output combined records.
//------------------------------------------------------------------------------
package main

import (
	"bufio"
	"flag"
	"fmt"
	"genometrics"
	"godb"
	"log"
	"os"
	"sample"
	"strings"
	"variant"
	"vcfmerge"
)

//------------------------------------------------
// file-scope vars, accessed by multiple funcs
//------------------------------------------------
var rsFilePath string
var vcfPathPref string
var threshold float64
var assayTypes string
var validAssaytypes = map[string]bool{}

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
		defaultvcfPathPref = ""
		vusage             = "default path prefix for vcf files"
		defaultThreshold   = 0.9
		thrusage           = "Prob threshold"
		defaultAssayTypes  = "affy,illumina,broad,metabo,exome"
		atusage            = "Assay types"
	)
	flag.StringVar(&rsFilePath, "rsfile", defaultRsFilePath, rsusage)
	flag.StringVar(&rsFilePath, "r", defaultRsFilePath, rsusage+" (shorthand)")
	flag.StringVar(&vcfPathPref, "vcfprfx", defaultvcfPathPref, vusage)
	flag.StringVar(&vcfPathPref, "v", defaultvcfPathPref, vusage+" (shorthand)")
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

func main() {
	fileRecords := make(chan string, 10000)

	f, err := os.Open(rsFilePath)
	check(err)
	defer f.Close()

	atList := strings.Split(assayTypes, ",")
	//fmt.Printf("%v\n", atList)
	for at := range atList {
		validAssaytypes[atList[at]] = true
	}

	rsidList := make([]string, 0, 1000)
	rsidCount := 0
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		rsid := scanner.Text()
		rsidCount++
		rsidList = append(rsidList, rsid)
		variants, filepaths := godb.Getvardbdata(vcfPathPref, rsid)
		for idx, variant := range variants {
			if _, ok := validAssaytypes[variant.Assaytype]; ok {
				godb.Getvarfiledata(filepaths[idx], variant, fileRecords)
			}
		}
	}
	check(err)

	close(fileRecords)

	// Map rsid's to their retrieved vcf file records
	rsids := make(map[string][][]string, rsidCount)
	// And to their vcf data
	rsidsData := make(map[string][]vcfmerge.Vcfdata, rsidCount)

	// What assaytypes do we have?
	assaytypes := make(map[string]int, 10) // 10 is a guess
	assaytypeList := make([]string, 0)

	// Read the channel of file records
	for record := range fileRecords {
		var recdata vcfmerge.Vcfdata
		fields := strings.Split(record, "\t")
		if _, ok := validAssaytypes[fields[0]]; !ok {
			continue
		}
		if validAssaytypes[fields[0]] != true {
			continue
		}
		prfx, _ := variant.GetVCFPrfxSfx(fields[1:])
		recdata.Probidx = variant.GetProbIdx(prfx)
		if _, ok := assaytypes[fields[0]]; !ok {
			assaytypes[fields[0]] = 1
			assaytypeList = append(assaytypeList, fields[0])
		}
		varid := variant.GetVarid(prfx)
		if _, ok := rsids[varid]; !ok {
			rsids[varid] = make([][]string, 0)
			rsidsData[varid] = make([]vcfmerge.Vcfdata, 0)
		}
		rsids[varid] = append(rsids[varid], fields)
		rsidsData[varid] = append(rsidsData[varid], recdata)
	}

	// Process sample data to:
	// - build header columns
	// - get maps to go from source column numbers to numbers in the combined version
	sampleNameMap, samplePosnMap := godb.GetSamplesByAssaytype()
	combocols := sample.GetCombinedSampleMapByAssaytypes(sampleNameMap, assaytypeList)
	// combocols := sample.GetCombinedSampleMap(sampleNameMap)

	// Get column headers
	colhdrStr, comboNames := vcfmerge.GetCombinedColumnHeaders(combocols)
	//fmt.Printf("%s\n", "combined"+"\t"+colhdrStr)
	fmt.Printf("%s\n", colhdrStr)

	//for _, atype := range assaytypeList {
	//	name_str := vcfmerge.GetColumnHeaders(samplePosnMap[atype])
	//fmt.Printf("%s\n", atype+"\t"+name_str)
	//}
	var genomet genometrics.AllMetrics

	// output the vcf records in input order, can also log the 'NOT FOUND's at this point
	//fmt.Printf("METRICS,platform,rsid,CR,RAF,AAF,MAF,HWEP,HET,COMMON,RARE,N,MISS,DOT,REFPAF,OK\n")
	for _, rsid := range rsidList {
		if records, ok := rsids[rsid]; ok {
			recStr := vcfmerge.CombineOne(records, rsidsData[rsid], rsid, samplePosnMap, combocols, comboNames, threshold, &genomet)
			//fmt.Printf("%s\n", "combined"+"\t"+recStr)
			fmt.Printf("%s\n", recStr)
			// output individual assay records
			//for _, rec := range records {
			//	assaytype := rec[0]
			//	cr, raf, aaf, maf, hwep, het, common, rare, n, miss, dot, refpaf := genometrics.Metrics_for_record(rec[1:], threshold)
			//	flag_str := ""
			//	if hwep < 0.05 {
			//		flag_str = "***"
			//	}
			//fmt.Printf("METRICS,%s,%s,%f,%f,%f,%f,%.6f,%d,%d,%d,%d,%d,%d,%f,%s\n",
			//	assaytype, rsid, cr, raf, aaf, maf, hwep, het, common, rare, n, miss, dot, refpaf, flag_str)
			//	recStr := strings.Join(rec, "\t")
			//fmt.Printf("%s\n", recStr)
			//}
		}
	}
	log.Printf("##METRICS (ALL),AllGenos=%d,UniqueGenos=%d,Alloverlap=%d,Two=%d,GTTwo=%d,Odiff=%d,OMiss=%d,OMissRes=%d,NoAssay=%d\n",
		genomet.AllGenoCount, genomet.UniqueGenoCount, genomet.OverlapTestCount,
		genomet.TwoOverlapCount, genomet.GtTwoOverlapCount, genomet.MismatchCount,
		genomet.MissTestCount, genomet.MissingCount, genomet.NoAssayCount)
}
