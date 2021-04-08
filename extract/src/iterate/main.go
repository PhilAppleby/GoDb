//------------------------------------------------------------------------------
// Test test application to iterate over a list of rsid's, to
// Get genotype data from the GoDb data store, main pupose is to time
// iterations to produce performance statistics.
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
	"time"
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
	start := time.Now()
	fileRecords := make(chan string, 1)

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
	var variants []godb.DBVariant
	var filepaths []string
	for scanner.Scan() {
		rsid := scanner.Text()
		rsidCount++
		rsidList = append(rsidList, rsid)
		loopStart := time.Now()
		for i := 0; i < 1000; i++ {
			variants, filepaths = godb.Getvardbdata(vcfPathPref, rsid)
		}
		elapsed := time.Since(loopStart)
		log.Printf("Mongo Iteration took %s", elapsed)
		loopStart = time.Now()
		for i := 0; i < 100; i++ {
			fileRecords = make(chan string, 100000)
			for idx, variant := range variants {
				if _, ok := validAssaytypes[variant.Assaytype]; ok {
					godb.Getvarfiledata(filepaths[idx], variant, fileRecords)
				}
			}
		}
		elapsed = time.Since(loopStart)
		log.Printf("VCF Iteration took %s", elapsed)
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
	// combocols := sample.GetCombinedSampleMap(sample_name_map)

	// Get column headers
	colhdrStr, comboNames := vcfmerge.GetCombinedColumnHeaders(combocols)
	//fmt.Printf("%s\n", "combined"+"\t"+colhdr_str)
	fmt.Printf("%s\n", colhdrStr)

	//for _, atype := range assaytype_list {
	//	name_str := vcfmerge.GetColumnHeaders(sample_posn_map[atype])
	//fmt.Printf("%s\n", atype+"\t"+name_str)
	//}
	var genomet genometrics.AllMetrics

	// output the vcf records in input order, can also log the 'NOT FOUND's at this point
	//fmt.Printf("METRICS,platform,rsid,CR,RAF,AAF,MAF,HWEP,HET,COMMON,RARE,N,MISS,DOT,REFPAF,OK\n")
	for _, rsid := range rsidList {
		if records, ok := rsids[rsid]; ok {
			loopStart := time.Now()
			for i := 0; i < 1000; i++ {
				_ = vcfmerge.CombineOne(records, rsidsData[rsid], rsid, samplePosnMap, combocols, comboNames, threshold, &genomet)
			}
			elapsed := time.Since(loopStart)
			log.Printf("Combo Iteration took %s", elapsed)
			//fmt.Printf("%s\n", "combined"+"\t"+rec_str)
			//fmt.Printf("%s\n", rec_str)
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
			//	rec_str := strings.Join(rec, "\t")
			//fmt.Printf("%s\n", rec_str)
			//}
		}
	}
	log.Printf("##METRICS (ALL),AllGenos=%d,UniqueGenos=%d,Alloverlap=%d,Two=%d,GTTwo=%d,Odiff=%d,OMiss=%d,OMissRes=%d,NoAssay=%d\n",
		genomet.AllGenoCount, genomet.UniqueGenoCount, genomet.OverlapTestCount,
		genomet.TwoOverlapCount, genomet.GtTwoOverlapCount, genomet.MismatchCount,
		genomet.MissTestCount, genomet.MissingCount, genomet.NoAssayCount)
	allElapsed := time.Since(start)
	log.Printf("END: %s", allElapsed)
}
