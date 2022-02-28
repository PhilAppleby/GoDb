//--------------------------------------------------------------------------------------
// Merge operation on 2 VCF files covering the same genomic range
// genofile contains genotypes, annotfile contains annotations
// Key match loop conditions:
// Keys match -> update the VCF data with the rsid from the annot file, write to vcf output,
//    read from both the genofile and the annotfile
// Annot file key low -> read from the annotfile
// TODO, to be completed
//--------------------------------------------------------------------------------------
// Author: P Appleby
//--------------------------------------------------------------------------------------
package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"genometrics"
	"log"
	"os"
	"sort"
	"strings"
	"variant"
	"vcfmerge"
)

// min_posn = 0
// maxPosn is intended to be greater than any value for genomic position
const maxPosn int64 = 999999999999

var emptyRecord = []string{}

//-----------------------------------------------
// global vars, accessed by multiple funcs
//-----------------------------------------------
var logFilePath string
var vcfFilePath string
var annotFilePath string
var chr string

//-----------------------------------------------
// main package routines
//-----------------------------------------------
func init() {
	const (
		defaultLogFilePath   = "./data/filemergevcf_output.log"
		lusage               = "Log file"
		defaultVcfFilePath   = "./data/test.vcf.gz"
		vusage               = "path for vcf file"
		defaultAnnotFilePath = "./data/annot.vcf.gz"
		ausage               = "path for annot file"
		defaultChr           = "22"
		chrusage             = "default chromosome (number as string)"
	)
	flag.StringVar(&logFilePath, "logfile", defaultLogFilePath, lusage)
	flag.StringVar(&logFilePath, "l", defaultLogFilePath, lusage+" (shorthand)")
	flag.StringVar(&vcfFilePath, "vcffile", defaultVcfFilePath, vusage)
	flag.StringVar(&vcfFilePath, "v", defaultVcfFilePath, vusage+" (shorthand)")
	flag.StringVar(&annotFilePath, "annotfile", defaultAnnotFilePath, ausage)
	flag.StringVar(&annotFilePath, "v", defaultAnnotFilePath, ausage+" (shorthand)")
	flag.StringVar(&chr, "chr", defaultChr, chrusage)
	flag.StringVar(&chr, "c", defaultChr, chrusage+" (shorthand)")
	flag.Parse()
}

func check(e error) {
	if e != nil {
		log.Println("err != nil")
		log.Fatal(e)
	}
}

// Main entry point to program – most of the following
// lines are for setting up and initializing the run
// before repeatedly calling the write and read functions
func main() {
	// set up logging to a file
	lf, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	check(err)
	defer lf.Close()

	log.SetOutput(lf)
	log.Printf("START merge %s, %s\n", vcfFilePath, annotFilePath)

	fhv, err := os.Open(vcfFilePath)
	check(err)
	defer fhv.Close()
	grv, err := gzip.NewReader(fhv)
	check(err)
	defer grv.Close()
	readerv := bufio.NewReader(grv)
	// print VCF file headers

	fha, err := os.Open(annotFilePath)
	check(err)
	defer fha.Close()
	gra, err := gzip.NewReader(fha)
	check(err)
	defer gra.Close()
	readerv := bufio.NewReader(gra)
	// go past annot file headers

}

//-------------------------------------------------------------
// Read a record from a single reader and split it to format a string slice
//-------------------------------------------------------------
func getNextRecordSlice(rdr *bufio.Reader) ([]string, int64, string) {
	posn := maxPosn
	data := emptyRecord
	varid := ""
	text, err := rdr.ReadString('\n')
	if err == nil {
		text = strings.TrimRight(text, "\n")
		data = strings.Split(text, "\t")
		posn = int64(variant.GetPosn(data))
		varid = variant.GetVarid(data)
	}
	return data, posn, varid
}

//-------------------------------------------------------------
// If all keys are high there are no records to be read
//-------------------------------------------------------------
func recordsRemain(keys map[string]int64) bool {
	for _, pos := range keys {
		if pos < maxPosn {
			return true
		}
	}
	return false
}

//-------------------------------------------------------------
// Write output from low-key records
//-------------------------------------------------------------
func outputFromLowKeyRecords(records map[string][]string, keys map[string]int64,
	samplePosnMap map[string]map[int]string,
	combocols map[string]int, comboNames []string, threshold float64, genomet *genometrics.AllMetrics) {
	//
	lowKeys := getLowKeys(keys)
	lowKeyAt := make([]string, 0)

	for at := range lowKeys {
		lowKeyAt = append(lowKeyAt, at)
	}
	sort.Sort(sort.Reverse(sort.StringSlice(lowKeyAt)))

	vcfrecords := make([][]string, 0, len(records))
	rsid := ""
	for _, at := range lowKeyAt {
		rsid = variant.GetVarid(records[at])
		rec := make([]string, 1, len(records[at])+1)
		rec[0] = at
		rec = append(rec, records[at]...)
		vcfrecords = append(vcfrecords, rec)
		//fmt.Printf("LOWKEY OUTPUT %s, %d\n", at, key)
	}
	var vcfd []vcfmerge.Vcfdata
	var rsidGenomet genometrics.AllMetrics
	recStr := vcfmerge.CombineOne(vcfrecords, vcfd, rsid, samplePosnMap, combocols, comboNames, threshold, &rsidGenomet)
	genometrics.Increment(genomet, &rsidGenomet)
	errorPct := 0.0
	if rsidGenomet.MismatchCount > 0 {
		errorPct = (float64(rsidGenomet.MismatchCount) / float64(rsidGenomet.OverlapTestCount)) * 100
	}
}

//-------------------------------------------------------------
// Inititiate the next cycle
//-------------------------------------------------------------
func readFromLowKeyRecords(records map[string][]string, keys map[string]int64, rdrs map[string]*bufio.Reader, varids map[string]string) (map[string][]string, map[string]int64, map[string]string) {
	lowKeys := getLowKeys(keys)
	for assaytype := range lowKeys {
		records[assaytype], keys[assaytype], varids[assaytype] = getNextRecordSlice(rdrs[assaytype])
	}
	return records, keys, varids
}

func getLowKeys(keys map[string]int64) map[string]int64 {
	lowKeys := make(map[string]int64)
	lowKey := maxPosn
	for _, pos := range keys {
		if pos < lowKey {
			lowKey = pos
		}
	}
	if lowKey < maxPosn {
		for at, pos := range keys {
			if pos == lowKey {
				lowKeys[at] = pos
			}
		}
	}
	return lowKeys
}
