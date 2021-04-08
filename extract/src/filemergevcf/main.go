//--------------------------------------------------------------------------------------
// Merge 2 to n VCF files covering the same genomic range
//-------------------------------------------------------
// The are three major aspects:
// 1) Handling file merge as an extension of a two-way merge - this can be regarded
// as an implemetation of 'direct k-way' merge operating on files
//
// 2) For low-key records sets, record combination takes place if set-size > 0
//
// args:
//  --tpltfile: a text file of template file paths for files to be merged
//  --paramfile: file of parameters for genotype resolution
//  --chr: chromosome
//  --logfile: full filepath for logging
//  --vcfprfx: directory root for vcf files
//
//  Author: P Appleby, University of Dundee
//--------------------------------------------------------------------------------------
package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"genometrics"
	"io"
	"log"
	"os"
	"sample"
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
var tpltFilePath string
var paramFilePath string
var logFilePath string
var vcfPathPref string
var chr string
var threshold float64
var errpctthr float64

//-----------------------------------------------
// main package routines
//-----------------------------------------------
func init() {
	const (
		defaultTpltFilePath  = "./data/vcf_file_template.txt"
		tusage               = "File template strings"
		defaultParamFilePath = "./data/params.cfg"
		pusage               = "QC Parameter file"
		defaultLogFilePath   = "./data/filemergevcf_output.log"
		lusage               = "Log file"
		defaultvcfPathPref   = "/var/data"
		vusage               = "default path prefix for vcf files"
		defaultErrPct        = 10.0
		epctusage            = "ErrPct threshold"
		defaultThreshold     = 0.9
		thrusage             = "Prob threshold"
		defaultChr           = "22"
		chrusage             = "default chromosome (number as string)"
	)
	flag.StringVar(&tpltFilePath, "tpltfile", defaultTpltFilePath, tusage)
	flag.StringVar(&tpltFilePath, "t", defaultTpltFilePath, tusage+" (shorthand)")
	flag.StringVar(&paramFilePath, "paramfile", defaultParamFilePath, pusage)
	flag.StringVar(&paramFilePath, "p", defaultParamFilePath, pusage+" (shorthand)")
	flag.StringVar(&logFilePath, "logfile", defaultLogFilePath, lusage)
	flag.StringVar(&logFilePath, "l", defaultLogFilePath, lusage+" (shorthand)")
	flag.StringVar(&vcfPathPref, "vcfprfx", defaultvcfPathPref, vusage)
	flag.StringVar(&vcfPathPref, "v", defaultvcfPathPref, vusage+" (shorthand)")
	flag.Float64Var(&threshold, "threshold", defaultThreshold, thrusage)
	flag.Float64Var(&threshold, "h", defaultThreshold, thrusage+" (shorthand)")
	flag.Float64Var(&errpctthr, "errpct", defaultErrPct, epctusage)
	flag.Float64Var(&errpctthr, "e", defaultErrPct, epctusage+" (shorthand)")
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
	log.Printf("START merge %s, %s\n", paramFilePath, tpltFilePath)

	// Load file templates
	f, err := os.Open(tpltFilePath)
	check(err)
	defer f.Close()
	scanner := bufio.NewScanner(f)
	assaytypeFilename := make(map[string]string)
	assaytypeList := make([]string, 0)
	fscanners := make(map[string]*bufio.Scanner)
	freaders := make(map[string]*bufio.Reader)

	for scanner.Scan() {
		text := scanner.Text()
		if !strings.HasPrefix(text, "#") {
			fields := strings.Split(text, "=")
			assaytypeFilename[fields[0]] = fmt.Sprintf(fields[1], vcfPathPref, chr)
			assaytypeList = append(assaytypeList, fields[0])
		}
	}
	// Load QC parameters (if present)
	fp, err := os.Open(paramFilePath)
	check(err)
	defer fp.Close()
	testnum := ""
	callrate := ""
	mafdelta := ""
	infoscore := ""
	scanner = bufio.NewScanner(fp)
	for scanner.Scan() {
		text := scanner.Text()
		if !strings.HasPrefix(text, "#") {
			fields := strings.Split(text, "=")
			if fields[0] == "TESTNUM" {
				testnum = fields[1]
			}
			if fields[0] == "CALLRATE" {
				callrate = fields[1]
			}
			if fields[0] == "MAFDELTA" {
				mafdelta = fields[1]
			}
			if fields[0] == "INFOSCORE" {
				infoscore = fields[1]
			}
		}
	}
	runParams := genometrics.GetRunParams(testnum, mafdelta, callrate, infoscore)
	log.Printf("Params: %v\n", runParams)

	for key, value := range assaytypeFilename {
		fh, err := os.Open(value)
		check(err)
		defer fh.Close()
		gr, err := gzip.NewReader(fh)
		check(err)
		defer gr.Close()
		scanner := bufio.NewScanner(gr)
		reader := bufio.NewReader(gr)
		fscanners[key] = scanner
		freaders[key] = reader
	}
	// handle file headers
	headers := make(map[string][]string)
	for assaytype, rdr := range freaders {
		headers[assaytype], _ = getSampleHeaders(rdr)
		//fmt.Printf("%s hdr len = %d\n", assaytype, len(headers[assaytype]))
	}
	// Headers and combined header map
	sampleNameMap, samplePosnMap := sample.MakeSamplesByAssaytype(headers)
	combocols := sample.GetCombinedSampleMap(sampleNameMap)
	// combocols := sample.GetCombinedSampleMapByAssaytypes(sampleNameMap, assaytypeList)
	colhdrStr, comboNames := vcfmerge.GetCombinedColumnHeaders(combocols)
	//fmt.Printf("%s\n", "combined"+"\t"+colhdrStr)
	printHeaders()
	fmt.Printf("%s\n", colhdrStr)

	// read first records and capture keys (genomic positions)
	records := make(map[string][]string)
	keys := make(map[string]int64)
	varids := make(map[string]string)
	for assaytype, rdr := range freaders {
		records[assaytype], keys[assaytype], varids[assaytype] = getNextRecordSlice(rdr)
	}
	var genomet genometrics.AllMetrics
	outctr := 0
	// process until all files exhausted
	for recordsRemain(keys) {
		outputFromLowKeyRecords(records, keys, samplePosnMap, combocols, comboNames, threshold, &genomet)
		outctr++
		records, keys, varids = readFromLowKeyRecords(records, keys, freaders, varids)
	}
	errorPct := (float64(genomet.MismatchCount) / float64(genomet.OverlapTestCount)) * 100
	log.Printf("EXIT,wrt=%d,AllGenos=%d,UniqueGenos=%d,Alloverlap=%d,Two=%d,GTTwo=%d\n",
		outctr, genomet.AllGenoCount, genomet.UniqueGenoCount, genomet.OverlapTestCount,
		genomet.TwoOverlapCount, genomet.GtTwoOverlapCount)
	log.Printf("EXIT,OverlapGenoDiffs=%d,DiffProbDiffs=%d, SameProbDiffs=%d,MissingGenoTested=%d,MissingUnresolved=%d,NoAssay=%d,ErrorPct=%.3f\n",
		genomet.MismatchCount, genomet.DiffProbDiffs, genomet.SameProbDiffs,
		genomet.MissTestCount, genomet.MissingCount, genomet.NoAssayCount, errorPct)
}

//-------------------------------------------------------------
// Get headers with column(sample) names
//-------------------------------------------------------------
func getSampleHeaders(rdr *bufio.Reader) ([]string, error) {
	var sfx []string

	eof := false
	for !eof {
		text, err := rdr.ReadString('\n')
		if err == io.EOF {
			return emptyRecord, err
		}
		text = strings.TrimRight(text, "\n")
		if strings.HasPrefix(text, "##") {
			//fmt.Printf("%s\n", text)
		} else {
			if strings.HasPrefix(text, "#") {
				_, sfx = variant.GetVCFPrfxSfx(strings.Split(text, "\t"))
				break
			}
		}
	}
	return sfx, nil
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
	if errorPct < errpctthr {
		fmt.Printf("%s\n", recStr)
	} else {
		genometrics.LogMetrics(1, rsid, 1, "##ERRPCT", &rsidGenomet)
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

func printHeaders() {
	fmt.Printf("%s\n", "##fileformat=VCFv4.2")
	fmt.Printf("%s\n", "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">")
	fmt.Printf("%s\n", "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">")
	fmt.Printf("%s\n", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
	fmt.Printf("%s\n", "##INFO=<ID=RefPanelAF,Number=A,Type=Float,Description=\"Allele frequency in imputation reference panel\">")
	fmt.Printf("%s\n", "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype dosage\">")
	fmt.Printf("%s\n", "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype posterior probabilities\">")
	fmt.Printf("%s\n", "##FORMAT=<ID=AT,Number=1,Type=String,Description=\"Assay Type\">")
	fmt.Printf("%s\n", "##INFO=<ID=TYPED,Number=0,Type=Flag,Description=\"Typed in input data\">")
}
