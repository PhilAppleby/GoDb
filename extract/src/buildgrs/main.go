//------------------------------------------------------------------------------
// Access GoDb via the API and use the results to build Genetic Risk scores
// Steps:
// 1) Read in a file of rsids and GRS information, including effect allele and
// wgt to apply
// 2) For each:
//------------------------------------------------------------------------------
package main

import (
	"bufio"
	"flag"
	"fmt"
	"godb"
	"grs"
	"log"
	"os"
	"strconv"
	"strings"
)

//------------------------------------------------
// file-scope vars, accessed by multiple funcs
//------------------------------------------------
var logFilePath string
var rsFilePath string
var vcfPathPref string
var threshold float64
var errpctthr float64
var assayTypes string
var logLevel int
var validAssaytypes = map[string]bool{}

//------------------------------------------------
// main package routines
//------------------------------------------------
//------------------------------------------------
// init() set up and parse cmd line flags
//------------------------------------------------
func init() {
	const (
		defaultLogFilePath = "./logs/buildgrs_output.log"
		lusage             = "Log file"
		defaultRsFilePath  = "./data/grs_list.txt"
		rsusage            = "File containing list of rsnumbers, effect allele, weight"
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

	grsList := make([]string, 0, 1000)
	grsInCount := 0

	f, err := os.Open(rsFilePath)
	check(err)
	defer f.Close()
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		grsLine := scanner.Text()
		grsInCount++
		grsList = append(grsList, grsLine)
	}

	atList := strings.Split(assayTypes, ",")
	for at := range atList {
		validAssaytypes[atList[at]] = true
	}

	rsidList, eaMap, eafMap, wgtMap := getGrsLists(grsList)

	_, _, genorecs := godb.Getallvardata(vcfPathPref, rsidList, validAssaytypes, threshold)

	grScores := grs.GetScores(genorecs, eaMap, eafMap, wgtMap)

	for _, gScore := range grScores {
		fmt.Printf("%s\n", gScore)
	}

}
func getGrsLists(rsIDLines []string) ([]string, map[string]string, map[string]float64, map[string]float64) {
	rsIDList := make([]string, 0, len(rsIDLines))
	eaMap := make(map[string]string)
	eafMap := make(map[string]float64)
	wgtMap := make(map[string]float64)
	colMap := make(map[string]int)

	hdrStr := rsIDLines[0]
	hdrData := strings.Split(hdrStr, ",")
	for i, col := range hdrData {
		colMap[col] = i
	}

	for _, rLine := range rsIDLines[1:] {
		grsData := strings.Split(rLine, ",")
		rsIDList = append(rsIDList, grsData[colMap["varid"]])
		eaMap[grsData[colMap["varid"]]] = grsData[colMap["ea"]]
		eafMap[grsData[colMap["varid"]]], _ = strconv.ParseFloat(grsData[colMap["eaf"]], 64)
		wgtMap[grsData[colMap["varid"]]], _ = strconv.ParseFloat(grsData[colMap["wgt"]], 64)
	}
	return rsIDList, eaMap, eafMap, wgtMap
}
