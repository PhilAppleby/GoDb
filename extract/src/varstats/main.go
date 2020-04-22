//--------------------------------------------------------------------------------------
// Filter vcf files based on common attributes
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
	"math"
	"os"
	"strings"
	"variant"
)

// max_posn is intended to be greater than any
// value for genomic position
const max_posn int64 = math.MaxInt64
var empty_record = []string{}
//-----------------------------------------------
// global vars, accessed by multiple funcs
//-----------------------------------------------
var logFilePath string
var vcfPath string
var threshold float64
//-----------------------------------------------
// main package routines
//-----------------------------------------------
func init() {
	const (
		defaultLogFilePath   = "./data/vcffilter_output.log"
		lusage               = "Log file"
		defaultvcfPath       = "./data/test.vcf.gz"
		vusage               = "Full path for vcf files"
		defaultThreshold     = 0.9
		thrusage             = "Prob threshold"
	)
	flag.StringVar(&logFilePath, "logfile", defaultLogFilePath, lusage)
	flag.StringVar(&logFilePath, "l", defaultLogFilePath, lusage+" (shorthand)")
	flag.StringVar(&vcfPath, "vcfpath", defaultvcfPath, vusage)
	flag.StringVar(&vcfPath, "v", defaultvcfPath, vusage+" (shorthand)")
	flag.Float64Var(&threshold, "threshold", defaultThreshold, thrusage)
	flag.Float64Var(&threshold, "t", defaultThreshold, thrusage+" (shorthand)")
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
	lf, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0640)
	check(err)
	defer lf.Close()

	log.SetOutput(lf)
	log.Printf("START pthr=%.2f\n", threshold)
  rcount := 0
  wcount := 0
	// Open VCF 
	f, err := os.Open(vcfPath)
	check(err)
	defer f.Close()
	gr, err := gzip.NewReader(f)
	check(err)
	//scanner := bufio.NewScanner(gr)
  reader := bufio.NewReader(gr)
  fmt.Printf("chr,rsid,posn,CR,MAF,HWEP,INFO,N,MISS\n")
  skip_headers(reader)
  for {
    rcount += 1
    _, data, _, _, err := get_next_record(reader)
    if err == io.EOF {
      break
    }
    chrom := variant.GetChrom(data)
    posn := variant.GetPosnStr(data)
    varid := variant.GetVarid(data)
    rec_cr, _, _, rec_maf, rec_hwe, _, _, _, n, miss, _, _ := genometrics.Metrics_for_record(data, threshold)
    rec_info := variant.GetInfoScore(data)
    wcount += 1
    fmt.Printf("%s,%s,%s,%.2f,%.6f,%.8f,%.6f,%d,%d\n", chrom, posn, varid, rec_cr, rec_maf, rec_hwe, rec_info, n, miss)
  }
	log.Printf("END filter Rd=%d, Wrt=%d\n", rcount, wcount)
}
//-------------------------------------------------------------
// skip VCF headers
//-------------------------------------------------------------
func skip_headers(rdr *bufio.Reader) {
  for {
    text, err := rdr.ReadString('\n')
    if err == io.EOF {
      return
    }
    text = strings.TrimRight(text, "\n")
    if strings.HasPrefix(text, "##") {
      continue
    } else {
      if strings.HasPrefix(text, "#") {
        break
      }
    }
  }
  return 
}
//-------------------------------------------------------------
// Read a record from a single reader and split it to format a string array
//-------------------------------------------------------------
func get_next_record(rdr *bufio.Reader) (string, []string, int64, string, error) {
  posn := max_posn
  data := empty_record
  varid := ""
  text, err := rdr.ReadString('\n')
  if err == nil {
    text = strings.TrimRight(text, "\n")
    data = strings.Split(text, "\t")
    posn = variant.GetPosn(data)
    varid = variant.GetVarid(data)
  }
  return text, data, posn, varid, err
}

