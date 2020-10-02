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
	"os"
	"strings"
	"variant"
)

// max_posn is intended to be greater than any
// value for genomic position
const max_posn int64 = 9999999999

var empty_record = []string{}

//-----------------------------------------------
// global vars, accessed by multiple funcs
//-----------------------------------------------
var logFilePath string
var vcfPath string
var exclPath string
var threshold float64
var maf float64
var maffactor float64
var hwe float64
var cr float64
var infoscore float64

//-----------------------------------------------
// main package routines
//-----------------------------------------------
func init() {
	const (
		defaultLogFilePath   = "./data/vcffilter_output.log"
		lusage               = "Log file"
		defaultvcfPath       = "./data/test.vcf.gz"
		vusage               = "Full path for vcf files"
		defaultexclPath      = ""
		eusage               = "Full path for excluded posn files"
		defaultThreshold     = 0.9
		thrusage             = "Prob threshold"
		defaultMaf           = 0.05
		mafusage             = "Minor allele frequency"
		defaultMaffact       = 0.0
		maffactusage         = "MAF mult / division factor"
		defaultHwe           = 0.000001
		hweusage             = "Hardy weinberg equilibrium via exact test"
		defaultCr            = 0.9
		crusage              = "Call Rate"
		defaultInfo          = 0.9
		infousage            = "Imputation INFO score"
	)
	flag.StringVar(&logFilePath, "logfile", defaultLogFilePath, lusage)
	flag.StringVar(&logFilePath, "l", defaultLogFilePath, lusage+" (shorthand)")
	flag.StringVar(&vcfPath, "vcfpath", defaultvcfPath, vusage)
	flag.StringVar(&vcfPath, "v", defaultvcfPath, vusage+" (shorthand)")
	flag.StringVar(&exclPath, "exclpath", defaultexclPath, eusage)
	flag.StringVar(&exclPath, "e", defaultexclPath, eusage+" (shorthand)")
	flag.Float64Var(&threshold, "threshold", defaultThreshold, thrusage)
	flag.Float64Var(&threshold, "t", defaultThreshold, thrusage+" (shorthand)")
	flag.Float64Var(&maf, "maf", defaultMaf, mafusage)
	flag.Float64Var(&maf, "m", defaultMaf, mafusage+" (shorthand)")
	flag.Float64Var(&maffactor, "maffactor", defaultMaffact, maffactusage)
	flag.Float64Var(&maffactor, "f", defaultMaffact, maffactusage+" (shorthand)")
	flag.Float64Var(&hwe, "hwe", defaultHwe, hweusage)
	flag.Float64Var(&hwe, "w", defaultHwe, hweusage+" (shorthand)")
	flag.Float64Var(&cr, "cr", defaultCr, crusage)
	flag.Float64Var(&cr, "c", defaultCr, crusage+" (shorthand)")
	flag.Float64Var(&infoscore, "info", defaultInfo, infousage)
	flag.Float64Var(&infoscore, "i", defaultInfo, infousage+" (shorthand)")
	flag.Parse()
}
//------------------------------------------------
// check(error) crude general error handler
//------------------------------------------------
func check(e error) {
	if e != nil {
		log.Println("err != nil")
		log.Fatal(e)
	}
}
//------------------------------------------------
// main() program entry point
//------------------------------------------------
func main() {
	// set up logging to a file
	lf, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0640)
	check(err)
	defer lf.Close()

	log.SetOutput(lf)
	//log.Printf("START filter MAF=%f,CR=%.2f,HWE=%f,INFO=%.2f,threshold=%.2f,maffactor=%.2f\n", maf, cr, hwe, infoscore, threshold, maffactor)
	log.Printf("START filter MAF=%f,CR=%.2f,HWE=%f,INFO=%.2f,threshold=%.2f\n", maf, cr, hwe, infoscore, threshold)
  rcount := 0
  wcount := 0
  pcount := 0
  mcount := 0
  mfcount := 0
  hcount := 0
  crcount := 0
  icount := 0
  dcount := 0
  var excludedPosns = map[string]bool{}
	// Open Excluded posn file if present 
  if exclPath != "" {
    fe, err := os.Open(exclPath)
    check(err)
    defer fe.Close()
    scanner := bufio.NewScanner(fe)
    for scanner.Scan() {
      text := scanner.Text()
      excludedPosns[text] = true
    }
  }
	// Open VCF 
	f, err := os.Open(vcfPath)
	check(err)
	defer f.Close()
	gr, err := gzip.NewReader(f)
	check(err)
	//scanner := bufio.NewScanner(gr)
  reader := bufio.NewReader(gr)
  write_headers(reader)
  found_error := false
  for {
    found_error = false
    rcount += 1
    text, data, _, _, err := get_next_record(reader)
    if err == io.EOF {
      break
    }
    posn := variant.GetPosnStr(data)
    varid := variant.GetVarid(data)
    if _, ok := excludedPosns[posn]; ok {
      pcount +=1
      found_error = true
      //continue
    }
    rec_cr, _, _, rec_maf, rec_hwe, _, _, _, _, _, _, _ := genometrics.Metrics_for_record(data, threshold)
    //log.Printf("VARID=%s\tCR=%f\tMAF=%f\tHWE=%f\tn=%d\n", varid, rec_cr, rec_maf, rec_hwe)

    if rec_maf < maf {
      mcount += 1
      found_error = true
    }
    if rec_hwe < hwe {
      hcount += 1
      found_error = true
    }
    if rec_cr < cr {
      crcount += 1
      found_error = true
    }
    rec_info := variant.GetInfoScore(data)
    if rec_info < infoscore {
      icount += 1
      found_error = true
    }
    if maffactor != 0.0 {
      rec_rpaf := variant.GetRefPanelAF(data)
      if (rec_maf < (rec_rpaf / maffactor)) || (rec_maf > (rec_rpaf * maffactor)) {
        mfcount += 1
        found_error = true
      }
    }
    if varid == "." {
      dcount += 1
      found_error = true
    }
    if found_error == false {
      wcount += 1
      fmt.Printf("%s\n", text)
    }
  }
	log.Printf("END filter Rd=%d, Wrt=%d, posn=%d, maf=%d, hwe=%d, cr=%d, info=%d, maffact=%d, dot=%d\n",
    rcount, wcount, pcount, mcount, hcount, crcount, icount, mfcount, dcount)
}
//-------------------------------------------------------------
// Write headers for the VCF
//-------------------------------------------------------------
func write_headers(rdr *bufio.Reader) {
  for {
    text, err := rdr.ReadString('\n')
    if err == io.EOF {
      return
    }
    text = strings.TrimRight(text, "\n")
    if strings.HasPrefix(text, "##") {
      fmt.Printf("%s\n", text)
    } else {
      if strings.HasPrefix(text, "#") {
        fmt.Printf("%s\n", text)
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
    posn = int64(variant.GetPosn(data))
    varid = variant.GetVarid(data)
  }
  return text, data, posn, varid, err
}

