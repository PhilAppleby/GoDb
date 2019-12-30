// 
// godb - a collection of methods for GoDb (MongoDb) access
// Including variant and sample retrieval and tabix indexed
// VCF file access
//
package godb

import (
	"fmt"
	"github.com/brentp/bix"
	"github.com/brentp/irelate/interfaces"
	"gopkg.in/mgo.v2"
	"gopkg.in/mgo.v2/bson"
	"log"
	"strings"
	"sync"
	"variant"
)

//-----------------------------------------------
// loc - struct and methods for tabix file access
//-----------------------------------------------
type loc struct {
	chrom string
	start int
	end   int
}

func (s loc) Chrom() string {
	return s.chrom
}
func (s loc) Start() uint32 {
	return uint32(s.start)
}
func (s loc) End() uint32 {
	return uint32(s.end)
}

//------------------------------------------------------
// Variant: struct for the mongodb variants collection
//------------------------------------------------------
type Variant struct {
	//ID           bson.ObjectId `bson:"_id,omitempty"`
	Assaytype string `bson:"assaytype,omitempty"`
	//Info         float64       `bson:"info,omitempty"`
	//Cohort_1_hwe float64       `bson:"cohort_1_hwe,omitempty"`
	Rsid       string `bson:"rsid,omitempty"`
	AlleleB    string `bson:"alleleB,omitempty"`
	Position   int    `bson:"position,omitempty"`
	AlleleA    string `bson:"alleleA,omitempty"`
	Chromosome string `bson:"chromosome,omitempty"`
	//All_maf      float64       `bson:"all_maf,omitempty"`
}
//------------------------------------------------------
// FilePath: struct for the mongodb filepaths collection
//------------------------------------------------------
type FilePath struct {
	Assaytype    string `bson:"assaytype,omitempty"`
	Filepath     string `bson:"filepath,omitempty"`
	Fpath_prefix string `bson:"fpath_prefix,omitempty"`
	Fpath_suffix string `bson:"fpath_suffix,omitempty"`
}
//------------------------------------------------------
// Sample: struct for the mongodb samples collection
//------------------------------------------------------
type Sample struct {
  Assaytype string `bson:"assaytype,omitempty"`
  List_posn int    `bson:"list_posn,omitempty"`
  Sample_id string `bson:"sample_id,omitempty"`
}

var Sem chan struct{}
var Wg sync.WaitGroup
var fopen_ctr int
var sglDigitChrom map[string]int
var geno_strings = []string{"0/0", "0/1", "1/1"}

const firstGenoIdx = 9
const chrIdx = 0
const posnIdx = 1
const varIdx = 2
const refIdx = 3
const altIdx = 4
const qcIdx = 5
const filtIdx = 6
const infoIdx = 7
const fmtIdx = 8

func init() {
	sglDigitChrom = make(map[string]int)
	sglDigitChrom["atest"] = 1
	sglDigitChrom["atest2"] = 1
	sglDigitChrom["atest3"] = 1
	sglDigitChrom["atest4"] = 1
	sglDigitChrom["atest5"] = 1
	sglDigitChrom["affy"] = 1
	sglDigitChrom["illumina"] = 1
	sglDigitChrom["broad"] = 1
}

func check(e error) {
	//  fmt.Println("check")
	if e != nil {
		fmt.Println("err != nil")
		log.Fatal(e)
	}
}
//---------------------------------------------------------------------
// Getvardata: get variants collection data for an rsid
// For each result:
//   Find the relevent filepath, access the file and place the record 
//   in the 'recs' channel
//---------------------------------------------------------------------
func Getvardata(par string, session *mgo.Session, gdb string,
                variant_collection string, filepath_collection string,
                vcfPathPref string, rsid string, validAssaytypes map[string]bool, recs chan string) {

	variants := session.DB(gdb).C(variant_collection)
	filepaths := session.DB(gdb).C(filepath_collection)

	variant := Variant{}
	fdata := FilePath{}
  variant_count := 0
	// TODO if rsid begins with 'rs' proceed as below
	// else if rsid is of the form chr:posn (eg 2:174000523)
	// make a range query

	// set up and execute a query on the variants collection 
	find := variants.Find(bson.M{"rsid": rsid})

	items := find.Iter()
	for items.Next(&variant) {
		// query the filepaths collection
    variant_count += 1
    if _, ok := validAssaytypes[variant.Assaytype]; !ok {
      continue
    }
    if validAssaytypes[variant.Assaytype] != true {
      continue
    }
		err := filepaths.Find(bson.M{"assaytype": variant.Assaytype}).One(&fdata)
		check(err)
		if len(variant.Chromosome) == 1 {
			variant.Chromosome = "0" + variant.Chromosome
		}
		filestr := fmt.Sprintf("chr%s.vcf.gz", variant.Chromosome)
		fullfilepath := fdata.Filepath + "/" + filestr
		if vcfPathPref != "" {
			fullfilepath = vcfPathPref + "/" + fdata.Fpath_suffix + "/" + filestr
		}
		if par == "Y" {
			// This is where a go routing is started to handle I/O in parallel
			go accessvcffile(fullfilepath, variant.Assaytype, variant.Chromosome, variant.Position, variant.AlleleA, variant.AlleleB, recs)
		} else {
			accessvcffile(fullfilepath, variant.Assaytype, variant.Chromosome, variant.Position, variant.AlleleA, variant.AlleleB, recs)
		}
	}
  if (variant_count == 0) {
    fmt.Printf("NOT FOUND %s\n", rsid)
  }
}
//---------------------------------------------------------------------
// accessvcffile: Internal function to read a single VCF record using 
// the tabix index for the file.
//---------------------------------------------------------------------
func accessvcffile(f string, assaytype string, chrom string, posn int, ref string, alt string, recs chan string) {
	Wg.Add(1)
	defer Wg.Done()
	Sem <- struct{}{}
	defer func() { <-Sem }()
	tbx, err := bix.New(f)
	fopen_ctr++
	check(err)
	start := posn - 1
	end := posn
	if _, ok := sglDigitChrom[assaytype]; ok {
		if strings.HasPrefix(chrom, "0") {
			chrom = chrom[1:]
		}
	}
	// Query returns an io.Reader
	rdr, err := tbx.Query(loc{chrom, start, end})
	check(err)
	for {
		v, err := rdr.Next()
		if err != nil {
			break
		}
		vrecord := v.(interfaces.IVariant).String()
		vrecarr := strings.Split(vrecord, "\t")
		recref, recalt := variant.GetAlleles(vrecarr)
		if (recref == ref) && (recalt == alt) {
			recs <- fmt.Sprintf("%s\t%s", assaytype, vrecord)
		}
	}
	tbx.Close()
	fopen_ctr--
}
//---------------------------------------------------------------------
// GatSamplesByAssaytype: Get all samplea, for all AssayTtypes
// and their array indexes in the relevant VCF data
// NOTE:  this returns two maps of map: 
//   assaytype to sample_name to index 
//   assaytype to inddex to sample_name 
//---------------------------------------------------------------------
func GetSamplesByAssaytype(session *mgo.Session, gdb string, sample_collection string) (map[string]map[string]int, map[string]map[int]string) {
  sampleNamePosn := make(map[string]map[string]int)
  samplePosnName := make(map[string]map[int]string)

  samp := session.DB(gdb).C(sample_collection)

  sample := Sample{}

  find := samp.Find(bson.M{})

  items := find.Iter()
  for items.Next(&sample) {
    if _, ok := sampleNamePosn[sample.Assaytype]; !ok {
      sampleNamePosn[sample.Assaytype] = make(map[string]int)
      samplePosnName[sample.Assaytype] = make(map[int]string)
    }
    sampleNamePosn[sample.Assaytype][sample.Sample_id] = sample.List_posn
    samplePosnName[sample.Assaytype][sample.List_posn] = sample.Sample_id
  }
  return sampleNamePosn, samplePosnName
}
