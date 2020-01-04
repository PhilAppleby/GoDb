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
type DBVariant struct {
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
type DBFilePath struct {
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
// Getvardbdata: get variants collection data for an rsid
// For each result:
//   Find the relevent filepath and return all co-ordinate data
//---------------------------------------------------------------------
func Getvardbdata(session *mgo.Session, gdb string,
                  variant_collection string, filepath_collection string,
                  vcfPathPref string, rsid string) ([]DBVariant, []string) {

	variants := session.DB(gdb).C(variant_collection)
	filepaths := session.DB(gdb).C(filepath_collection)

	dbvariant := DBVariant{}
	var variant_list = make([]DBVariant, 0, 10)
	fdata := DBFilePath{}
  var filepath_list = make([]string, 0, 10)
  variant_count := 0
	// TODO if rsid begins with 'rs' proceed as below
	// set up and execute a query on the variants collection 
	find := variants.Find(bson.M{"rsid": rsid})

	items := find.Iter()
	for items.Next(&dbvariant) {
		// query the filepaths collection
    variant_count += 1
    variant_list = append(variant_list, dbvariant)
		err := filepaths.Find(bson.M{"assaytype": dbvariant.Assaytype}).One(&fdata)
		check(err)
		if len(dbvariant.Chromosome) == 1 {
			dbvariant.Chromosome = "0" + dbvariant.Chromosome
		}
		filestr := fmt.Sprintf("chr%s.vcf.gz", dbvariant.Chromosome)
		fullfilepath := fdata.Filepath + "/" + filestr
		if vcfPathPref != "" {
			fullfilepath = vcfPathPref + "/" + fdata.Fpath_suffix + "/" + filestr
		}
    filepath_list = append(filepath_list, fullfilepath)
	}

  if (variant_count == 0) {
    fmt.Printf("NOT FOUND %s\n", rsid)
  }
  return variant_list, filepath_list
}
//---------------------------------------------------------------------
// Getvarfiledata: Exported function to read a single VCF record using 
// the tabix index for the file.
//---------------------------------------------------------------------
func Getvarfiledata(f string, dbv DBVariant, recs chan string, wg *sync.WaitGroup) {
	Sem <- struct{}{}
	defer func() { <-Sem }()

	tbx, err := bix.New(f)
	fopen_ctr++
	check(err)
	start := dbv.Position - 1
	end := dbv.Position
	if _, ok := sglDigitChrom[dbv.Assaytype]; ok {
		if strings.HasPrefix(dbv.Chromosome, "0") {
			dbv.Chromosome = dbv.Chromosome[1:]
		}
	}
	// Query returns an io.Reader
	rdr, err := tbx.Query(loc{dbv.Chromosome, start, end})
	check(err)
	for {
		v, err := rdr.Next()
		if err != nil {
			break
		}
		vrecord := v.(interfaces.IVariant).String()
		vrecarr := strings.Split(vrecord, "\t")
		recref, recalt := variant.GetAlleles(vrecarr)
		if (recref == dbv.AlleleA) && (recalt == dbv.AlleleB) {
			recs <- fmt.Sprintf("%s\t%s", dbv.Assaytype, vrecord)
		}
	}
	tbx.Close()
	fopen_ctr--
  wg.Done()
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
