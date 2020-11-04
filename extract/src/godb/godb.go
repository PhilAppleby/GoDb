// 
// godb - a collection of methods for GoDb (MongoDb) access
// Including variant and sample retrieval and tabix indexed
// VCF file access
//
package godb

import (
  "encoding/json"
	"fmt"
  "genometrics"
	"github.com/brentp/bix"
	"github.com/brentp/irelate/interfaces"
	"gopkg.in/mgo.v2"
	"gopkg.in/mgo.v2/bson"
	"log"
	"os"
  "sample"
  "sync"
	"strings"
	"variant"
	"vcfmerge"
)

//-----------------------------------------------
// dbconfig struct for db access
//-----------------------------------------------
type dbconfig struct {
  Dbhost          string
  Dbname          string
  Var_collection  string
  Fp_collection   string
  Samp_collection string
}

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
// DBVariant: struct for the mongodb variants collection
//------------------------------------------------------
type DBVariant struct {
  Assaytype string `bson:"assaytype,omitempty"`
  Rsid       string `bson:"rsid,omitempty"`
  AlleleB    string `bson:"alleleB,omitempty"`
  Position   int    `bson:"position,omitempty"`
  AlleleA    string `bson:"alleleA,omitempty"`
  Chromosome string `bson:"chromosome,omitempty"`
  Ref_AF     float64
  Alt_AF     float64
  MAF        float64
  HWEP       float64
  Infoscore  float64
  CR         float64
  Missing    int
	NumSamples int
  Errpct     float64
}
//------------------------------------------------------
// DBFilePath: struct for the mongodb filepaths collection
//------------------------------------------------------
type DBFilePath struct {
	Assaytype    string `bson:"assaytype,omitempty"`
	Filepath     string `bson:"filepath,omitempty"`
	Fpath_prefix string `bson:"fpath_prefix,omitempty"`
	Fpath_suffix string `bson:"fpath_suffix,omitempty"`
}
//------------------------------------------------------
// DBSample: struct for the mongodb samples collection
//------------------------------------------------------
type DBSample struct {
  Assaytype string `bson:"assaytype,omitempty"`
  List_posn int    `bson:"list_posn,omitempty"`
  Sample_id string `bson:"sample_id,omitempty"`
}
//------------------------------------------------------
// DBGeneMap: struct for the mongodb genemap collection
//------------------------------------------------------
type DBGeneMap struct {
	Genename   string `bson:"genename,omitempty"`
	Name       string `bson:"name,omitempty"`
	Chrom      string `bson:"chrom,omitempty"`
	Strand     string `bson:"strand,omitempty"`
	txStart    int    `bson:"txStart,omitempty"`
	txEnd      int    `bson:"txEnd,omitempty"`
}

var fopen_ctr int
var sglDigitChrom map[string]int
var geno_strings = []string{"0/0", "0/1", "1/1"}
var dbconf dbconfig
var session *mgo.Session

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

var fsem chan struct{}

func init() {
  loadDbConfig()
  sess, err := mgo.Dial(dbconf.Dbhost)
  if err != nil {
    log.Fatal("Cannot connect to mongodb at ", dbconf.Dbhost, ", ", err)
  }
  session = sess
  fmt.Printf("Connected to mongodb at %s [%s]\n", dbconf.Dbhost, dbconf.Dbname)
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

func loadDbConfig() {
  config_file := os.Getenv("DBCONFIGFILE")
  file, err := os.Open(config_file)
  if err != nil {
    log.Fatal("Cannot open dbconfig file", err)
  }
  decoder := json.NewDecoder(file)
  dbconf = dbconfig{}
  err = decoder.Decode(&dbconf)
  if err != nil {
    log.Fatalln("Cannot get configuration from file", err)
  }
}

func check(e error) {
	//  fmt.Println("check")
	if e != nil {
		fmt.Println("err != nil")
		log.Fatal(e)
	}
}
//---------------------------------------------------------------------
// Getallvardata: get all var data 
//---------------------------------------------------------------------
func Getallvardata(vcfPathPref string, rsid_list []string, requested_assaytypes map[string]bool, pthr float64) ([]DBVariant) {

  var wg sync.WaitGroup
  // Control # of active goroutines (in this case also the # of open files)
  fsem = make(chan struct{}, 64)
  // 10,000 here is arbitrary, could obtain a count from the Db
  file_records := make(chan string, 10000)
  rsid_count := 0

  for _, rsid := range rsid_list {
    rsid_count += 1
    // For each rsid, access godb and get the lists of variants vs filepaths
    variants, filepaths := Getvardbdata(vcfPathPref, rsid)
    // For each variant, filepath combination get a file record
    for idx, variant := range variants {
      if _, ok := requested_assaytypes[variant.Assaytype]; ok {
        wg.Add(1)
        go getvarfiledata(filepaths[idx], variant, file_records, &wg)
      }
    }
  }

  // Wait for the file-reading go routines (defined in the godb package) to complete
  wg.Wait()
  close(fsem)
  close(file_records)
   // Map rsid's to their retrieved vcf file records
  rsids := make(map[string][][]string, rsid_count)
  // And to their vcf data
  rsids_data := make(map[string][]vcfmerge.Vcfdata, rsid_count)

  var variant_list = make([]DBVariant, 0, 10)

   // What assaytypes do we have in the results?
  assaytypes := make(map[string]bool, 10) // 10 is a guess
  assaytype_list := make([]string, 0)

  // Read the channel of file records
  for record := range file_records {
    var recdata vcfmerge.Vcfdata
    var dbvar DBVariant
    fields := strings.Split(record, "\t")
    // Do we want to output this assaytype?
    if _, ok := requested_assaytypes[fields[0]]; !ok {
      continue
    }
    if requested_assaytypes[fields[0]] != true {
      continue
    }
    prfx, _ := variant.GetVCFPrfx_Sfx(fields[1:])
    dbvar.Assaytype = fields[0]
    dbvar.Rsid = variant.GetVarid(prfx)
    dbvar.Chromosome = variant.GetChrom(prfx)
    dbvar.Position = variant.GetPosn(prfx)
    dbvar.AlleleA, dbvar.AlleleB = variant.GetAlleles(prfx)
    dbvar.CR, dbvar.Ref_AF, dbvar.Alt_AF, dbvar.MAF, dbvar.HWEP, _, _, _, dbvar.NumSamples, dbvar.Missing, _, _ = genometrics.Metrics_for_record(fields[1:], pthr)
    dbvar.Infoscore = variant.GetInfoScore(prfx)
    recdata.Probidx = variant.GetProbidx(prfx)
    // for determining combined rec size
    if _, ok := assaytypes[fields[0]]; !ok {
      assaytypes[fields[0]] = true
      assaytype_list = append(assaytype_list, fields[0])
    }
    if _, ok := rsids[dbvar.Rsid]; !ok {
      rsids[dbvar.Rsid] = make([][]string, 0)
      rsids_data[dbvar.Rsid] = make([]vcfmerge.Vcfdata, 0)
    }
    rsids[dbvar.Rsid] = append(rsids[dbvar.Rsid], fields)
    rsids_data[dbvar.Rsid] = append(rsids_data[dbvar.Rsid], recdata)
    variant_list = append(variant_list, dbvar)
  }
  // get all sample data from godb and organise into maps of maps:
  // assaytype -> sample name -> sample posn (sample_name_map)
  // assaytype -> sample posn -> sample name (sample_posn_map)
  sample_name_map, sample_posn_map := GetSamplesByAssaytype()
  // Condense all sample_names into a combined map samplename -> record position
  combocols := sample.GetCombinedSampleMapByAssaytypes(sample_name_map, assaytype_list)
  _, combo_names := vcfmerge.GetCombinedColumnHeaders(combocols)
  // Get column headers as a single tab delimited string, with prefix in place, and as a list, both in postion order

  // output the vcf records in input order, can also log the 'NOT FOUND's at this point
  for _, rsid := range rsid_list {
    if records, ok := rsids[rsid]; ok {
      var rsid_genomet genometrics.AllMetrics
      var dbvar DBVariant
      rec_str := vcfmerge.Combine_one(records, rsids_data[rsid], rsid, sample_posn_map, combocols, combo_names, pthr, &rsid_genomet)
      fields := strings.Split(rec_str, "\t")
      prfx, sfx := variant.GetVCFPrfx_Sfx(fields)
      //log.Printf("SFX = %v\n", sfx)
      dot := 0
      dbvar.Assaytype = "combined"
      dbvar.Rsid = variant.GetVarid(prfx)
      dbvar.Chromosome = variant.GetChrom(prfx)
      dbvar.Position = variant.GetPosn(prfx)
      dbvar.AlleleA, dbvar.AlleleB = variant.GetAlleles(prfx)
      dbvar.CR, dbvar.Ref_AF, dbvar.Alt_AF, dbvar.MAF, dbvar.HWEP, _, _, _, dbvar.NumSamples, dbvar.Missing, dot, _ = genometrics.Metrics_for_record(fields, pthr)
      dbvar.Infoscore = variant.GetInfoScore(prfx)
      dbvar.Errpct = (float64(rsid_genomet.MismatchCount) / float64(rsid_genomet.OverlapTestCount)) * 100
      log.Printf("%s combined, SFX len = %d, samples=%d, miss=%d, dot=%d\n", rsid, len(sfx), dbvar.NumSamples, dbvar.Missing, dot)
      //log.Printf("%s, SFX = %v\n", rsid, sfx)
      variant_list = append(variant_list, dbvar)
    }
  }
  return variant_list
}
//------------------------------------------------------------------------------
// wrap the godb.Getvarfiledata func, for use as a goroutine
//------------------------------------------------------------------------------
func getvarfiledata(f string, dbv DBVariant, recs chan string, wg *sync.WaitGroup) {
  fsem <- struct{}{}
  defer func() { <-fsem }()
  Getvarfiledata(f, dbv, recs)
  wg.Done()
}
//---------------------------------------------------------------------
// Getvardbdata: get variants collection data for an rsid
// For each result:
//   Find the relevent filepath and return all co-ordinate data
//---------------------------------------------------------------------
func Getvardbdata(vcfPathPref string, rsid string) ([]DBVariant, []string) {

	variants := session.DB(dbconf.Dbname).C(dbconf.Var_collection)
	filepaths := session.DB(dbconf.Dbname).C(dbconf.Fp_collection)

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
    log.Printf("##NOT FOUND %s\n", rsid)
  }
  return variant_list, filepath_list
}
//---------------------------------------------------------------------
// Getvarfiledata: Exported function to read a single VCF record using 
// the tabix index for the file.
//---------------------------------------------------------------------
func Getvarfiledata(f string, dbv DBVariant, recs chan string) {
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
}
//---------------------------------------------------------------------
// GatSamplesByAssaytype: Get all samplea, for all AssayTtypes
// and their array indexes in the relevant VCF data
// NOTE:  this returns two maps of map: 
//   assaytype to sample_name to index 
//   assaytype to index to sample_name 
//---------------------------------------------------------------------
func GetSamplesByAssaytype() (map[string]map[string]int, map[string]map[int]string) {
  sampleNamePosn := make(map[string]map[string]int)
  samplePosnName := make(map[string]map[int]string)

  samp := session.DB(dbconf.Dbname).C(dbconf.Samp_collection)

  sample := DBSample{}

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
