// Package godb ...
// A collection of methods for GoDb (MongoDb) access
// Including variant and sample retrieval and tabix indexed
// VCF file access
//
package godb

import (
	"encoding/json"
	"fmt"
	"genometrics"
	"log"
	"os"
	"sample"
	"strings"
	"sync"
	"variant"
	"vcfmerge"

	"github.com/brentp/bix"
	"github.com/brentp/irelate/interfaces"
	"gopkg.in/mgo.v2"
	"gopkg.in/mgo.v2/bson"
)

//-----------------------------------------------
// dbconfig struct for db access
//-----------------------------------------------
type dbconfig struct {
	Dbhost         string
	Dbname         string
	VarCollection  string
	FpCollection   string
	SampCollection string
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

// DBVariant ...
// struct for the mongodb variants collection
type DBVariant struct {
	Assaytype     string `bson:"assaytype,omitempty"`
	Rsid          string `bson:"rsid,omitempty"`
	AlleleB       string `bson:"alleleB,omitempty"`
	StartPosition int    `bson:"position,omitempty"`
	EndPosition   int
	AlleleA       string `bson:"alleleA,omitempty"`
	Chromosome    string `bson:"chromosome,omitempty"`
	RefAF         float64
	AltAF         float64
	MAF           float64
	HWEP          float64
	Infoscore     float64
	CR            float64
	Missing       int
	NumSamples    int
	Errpct        float64
	LineNum       int
}

// DBFilePath ...
// struct for the mongodb filepaths collection
type DBFilePath struct {
	Assaytype   string `bson:"assaytype,omitempty"`
	Filepath    string `bson:"filepath,omitempty"`
	FpathPrefix string `bson:"fpath_prefix,omitempty"`
	FpathSuffix string `bson:"fpath_suffix,omitempty"`
}

// DBSample ...
// struct for the mongodb samples collection
type DBSample struct {
	Assaytype string `bson:"assaytype,omitempty"`
	ListPosn  int    `bson:"list_posn,omitempty"`
	SampleID  string `bson:"sample_id,omitempty"`
}

// DBGeneMap ...
// struct for the mongodb genemap collection
type DBGeneMap struct {
	Genename string `bson:"genename,omitempty"`
	Name     string `bson:"name,omitempty"`
	Chrom    string `bson:"chrom,omitempty"`
	Strand   string `bson:"strand,omitempty"`
	txStart  int    `bson:"txStart,omitempty"`
	txEnd    int    `bson:"txEnd,omitempty"`
}

var fopenCtr int
var sglDigitChrom map[string]int
var genoStrings = []string{"0/0", "0/1", "1/1"}
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
	log.Printf("Godb connected to mongodb at %s [%s]\n", dbconf.Dbhost, dbconf.Dbname)
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
	configFile := os.Getenv("DBCONFIGFILE")
	file, err := os.Open(configFile)
	check("godb: Cannot open dbconfig file", err)

	decoder := json.NewDecoder(file)
	dbconf = dbconfig{}
	err = decoder.Decode(&dbconf)
	check("Cannot read dbconfig file", err)
}

func check(msg string, e error) {
	if e != nil {
		log.Println(msg)
		log.Fatal(e)
	}
}

// Getallvardata ...
// get all variant and geno data for a list of variants (rsids)
// NOTE: this function uses goroutines for parallel access to file resources
//---------------------------------------------------------------------
func Getallvardata(vcfPathPref string, rsidList []string, requestedAssaytypes map[string]bool, pthr float64) ([]DBVariant, []DBVariant, []string) {

	var wg sync.WaitGroup
	// Control # of active goroutines (in this case also the # of open files)
	fsem = make(chan struct{}, 64)
	// 10,000 here is arbitrary, could obtain a count from the Db
	fileRecords := make(chan string, 10000)
	rsidCount := 0

	for _, rsid := range rsidList {
		rsidCount++
		// For each rsid, access godb and get the lists of variants vs filepaths
		variants, filepaths := Getvardbdata(vcfPathPref, rsid)
		// For each variant, filepath combination get a file record
		for idx, variant := range variants {
			if _, ok := requestedAssaytypes[variant.Assaytype]; ok {
				wg.Add(1)
				go getvarfiledata(filepaths[idx], variant, fileRecords, &wg)
			}
		}
	}

	// Wait for the file-reading go routines (defined in the godb package) to complete
	wg.Wait()
	close(fsem)
	close(fileRecords)
	// Map rsid's to their retrieved vcf file records
	rsids := make(map[string][][]string, rsidCount)
	// And to their vcf data
	rsidsData := make(map[string][]vcfmerge.Vcfdata, rsidCount)

	var variantList = make([]DBVariant, 0, 10)
	var combinedVariantList = make([]DBVariant, 0, 10)
	var combinedRecords = make([]string, 0, 10)

	// What assaytypes do we have in the results?
	assaytypes := make(map[string]bool, 10) // 10 is a guess
	assaytypeList := make([]string, 0)

	lineCount := 0

	// Read the channel of file records
	for record := range fileRecords {
		var recdata vcfmerge.Vcfdata
		var dbvar DBVariant
		lineCount++
		fields := strings.Split(record, "\t")
		// Do we want to output this assaytype?
		if _, ok := requestedAssaytypes[fields[0]]; !ok {
			continue
		}
		if requestedAssaytypes[fields[0]] != true {
			continue
		}
		prfx, _ := variant.GetVCFPrfxSfx(fields[1:])
		dbvar.Assaytype = fields[0]
		dbvar.Rsid = variant.GetVarid(prfx)
		dbvar.Chromosome = variant.GetChrom(prfx)
		dbvar.StartPosition = variant.GetPosn(prfx)
		dbvar.EndPosition = dbvar.StartPosition
		dbvar.AlleleA, dbvar.AlleleB = variant.GetAlleles(prfx)
		dbvar.CR, dbvar.RefAF, dbvar.AltAF, dbvar.MAF, dbvar.HWEP, _, _, _, dbvar.NumSamples, dbvar.Missing, _, _ = genometrics.MetricsForRecord(fields[1:], pthr)
		dbvar.Infoscore = variant.GetInfoScore(prfx)
		dbvar.LineNum = lineCount
		recdata.Probidx = variant.GetProbIdx(prfx)
		// for determining combined rec size
		if _, ok := assaytypes[fields[0]]; !ok {
			assaytypes[fields[0]] = true
			assaytypeList = append(assaytypeList, fields[0])
		}
		if _, ok := rsids[dbvar.Rsid]; !ok {
			rsids[dbvar.Rsid] = make([][]string, 0)
			rsidsData[dbvar.Rsid] = make([]vcfmerge.Vcfdata, 0)
		}
		rsids[dbvar.Rsid] = append(rsids[dbvar.Rsid], fields)
		rsidsData[dbvar.Rsid] = append(rsidsData[dbvar.Rsid], recdata)
		variantList = append(variantList, dbvar)
	}
	// get all sample data from godb and organise into maps of maps:
	// assaytype -> sample name -> sample posn (sampleNameMap)
	// assaytype -> sample posn -> sample name (samplePosnMap)
	sampleNameMap, samplePosnMap := GetSamplesByAssaytype()
	// Condense all sample_names into a combined map samplename -> record position
	combocols := sample.GetCombinedSampleMapByAssaytypes(sampleNameMap, assaytypeList)
	// Get column headers as a single tab delimited string, with prefix in place, and as a list, both in postion order
	comboStr, comboNames := vcfmerge.GetCombinedColumnHeaders(combocols)
	combinedRecords = append(combinedRecords, comboStr)

	// output the vcf records in input order, can also log the 'NOT FOUND's at this point
	lineCount = 0
	for _, rsid := range rsidList {
		if records, ok := rsids[rsid]; ok {
			var rsidGenomet genometrics.AllMetrics
			var dbvar DBVariant
			recStr := vcfmerge.CombineOne(records, rsidsData[rsid], rsid, samplePosnMap, combocols, comboNames, pthr, &rsidGenomet)
			combinedRecords = append(combinedRecords, recStr)
			fields := strings.Split(recStr, "\t")
			prfx, _ := variant.GetVCFPrfxSfx(fields)
			//log.Printf("SFX = %v\n", sfx)
			//dot := 0
			lineCount++
			dbvar.Assaytype = "combined"
			dbvar.Rsid = variant.GetVarid(prfx)
			dbvar.Chromosome = variant.GetChrom(prfx)
			dbvar.StartPosition = variant.GetPosn(prfx)
			dbvar.EndPosition = dbvar.StartPosition
			dbvar.AlleleA, dbvar.AlleleB = variant.GetAlleles(prfx)
			dbvar.CR, dbvar.RefAF, dbvar.AltAF, dbvar.MAF, dbvar.HWEP, _, _, _, dbvar.NumSamples, dbvar.Missing, _, _ = genometrics.MetricsForRecord(fields, pthr)
			dbvar.Infoscore = variant.GetInfoScore(prfx)
			dbvar.Errpct = (float64(rsidGenomet.MismatchCount) / float64(rsidGenomet.OverlapTestCount)) * 100.0
			dbvar.LineNum = lineCount
			//log.Printf("%s combined, SFX len = %d, samples=%d, miss=%d, dot=%d\n", rsid, len(sfx), dbvar.NumSamples, dbvar.Missing, dot)
			log.Printf("%s combined, mismatch=%d, overlaps=%d, ErrPct=%.5f\n", rsid, rsidGenomet.MismatchCount, rsidGenomet.OverlapTestCount, dbvar.Errpct)

			//log.Printf("%s, SFX = %v\n", rsid, sfx)
			combinedVariantList = append(combinedVariantList, dbvar)
		}
	}
	return variantList, combinedVariantList, combinedRecords
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

// Getvardbdata ...
// get variants collection data for an rsid
// For each result:
//   Find the relevent filepath and return all co-ordinate data
func Getvardbdata(vcfPathPref string, rsid string) ([]DBVariant, []string) {

	variants := session.DB(dbconf.Dbname).C(dbconf.VarCollection)
	filepaths := session.DB(dbconf.Dbname).C(dbconf.FpCollection)

	dbvariant := DBVariant{}
	var variantList = make([]DBVariant, 0, 10)
	fdata := DBFilePath{}
	var filepathList = make([]string, 0, 10)
	variantCount := 0
	// TODO if rsid begins with 'rs' proceed as below
	// set up and execute a query on the variants collection
	//log.Printf("##SEARCH %s\n", rsid)
	find := variants.Find(bson.M{"rsid": rsid})

	items := find.Iter()
	for items.Next(&dbvariant) {
		// query the filepaths collection
		dbvariant.EndPosition = dbvariant.StartPosition
		variantCount++
		variantList = append(variantList, dbvariant)
		err := filepaths.Find(bson.M{"assaytype": dbvariant.Assaytype}).One(&fdata)
		check("Filepaths dbaccess error", err)
		if len(dbvariant.Chromosome) == 1 {
			dbvariant.Chromosome = "0" + dbvariant.Chromosome
		}
		filestr := fmt.Sprintf("chr%s.vcf.gz", dbvariant.Chromosome)
		fullfilepath := fdata.Filepath + "/" + filestr
		if vcfPathPref != "" {
			fullfilepath = vcfPathPref + "/" + fdata.FpathSuffix + "/" + filestr
		}
		filepathList = append(filepathList, fullfilepath)
	}

	if variantCount == 0 {
		log.Printf("##NOT FOUND %s\n", rsid)
	}
	return variantList, filepathList
}

// Getvarfiledata ...
// Exported function to read a single VCF record using
// the tabix index for the file.
func Getvarfiledata(f string, dbv DBVariant, recs chan string) {
	tbx, err := bix.New(f)
	fopenCtr++
	check("Tabix file open error: "+f, err)
	start := dbv.StartPosition - 1
	end := dbv.StartPosition
	if _, ok := sglDigitChrom[dbv.Assaytype]; ok {
		if strings.HasPrefix(dbv.Chromosome, "0") {
			dbv.Chromosome = dbv.Chromosome[1:]
		}
	}
	// Query returns an interfaces.RelatableIterator
	rdr, err := tbx.Query(loc{dbv.Chromosome, start, end})
	check("Tabix file read error: "+f, err)
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
	fopenCtr--
}

// GetvarfiledataByRange ...
// Exported function to find and return an io Reader over a range of records
// File access is by Tabix index
func GetvarfiledataByRange(f string, dbv DBVariant) interfaces.RelatableIterator {
	tbx, err := bix.New(f)
	check("Tabix file open error: "+f, err)
	start := dbv.StartPosition - 1
	end := dbv.StartPosition
	if _, ok := sglDigitChrom[dbv.Assaytype]; ok {
		if strings.HasPrefix(dbv.Chromosome, "0") {
			dbv.Chromosome = dbv.Chromosome[1:]
		}
	}
	// Query returns an io.Reader
	rdr, err := tbx.Query(loc{dbv.Chromosome, start, end})
	check("Tabix file read error: "+f, err)
	return rdr
}

// GetSamplesByAssaytype ...
// Get all samplea, for all AssayTypes
// and their array indexes in the relevant VCF data
// NOTE:  this returns two maps of map:
//   assaytype to sample_name to index
//   assaytype to index to sample_name
func GetSamplesByAssaytype() (map[string]map[string]int, map[string]map[int]string) {
	sampleNamePosn := make(map[string]map[string]int)
	samplePosnName := make(map[string]map[int]string)

	samp := session.DB(dbconf.Dbname).C(dbconf.SampCollection)

	sample := DBSample{}

	find := samp.Find(bson.M{})

	items := find.Iter()
	for items.Next(&sample) {
		if _, ok := sampleNamePosn[sample.Assaytype]; !ok {
			sampleNamePosn[sample.Assaytype] = make(map[string]int)
			samplePosnName[sample.Assaytype] = make(map[int]string)
		}
		sampleNamePosn[sample.Assaytype][sample.SampleID] = sample.ListPosn
		samplePosnName[sample.Assaytype][sample.ListPosn] = sample.SampleID
	}
	return sampleNamePosn, samplePosnName
}

// FormatOutput ...
// Format an array of VCF lines for Output, assume "", "vcf" or "csv" for "option"
//func FormatOutput(records []string, option string) (output string, outFmt string) {
//	if option == "" || option == "vcf" {
//		return strings.Join(records, "\n"), "vcf"
//	}
// first entry in records is vcf samples header
//	_, sampleCols := variant.GetVCFPrfxSfx(strings.Split(records[0], "\t"))
//	sampleMap := make(map[string][]int)

//	for _, sample := range sampleCols {
//		sampleMap[sample] = make([]int, 0, len(records))
//	}

//	for _, record := range records[1:] {
//		_, genos := variant.GetVCFPrfxSfx(strings.Split(record, "\t"))
//		for _, geno := range genos {
//
//		}
//	}

//	return "", option
//}
