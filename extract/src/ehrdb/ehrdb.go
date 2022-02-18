//
// ehrdb - a collection of methods for EhrDb (MongoDb) access
// - add data to and retrieve from pheno, pheno_meta
//
package ehrdb

import (
	"encoding/json"
	"fmt"
	"log"
	"os"

	"gopkg.in/mgo.v2"
	"gopkg.in/mgo.v2/bson"
)

//-----------------------------------------------
// dbconfig struct for db access
//-----------------------------------------------
type dbconfig struct {
	Dbhost                string
	Dbname                string
	PhenoCollection       string
	PhenoMetaCollection   string
	GrsInputCollection    string
	GrsMetaCollection     string
	GrsScoreCollection    string
	VarlistMetaCollection string
	VarlistCollection     string
}

// DBPheno ...
// struct for the mongodb pheno collection
//------------------------------------------------------
type DBPheno struct {
	IId   string `bson:"iid,omitempty"`
	Name  string `bson:"name,omitempty"`
	Value string `bson:"value,omitempty"`
}

// DBPhenoMeta ...
// struct for the mongodb pheno meta collection
//------------------------------------------------------
type DBPhenoMeta struct {
	Name        string `bson:"name,omitempty"`
	Source      string `bson:"source,omitempty"`
	Description string `bson:"description,omitempty"`
	PhenoClass  string `bson:"phenoclass,omitempty"`
}

// DBGrsInput ...
// struct for the mongodb grsinput collection
//------------------------------------------------------
type DBGrsInput struct {
	Name   string `bson:"name,omitempty"`
	Varid  string `bson:"varid,omitempty"`
	Ea     string `bson:"ea,omitempty"`
	Eaf    string `bson:"eaf,omitempty"`
	Weight string `bson:"wgt,omitempty"`
}

// DBGrsMeta ...
// struct for the mongodb grsinput collection
//------------------------------------------------------
type DBGrsMeta struct {
	Name        string `bson:"name,omitempty"`
	Description string `bson:"description,omitempty"`
}

// DBGrsScore ...
// struct for the mongodb grsscore collection
//------------------------------------------------------
type DBGrsScore struct {
	Name     string `bson:"name,omitempty"`
	IID      string `bson:"iid,omitempty"`
	SnpCount int    `bson:"snpcount,omitempty"`
	Score    string `bson:"score,omitempty"`
}

// DBVarlistMeta ...
// struct for the mongodb variantlist meta collection
//------------------------------------------------------
type DBVarlistMeta struct {
	Name        string `bson:"name,omitempty"`
	Description string `bson:"description,omitempty"`
}

// DBListVar ...
// struct for the mongodb variant in a list of variants collection
//------------------------------------------------------
type DBListVar struct {
	Name  string `bson:"name,omitempty"`
	VarID string `bson:"varid,omitempty"`
}

var dbconf dbconfig
var session *mgo.Session

func init() {
	loadDbConfig()
	sess, err := mgo.Dial(dbconf.Dbhost)
	if err != nil {
		log.Fatal("Cannot connect to mongodb at ", dbconf.Dbhost, ", ", err)
	}
	session = sess
	log.Printf("Ehrdb connected to mongodb at %s [%s]\n", dbconf.Dbhost, dbconf.Dbname)
}

func loadDbConfig() {
	configFile := os.Getenv("DBEHRCONFIGFILE")
	file, err := os.Open(configFile)
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

// ******* Stored phenotype section ***********************************

// GetPhenoByName ...
// Get all entries for a phenotype by name
//---------------------------------------------------------------------
func GetPhenoByName(phenoname string) (map[string]string, int) {
	phenoIDValue := make(map[string]string)

	phenoColl := session.DB(dbconf.Dbname).C(dbconf.PhenoCollection)

	phenotype := DBPheno{}

	find := phenoColl.Find(bson.M{"name": phenoname})

	items := find.Iter()
	count := 0
	for items.Next(&phenotype) {
		phenoIDValue[phenotype.IId] = phenotype.Value
		count++
	}
	return phenoIDValue, count
}

// GetPhenoMetaByName ...
// Get all entries for a phenotype by name
//---------------------------------------------------------------------
func GetPhenoMetaByName(phenoname string) DBPhenoMeta {
	phenoMetaColl := session.DB(dbconf.Dbname).C(dbconf.PhenoMetaCollection)

	phenoMeta := DBPhenoMeta{}

	find := phenoMetaColl.Find(bson.M{"name": phenoname})

	items := find.Iter()
	for items.Next(&phenoMeta) {
	}
	return phenoMeta
}

// GetPhenoMetaNames ...
// Get all pheno names from the pheno meta collection
//---------------------------------------------------------------------
func GetPhenoMetaNames() []string {
	var phenoNameList = make([]string, 0, 10)

	phenoMetaColl := session.DB(dbconf.Dbname).C(dbconf.PhenoMetaCollection)

	phenoMeta := DBPhenoMeta{}

	find := phenoMetaColl.Find(bson.M{})

	items := find.Iter()
	for items.Next(&phenoMeta) {
		phenoNameList = append(phenoNameList, phenoMeta.Name)
	}
	return phenoNameList
}

// InsertPhenoDataWithCheck ...
// Insert pheno data from a list and a phenoMeta record
// but first check for existence (in the phenoMetaColl)
//---------------------------------------------------------------------
func InsertPhenoDataWithCheck(phenoname string, phenosource string, phenodesc string,
	phenoclass string, phenoitems map[string]string) (bool, string) {
	msg := ""
	phenoMetaColl := session.DB(dbconf.Dbname).C(dbconf.PhenoMetaCollection)

	find := phenoMetaColl.Find(bson.M{"name": phenoname})

	items := find.Iter()
	phenoMeta := DBPhenoMeta{}

	if items.Next(&phenoMeta) != false {
		return false, "Phenotype already exists"
	}
	res := InsertPhenoMetaData(phenoname, phenosource, phenodesc, phenoclass)
	if res != true {
		msg = "Failed to insert pheno meta data"
		return res, msg
	}

	res = InsertPhenoData(phenoname, phenoitems)
	if res != true {
		msg = "Failed to insert pheno data"
	}
	return res, msg
}

// InsertPhenoMetaData ...
// Insert pheno meta data from string arguments
//---------------------------------------------------------------------
func InsertPhenoMetaData(name string, source string, desc string, pclass string) bool {
	phenoMetaColl := session.DB(dbconf.Dbname).C(dbconf.PhenoMetaCollection)
	dbdata := bson.M{"name": name, "source": source, "description": desc, "phenoclass": pclass}

	phenoMetaColl.Insert(dbdata)

	return true
}

// InsertPhenoData ...
// Insert pheno data from a list
//---------------------------------------------------------------------
func InsertPhenoData(phenoname string, phenoitems map[string]string) bool {

	phenoColl := session.DB(dbconf.Dbname).C(dbconf.PhenoCollection)

	for iid, value := range phenoitems {
		dbdata := bson.M{"name": phenoname, "iid": iid, "value": value}
		phenoColl.Insert(dbdata)
	}

	return true
}

// ******* Stored variantlist section ***********************************

// GetVarlistByName ...
// Get all entries for a variantlist by name
//---------------------------------------------------------------------
func GetVarlistByName(name string) ([]string, int) {
	varlist := make([]string, 0, 100)

	varlistColl := session.DB(dbconf.Dbname).C(dbconf.VarlistCollection)

	variant := DBListVar{}

	find := varlistColl.Find(bson.M{"name": name})

	items := find.Iter()
	count := 0
	for items.Next(&variant) {
		varlist = append(varlist, variant.VarID)
		count++
	}
	return varlist, count
}

// GetVarlistMetaByName ...
// Get single entry for the varlistMeta collection by name
//---------------------------------------------------------------------
func GetVarlistMetaByName(name string) DBVarlistMeta {
	varlistMetaColl := session.DB(dbconf.Dbname).C(dbconf.VarlistMetaCollection)

	varlistMeta := DBVarlistMeta{}

	find := varlistMetaColl.Find(bson.M{"name": name})

	items := find.Iter()
	for items.Next(&varlistMeta) {
	}
	return varlistMeta
}

// GetVarlistMetaNames ...
// Get all varlist names from the varlist meta collection
//---------------------------------------------------------------------
func GetVarlistMetaNames() []string {
	var varlistNameList = make([]string, 0, 10)

	varlistMetaColl := session.DB(dbconf.Dbname).C(dbconf.VarlistMetaCollection)

	varlistMeta := DBVarlistMeta{}

	find := varlistMetaColl.Find(bson.M{})

	items := find.Iter()
	for items.Next(&varlistMeta) {
		varlistNameList = append(varlistNameList, varlistMeta.Name)
	}
	return varlistNameList
}

// InsertVarlistDataWithCheck ...
// Insert variant list  data
// but first check for existence (in the varlistMetaColl)
//---------------------------------------------------------------------
func InsertVarlistDataWithCheck(name string, desc string,
	items []string) (bool, string) {
	msg := ""
	varlistMetaColl := session.DB(dbconf.Dbname).C(dbconf.VarlistMetaCollection)

	find := varlistMetaColl.Find(bson.M{"name": name})

	metaitems := find.Iter()
	varlistMeta := DBVarlistMeta{}

	if metaitems.Next(&varlistMeta) != false {
		return false, "Variant list already exists"
	}
	res := InsertVarlistMetaData(name, desc)
	if res != true {
		msg = "Failed to insert variant list meta data"
		return res, msg
	}

	res = InsertVarlistData(name, items)
	if res != true {
		msg = "Failed to insert variant list data"
	}
	return res, msg
}

// InsertVarlistMetaData ...
// Insert pheno meta data from string arguments
//---------------------------------------------------------------------
func InsertVarlistMetaData(name string, desc string) bool {
	varlistMetaColl := session.DB(dbconf.Dbname).C(dbconf.VarlistMetaCollection)
	dbdata := bson.M{"name": name, "description": desc}

	varlistMetaColl.Insert(dbdata)

	return true
}

// InsertVarlistData ...
// Insert variant list data from a list
//---------------------------------------------------------------------
func InsertVarlistData(name string, items []string) bool {

	varlistColl := session.DB(dbconf.Dbname).C(dbconf.VarlistCollection)

	for _, varid := range items {
		dbdata := bson.M{"name": name, "varid": varid}
		varlistColl.Insert(dbdata)
	}

	return true
}

// ******* Grs input section ***********************************

// GetGrsInputByName ...
// Get all entries for a phenotype by name
//---------------------------------------------------------------------
func GetGrsInputByName(name string) (map[string][]string, int) {
	grsInputData := make(map[string][]string)

	grsInputColl := session.DB(dbconf.Dbname).C(dbconf.GrsInputCollection)

	grsInput := DBGrsInput{}

	find := grsInputColl.Find(bson.M{"name": name})

	items := find.Iter()
	count := 0
	for items.Next(&grsInput) {
		var theRest = make([]string, 0, 3)
		theRest[0] = grsInput.Ea
		theRest[1] = grsInput.Eaf
		theRest[2] = grsInput.Weight
		grsInputData[grsInput.Varid] = theRest
		count++
	}
	return grsInputData, count
}

// GetGrsInputAsStringArrayByName ...
// Get all entries for a phenotype by name
//---------------------------------------------------------------------
func GetGrsInputAsStringArrayByName(name string) ([]string, int) {
	grsInputData := make([]string, 0)

	grsInputColl := session.DB(dbconf.Dbname).C(dbconf.GrsInputCollection)

	grsInput := DBGrsInput{}

	find := grsInputColl.Find(bson.M{"name": name})

	items := find.Iter()
	count := 0
	grsInputData = append(grsInputData, "varid,ea,eaf,wgt")
	for items.Next(&grsInput) {
		grsInputData = append(grsInputData, fmt.Sprintf("%s,%s,%s,%s", grsInput.Varid, grsInput.Ea, grsInput.Eaf, grsInput.Weight))
		count++
	}
	return grsInputData, count
}

// GetGrsMetaByName ...
// Get all entries for a phenotype by name
//---------------------------------------------------------------------
func GetGrsMetaByName(name string) DBGrsMeta {
	grsMetaColl := session.DB(dbconf.Dbname).C(dbconf.GrsMetaCollection)

	grsMeta := DBGrsMeta{}

	find := grsMetaColl.Find(bson.M{"name": name})

	items := find.Iter()
	for items.Next(&grsMeta) {
	}
	return grsMeta
}

// GetGrsMetaNames ...
// Get all pheno names from the pheno meta collection
//---------------------------------------------------------------------
func GetGrsMetaNames() []string {
	var grsNameList = make([]string, 0, 10)

	grsMetaColl := session.DB(dbconf.Dbname).C(dbconf.GrsMetaCollection)

	grsMeta := DBGrsMeta{}

	find := grsMetaColl.Find(bson.M{})

	items := find.Iter()
	for items.Next(&grsMeta) {
		grsNameList = append(grsNameList, grsMeta.Name)
	}
	return grsNameList
}

// InsertGrsMetaData ...
// Insert grs meta data from string arguments
//---------------------------------------------------------------------
func InsertGrsMetaData(name string, desc string) bool {
	grsMetaColl := session.DB(dbconf.Dbname).C(dbconf.GrsMetaCollection)
	dbdata := bson.M{"name": name, "description": desc}

	grsMetaColl.Insert(dbdata)

	return true
}

// InsertGrsInputDataWithCheck ...
// Insert GRS input data, usually a SNP list with effect-allele information and
// weight calculated in a previous study
//---------------------------------------------------------------------
func InsertGrsInputDataWithCheck(grsname string, grsdesc string, grsinputitems map[string][]string) (bool, string) {
	msg := ""
	grsMetaColl := session.DB(dbconf.Dbname).C(dbconf.GrsMetaCollection)

	find := grsMetaColl.Find(bson.M{"name": grsname})

	items := find.Iter()
	grsInput := DBGrsInput{}

	if items.Next(&grsInput) != false {
		return false, "GRS input name already exists"
	}

	res := InsertGrsMetaData(grsname, grsdesc)
	if res != true {
		msg = "Failed to insert grs meta data"
		return res, msg
	}

	res = InsertGrsInputData(grsname, grsinputitems)
	if res != true {
		msg = "Failed to insert grs input data"
	}
	return res, msg
}

// InsertGrsInputData ...
// Insert grs input data from a list of items
//---------------------------------------------------------------------
func InsertGrsInputData(name string, items map[string][]string) bool {

	grsInputColl := session.DB(dbconf.Dbname).C(dbconf.GrsInputCollection)

	for varid, value := range items {
		dbdata := bson.M{"name": name, "varid": varid, "ea": value[0], "eaf": value[1], "wgt": value[2]}
		grsInputColl.Insert(dbdata)
	}

	return true
}
