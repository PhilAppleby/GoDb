package main

import (
	"ehrdb"
	"godb"
	"html/template"
	"log"
	"net/http"
	"os/exec"
	"strconv"
	"strings"
	"time"
)

// VariantData ...
type VariantData struct {
	Variant       string
	DataList      []godb.DBVariant
	ComboDataList []godb.DBVariant
	Pthr          float64
	PhenoName     string
	PhenoSource   string
	PhenoDesc     string
	PhenoClass    string
	PhenoCount    int
	AssocResults  []string
}

// PhenoData ...
type PhenoData struct {
	PhenoList []string
	Pthr      float64
}

// GET /err?msg=
// shows the error message page
func err(w http.ResponseWriter, r *http.Request) {
	t := template.Must(template.ParseFiles(
		config.Templates+"/error.html",
		config.Templates+"/navigation.html"))
	data := r.URL.Query()["msg"][0]
	t.ExecuteTemplate(w, "error", data)
}

// Registered handler for bare "index" and below
func index(w http.ResponseWriter, r *http.Request) {
	var data PhenoData
	data.PhenoList = ehrdb.GetPhenoMetaNames()
	t := template.Must(template.ParseFiles(
		config.Templates+"/index.html",
		config.Templates+"/navigation.html"))
	data.Pthr, _ = strconv.ParseFloat(config.Pthr, 64)
	t.ExecuteTemplate(w, "index", data)
}

// Registered handler for "notes"
func notes(w http.ResponseWriter, r *http.Request) {
	var data string
	t := template.Must(template.ParseFiles(
		config.Templates+"/notes.html",
		config.Templates+"/navigation.html"))
	t.ExecuteTemplate(w, "notes", data)
}

// Registered handler for bare "result" and below.
func results(w http.ResponseWriter, r *http.Request) {

	if len(r.URL.Query()) != 0 {
		var data VariantData
		rsidList := strings.Split(r.URL.Query()["variant"][0], ",")
		data.Variant = r.URL.Query()["variant"][0]
		data.Pthr, _ = strconv.ParseFloat(config.Pthr, 64)
		tmppthr, err := strconv.ParseFloat(r.URL.Query()["pthr"][0], 64)
		if err == nil {
			data.Pthr = tmppthr
		}
		start := time.Now()
		variants, combinedvariants, genorecs := godb.Getallvardata(config.VcfPrfx, rsidList, getAssaytypes(), data.Pthr)
		elapsed := time.Since(start)
		log.Printf("res: dbaccess took %s", elapsed)
		data.DataList = variants
		data.ComboDataList = combinedvariants
		phenoName := r.URL.Query()["pheno"][0]

		if phenoName == "None" {
			t := template.Must(template.ParseFiles(
				config.Templates+"/results.html",
				config.Templates+"/vartables.html",
				config.Templates+"/navigation.html"))
			t.ExecuteTemplate(w, "results", data)
		} else {
			data.PhenoName = phenoName
			phenoMeta := ehrdb.GetPhenoMetaByName(phenoName)
			data.PhenoSource = phenoMeta.Source
			data.PhenoDesc = phenoMeta.Description
			data.PhenoClass = phenoMeta.PhenoClass
			// Get phenotype data from ehrdb
			phenoData, pCount := ehrdb.GetPhenoByName(phenoName)
			data.PhenoCount = pCount
			phenoFileName := config.PhenofilePath + "/" + phenoName + ".csv"
			writeBufferedFile(phenoFileName, convertStringMapToCSV(phenoData))
			genoFileName := config.OutfilePath + "/temp.vcf"
			writeBufferedFile(genoFileName, genorecs)
			cmd := exec.Command(config.AssocBinaryCmd, genoFileName, phenoFileName)
			if data.PhenoClass != "Binary" {
				cmd = exec.Command(config.AssocCmd, genoFileName, phenoFileName)
			}
			//err = cmd.Run()
			res, err := cmd.Output()
			check(err)
			data.AssocResults = assocReformat(string(res))
			//fmt.Printf("%s\n", "combined"+"\t"+colhdr_str)
			t := template.Must(template.ParseFiles(
				config.Templates+"/assocresults.html",
				config.Templates+"/vartables.html",
				config.Templates+"/navigation.html"))
			elapsed := time.Since(start)
			log.Printf("res: dbaccess + assoctest took %s", elapsed)
			t.ExecuteTemplate(w, "assocresults", data)
		}
	} else {
		url := []string{"/index"}
		http.Redirect(w, r, strings.Join(url, ""), 302)
	}
}
