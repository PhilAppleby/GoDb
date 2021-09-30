package main

import (
	"bytes"
	"compress/gzip"
	"godb"
	"io/ioutil"
	"log"
	"net/http"
	"strconv"
	"strings"
	"time"
	"variant"
)

func dataDownload(w http.ResponseWriter, r *http.Request) {
	log.Printf("download Query len: %d\n", len(r.URL.Query()))
	log.Printf("download Query: %v\n", r.URL.Query())
	fmtChoice := "vcf"

	if len(r.URL.Query()) != 0 {
		log.Printf("download resbtn: %v\n", r.URL.Query()["resbtn"][0])
		if r.URL.Query()["resbtn"][0] == "Download" {
			if _, ok := r.URL.Query()["variant"]; !ok {
				err(w, r) // Internal error handling
			}
			if _, ok := r.URL.Query()["pthr"]; !ok {
				err(w, r) // Internal error handling
			}
			if _, ok := r.URL.Query()["ffmt"]; ok {
				log.Printf("file fmt: %s", r.URL.Query()["ffmt"][0])
				fmtChoice = r.URL.Query()["ffmt"][0]
			}
			start := time.Now()
			var variantList = make([]string, 0)
			variantList = append(variantList, r.URL.Query()["variant"][0])
			pthr, _ := strconv.ParseFloat(r.URL.Query()["pthr"][0], 64)
			_, _, comborecs := godb.Getallvardata(config.VcfPrfx, variantList, getAssaytypes(), pthr)
			// content, outFmt := godb.FormatOutput(comborecs, fmtChoice)
			outFileName := config.OutfilePath + "/" + variantList[0] + "." + fmtChoice + ".gz"
			dnldFileName := variantList[0] + "." + fmtChoice + ".gz"

			// Initialize gzip
			buf := &bytes.Buffer{}
			gzWriter := gzip.NewWriter(buf)
			outrecs := makeOutputRecords(comborecs, fmtChoice)
			gzWriter.Write([]byte(strings.Join(outrecs, "\n") + "\n"))
			gzWriter.Close()
			ioutil.WriteFile(outFileName, buf.Bytes(), 0644)
			w.Header().Set("Content-Disposition", "attachment; filename="+strconv.Quote(dnldFileName))
			w.Header().Set("Content-Type", "application/octet-stream")
			elapsed := time.Since(start)
			log.Printf("build file for download took %s", elapsed)
			http.ServeFile(w, r, outFileName)
		} else {
			url := []string{"/index"}
			http.Redirect(w, r, strings.Join(url, ""), 302)
		}
	} else {
		url := []string{"/index"}
		http.Redirect(w, r, strings.Join(url, ""), 302)
	}
}
func makeOutputRecords(recs []string, fmtChoice string) []string {
	if fmtChoice != "vcf" {
		recs = makeCsvData(recs)
	}
	return recs
}
func makeCsvData(recs []string) []string {
	rtnrecs := make([]string, len(recs))
	if len(recs) > 0 {
		_, samples := variant.GetVCFPrfxSfx(strings.Split(recs[0], "\t"))

	}
	return rtnrecs
}
