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
			outrecs := makeOutputRecords(comborecs, fmtChoice, pthr)
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
func makeOutputRecords(recs []string, fmtChoice string, pthr float64) []string {
	if fmtChoice != "vcf" {
		recs = makeCsvData(recs, pthr)
	}
	return recs
}

// This is a transposition function ...
func makeCsvData(recs []string, pthr float64) []string {
	samplerecs := make(map[string][]string)
	nulrecs := make([]string, 0)
	if len(recs) > 0 {
		_, samples := variant.GetVCFPrfxSfx(strings.Split(recs[0], "\t"))
		rtnrecs := make([]string, len(samples)+1)
		hdrdata := make([]string, len(recs)+1)
		hdrdata[0] = "FID"
		hdrdata[1] = "IID"
		for _, sample := range samples {
			samplerecs[sample] = make([]string, len(recs)+1)
		}
		for j, rec := range recs[1:] {
			fields := strings.Split(rec, "\t")
			prfx, genotypes := variant.GetVCFPrfxSfx(fields)
			hdrdata[j+2] = variant.GetVarid(prfx)
			pindex := variant.GetProbIdx(prfx)
			for i, sample := range samples {
				samplerecs[sample][0] = sample
				samplerecs[sample][1] = sample
				samplerecs[sample][j+2] = variant.GetGenoAsIntStr(genotypes[i], pthr, pindex)
			}
		}
		rtnrecs[0] = strings.Join(hdrdata, ",")
		j := 1
		for _, v := range samplerecs {
			rtnrecs[j] = strings.Join(v, ",")
			j++
		}
		return rtnrecs
	}
	// no data return empty string array
	return nulrecs
}
