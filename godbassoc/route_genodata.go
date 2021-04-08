package main

import (
	"bufio"
	"godb"
	"log"
	"net/http"
	"os"
	"strconv"
	"strings"
)

func dataDownload(w http.ResponseWriter, r *http.Request) {
	log.Printf("download Query len: %d\n", len(r.URL.Query()))
	log.Printf("download Query: %v\n", r.URL.Query())

	if len(r.URL.Query()) != 0 {
		log.Printf("download resbtn: %v\n", r.URL.Query()["resbtn"][0])
		if r.URL.Query()["resbtn"][0] == "Download" {
			if _, ok := r.URL.Query()["variant"]; !ok {
				err(w, r) // Internal error handling
			}
			if _, ok := r.URL.Query()["pthr"]; !ok {
				err(w, r) // Internal error handling
			}
			var variantList = make([]string, 0)
			variantList = append(variantList, r.URL.Query()["variant"][0])
			pthr, _ := strconv.ParseFloat(r.URL.Query()["pthr"][0], 64)
			_, _, comborecs := godb.Getallvardata(config.VcfPrfx, variantList, getAssaytypes(), pthr)
			outFileName := config.OutfilePath + "/" + variantList[0] + ".vcf"
			f2, _ := os.Create(outFileName)

			defer f2.Close()
			wf := bufio.NewWriter(f2)
			log.Printf("Number vcf records: %v\n", len(comborecs))
			for _, line := range comborecs {
				wf.WriteString(line + "\n")
			}
			wf.Flush()
		}
	}
	url := []string{"/index"}
	http.Redirect(w, r, strings.Join(url, ""), 302)
}
