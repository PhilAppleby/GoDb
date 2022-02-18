//
// Search / retrieval for a godb instance
// Also includes single SNP assocation testing And
// phenotype upload
//
package main

import (
	"fmt"
	"log"
	"net/http"
	"time"
)

func init() {
}

func check(e error, source string) {
	if e != nil {
		fmt.Println("main: check routine source: ", source)
		fmt.Println("err != nil")
		log.Fatal(e)
	}
}

func main() {
	p("godbassoc", version(), "started at", config.Address)

	// handle static assets
	mux := http.NewServeMux()
	files := http.FileServer(http.Dir(config.Static))
	mux.Handle("/public/", http.StripPrefix("/public/", files))
	// route pattern matching
	// index, defined in route_main.go
	mux.HandleFunc("/", index)
	mux.HandleFunc("/index", index)
	mux.HandleFunc("/grsrun", grsrun)
	mux.HandleFunc("/results/", results)
	mux.HandleFunc("/grsresults/", grsresults)

	mux.HandleFunc("/notes", notes)
	mux.HandleFunc("/err", err)

	// defined in route_pheno.go
	mux.HandleFunc("/phenoupload", phenoUpload)
	mux.HandleFunc("/phenoprocess", phenoFileProcess)

	// defined in route_varlist.go
	mux.HandleFunc("/varlistupload", varlistUpload)
	mux.HandleFunc("/varlistprocess", varlistFileProcess)

	// defined in route_grs.go
	mux.HandleFunc("/grsupload", grsUpload)
	mux.HandleFunc("/grsprocess", grsFileProcess)

	// defined in route_genodata.go
	mux.HandleFunc("/download", dataDownload)

	// starting up the server
	server := &http.Server{
		Addr:           config.Address,
		Handler:        mux,
		ReadTimeout:    time.Duration(config.ReadTimeout * int64(time.Second)),
		WriteTimeout:   time.Duration(config.WriteTimeout * int64(time.Second)),
		MaxHeaderBytes: 1 << 20,
	}
	server.ListenAndServe()
}
