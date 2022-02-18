package main

import (
	"bufio"
	"bytes"
	"ehrdb"
	"fmt"
	"html/template"
	"io/ioutil"
	"log"
	"net/http"
	"strings"
)

func varlistUpload(w http.ResponseWriter, r *http.Request) {
	var data []string
	t := template.Must(template.ParseFiles(
		config.Templates+"/varlistupload.html",
		config.Templates+"/navigation.html"))
	t.ExecuteTemplate(w, "varlistupload", data)
}

func varlistFileProcess(w http.ResponseWriter, r *http.Request) {
	log.Printf("Varlist File Upload Endpoint Hit\n")
	log.Printf("%v\n", r)

	// Parse the multipart form, 10 << 20 specifies a maximum
	// upload of 10 MB files.
	r.ParseMultipartForm(10 << 20)
	// FormFile returns the first file for the given key `phenoFile`
	// it also returns the FileHeader so we can get the Filename,
	// the Header and the size of the file
	file, handler, err := r.FormFile("varlistFile")
	if err != nil {
		log.Printf("File retrieve err %v\n", err)
		return
	}
	defer file.Close()
	log.Printf("Uploaded File: %+v\n", handler.Filename)
	log.Printf("File Size: %+v\n", handler.Size)
	log.Printf("MIME Header: %+v\n", handler.Header)

	log.Printf("Form value (vlname): %+v", r.FormValue("vlname"))
	vlname := r.FormValue("vlname")
	vldesc := r.FormValue("vldesc")
	// read all of the contents of the uploaded file into a
	// byte array
	fileBytes, err := ioutil.ReadAll(file)
	if err != nil {
		log.Printf("File readall err %v\n", err)
	}
	// write this byte array to our temporary file
	err = ioutil.WriteFile(config.VarlistfilePath+"/"+handler.Filename, fileBytes, 0644)
	if err != nil {
		log.Printf("Write err %v\n", err)
		errorMessage(w, r, fmt.Sprintf("Write err %v\n", err))
	}
	scanner := bufio.NewScanner(bytes.NewReader(fileBytes))
	scanner.Scan()
	hdr := scanner.Text()
	hdrData := strings.Split(hdr, ",")
	log.Printf("Varlist HDR=%v", hdrData)
	if vlhdrIsValid(hdrData) == false {
		errorMessage(w, r, handler.Filename+": Invalid file header detected expected: "+config.VarlistColumns+" got: "+hdr)
	} else {
		varlist := make([]string, 0)
		for scanner.Scan() {
			lineData := strings.Split(scanner.Text(), ",")
			varlist = append(varlist, lineData[0])
		}
		res, msg := ehrdb.InsertVarlistDataWithCheck(vlname, vldesc, varlist)
		if res != true {
			errorMessage(w, r, "Varlist upload failed for "+vlname+" "+msg)
		} else {
			url := []string{"/index"}
			http.Redirect(w, r, strings.Join(url, ""), 302)
		}
	}
}

//  hdrIsValid ...
//  Does the supplied header record contain all required fields?
//  TODO - this needs to be refactored and combined with similar methods
func vlhdrIsValid(hdr []string) bool {
	colmap := getVarlistColMap()
	numcols := len(colmap)
	log.Printf("varlist numcols: %d\n", numcols)
	matchcols := 0

	for colnum := range hdr {
		//log.Printf("colname: %s\n", hdr[colnum])
		if _, ok := colmap[hdr[colnum]]; ok {
			matchcols++
		}
	}
	log.Printf("matchcols: %d\n", matchcols)
	if numcols == matchcols {
		return true
	}
	return false
}
