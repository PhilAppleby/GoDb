package main

import (
	"bufio"
	"bytes"
	"ehrdb"
	"fmt"
	"godb"
	"grs"
	"html/template"
	"io/ioutil"
	"log"
	"net/http"
	"strings"
	"time"
)

func grsUpload(w http.ResponseWriter, r *http.Request) {
	var data []string
	t := template.Must(template.ParseFiles(
		config.Templates+"/grsupload.html",
		config.Templates+"/navigation.html"))
	t.ExecuteTemplate(w, "grsupload", data)
}

func grsFileProcess(w http.ResponseWriter, r *http.Request) {
	start := time.Now()
	log.Printf("File Upload Endpoint Hit\n")
	log.Printf("%v\n", r)

	// Parse the multipart form, 10 << 20 specifies a maximum
	// upload of 10 MB files.
	r.ParseMultipartForm(10 << 20)
	// FormFile returns the first file for the given key `phenoFile`
	// it also returns the FileHeader so we can get the Filename,
	// the Header and the size of the file
	file, handler, err := r.FormFile("grsFile")
	if err != nil {
		log.Printf("File retrieve err %v\n", err)
		return
	}
	defer file.Close()
	log.Printf("Uploaded File: %+v\n", handler.Filename)
	log.Printf("File Size: %+v\n", handler.Size)
	log.Printf("MIME Header: %+v\n", handler.Header)

	log.Printf("Form value (gname): %+v", r.FormValue("gname"))
	gname := r.FormValue("gname")
	gdesc := r.FormValue("gdesc")
	// read all of the contents of the uploaded file into a
	// byte array
	fileBytes, err := ioutil.ReadAll(file)
	if err != nil {
		log.Printf("File readall err %v\n", err)
	}
	// write this byte array to our temporary file
	err = ioutil.WriteFile(config.GrsfilePath+"/"+handler.Filename, fileBytes, 0644)
	if err != nil {
		log.Printf("Write err %v\n", err)
		errorMessage(w, r, fmt.Sprintf("Write err %v\n", err))
	}
	scanner := bufio.NewScanner(bytes.NewReader(fileBytes))
	scanner.Scan()
	hdr := scanner.Text()
	hdrData := strings.Split(hdr, ",")
	log.Printf("HDR=%v", hdrData)
	if grsHdrIsValid(hdrData) == false {
		errorMessage(w, r, handler.Filename+": Invalid file header detected expected: "+config.GrsColumns+" got: "+hdr)
	} else {
		var grsLines []string
		grsLines = append(grsLines, hdr)
		grsitems := make(map[string][]string)
		for scanner.Scan() {
			line := scanner.Text()
			grsLines = append(grsLines, line)
			lineData := strings.Split(line, ",")
			grsitems[lineData[0]] = lineData[1:]
		}
		res, msg := ehrdb.InsertGrsInputDataWithCheck(gname, gdesc, grsitems)
		if res != true {
			errorMessage(w, r, "GRS file upload failed for "+gname+" "+msg)
			elapsed := time.Since(start)
			log.Printf("res: Failed time = %s", elapsed)
		} else {
			rsidList, eaMap, eafMap, wgtMap := grs.GetGrsMaps(grsLines)

			_, _, genorecs := godb.Getallvardata(config.VcfPrfx, rsidList, getAssaytypes(), getThresholdAsFloat())

			_, grScoresFlip := grs.GetScores(genorecs, eaMap, eafMap, wgtMap)

			for _, gScore := range grScoresFlip {
				log.Printf("%s\n", gScore)
			}
			elapsed := time.Since(start)
			log.Printf("res: GRS save / calc time = %s", elapsed)
			url := []string{"/index"}
			http.Redirect(w, r, strings.Join(url, ""), 302)
		}
	}
}

//  hdrIsValid ...
//  Does the supplied header record contain all required fields?
func grsHdrIsValid(hdr []string) bool {
	colmap := getGrsColMap()
	numcols := len(colmap)
	log.Printf("numcols: %d\n", numcols)
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
