package main

import (
	"encoding/json"
	"log"
	"os"
	"strconv"
	"strings"
)

// Configuration ...
// For use throughout the app
type Configuration struct {
	Address        string `json:"addr"`
	ReadTimeout    int64  `json:"readto"`
	WriteTimeout   int64  `json:"writeto"`
	Static         string `json:"static"`
	Templates      string `json:"templates"`
	VcfPrfx        string `json:"vcfprfx"`
	Assaytypes     string `json:"assaytypes"`
	Pthr           string `json:"probthr"`
	AssocCmd       string `json:"assoccmd"`
	AssocBinaryCmd string `json:"assoccmdbin"`
	OutfilePath    string `json:"outfilepath"`
	PhenofilePath  string `json:"phenofilepath"`
	PhenoColumns   string `json:"phenocolumns"`
	GrsfilePath    string `json:"grsfilepath"`
	GrsColumns     string `json:"grscolumns"`
	Uploads        string `json:"uploads"`
}

var config Configuration
var logger *log.Logger
var requestedAssaytypes = map[string]bool{}
var validPhenoColumns = map[string]bool{}
var validGrsColumns = map[string]bool{}

func init() {
	loadConfig()
	logFile := os.Getenv("LOGFILE")
	file, err := os.OpenFile(logFile, os.O_CREATE|os.O_WRONLY|os.O_APPEND, 0666)
	if err != nil {
		log.Fatalln("Failed to open log file", err)
	}
	logger = log.New(file, "INFO ", log.Ldate|log.Ltime|log.Lshortfile)
}

func loadConfig() {
	configFile := os.Getenv("CONFIGFILE")
	file, err := os.Open(configFile)
	if err != nil {
		log.Fatalln("Cannot open config file", err)
	}
	decoder := json.NewDecoder(file)
	config = Configuration{}
	err = decoder.Decode(&config)
	if err != nil {
		log.Fatalln("Cannot get configuration from file", err)
	}
	atList := strings.Split(config.Assaytypes, ",")
	for at := range atList {
		requestedAssaytypes[atList[at]] = true
	}
	log.Printf("Req ATs: %v\n", requestedAssaytypes)
	colList := strings.Split(config.PhenoColumns, ",")
	for colname := range colList {
		validPhenoColumns[colList[colname]] = true
	}
	log.Printf("Valid pheno input cols: %v\n", validPhenoColumns)
	colList = strings.Split(config.GrsColumns, ",")
	for colname := range colList {
		validGrsColumns[colList[colname]] = true
	}
	log.Printf("Valid grs input cols: %v\n", validGrsColumns)
}

func getAssaytypes() map[string]bool {
	return requestedAssaytypes
}

func getPhenoColMap() map[string]bool {
	return validPhenoColumns
}

func getGrsColMap() map[string]bool {
	return validGrsColumns
}

func getThresholdAsFloat() float64 {
	threshold, err := strconv.ParseFloat(config.Pthr, 64)
	if err != nil {
		return 0.9 // conservative default
	}
	return threshold
}
