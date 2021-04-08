package main

import (
	"encoding/json"
	"log"
	"os"
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
	Uploads        string `json:"uploads"`
}

var config Configuration
var logger *log.Logger
var requestedAssaytypes = map[string]bool{}
var validPhenoColumns = map[string]bool{}

func init() {
	loadConfig()
	// logFile := os.Getenv("LOGFILE")
	logFile := "./logs/godbassoc.log"
	file, err := os.OpenFile(logFile, os.O_CREATE|os.O_WRONLY|os.O_APPEND, 0666)
	if err != nil {
		log.Fatalln("Failed to open log file", err)
	}
	logger = log.New(file, "INFO ", log.Ldate|log.Ltime|log.Lshortfile)
}

func loadConfig() {
	// configFile := os.Getenv("CONFIGFILE")
	configFile := "./cfg/config.json"
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
	log.Printf("Valid input cols: %v\n", validPhenoColumns)
}
func getAssaytypes() map[string]bool {
	return requestedAssaytypes
}
func getPhenoColMap() map[string]bool {
	return validPhenoColumns
}
