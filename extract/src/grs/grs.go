// Package grs ...
//
package grs

import (
	"fmt"
	"log"
	"strings"
	"variant"
)

var firstSampleCol int

// GetScores ...
//---------------------------------------------------------------------
func GetScores(records []string, eas map[string]string, eafs map[string]float64, wgts map[string]float64) []string {
	firstSampleCol = 9

	genoScores := make(map[string]float64)
	varCounts := make(map[string]int)

	sampleData := strings.Split(records[0], "\t")[firstSampleCol:]
	scoreList := make([]string, 0, len(sampleData))

	for _, record := range records[1:] {
		recData := strings.Split(record, "\t")
		genoData := recData[firstSampleCol:]
		varid := variant.GetVarid(recData)
		_, altAllele := variant.GetAlleles(recData)
		if eas[varid] != altAllele {
			log.Printf("FLIP for: %s [%s,%s]\n", varid, altAllele, eas[varid])
		}
		for i, geno := range genoData {
			if _, ok := genoScores[sampleData[i]]; !ok {
				genoScores[sampleData[i]] = 0.0
				varCounts[sampleData[i]] = 0
			}
			if !strings.HasPrefix(geno, ".") {
				genoScores[sampleData[i]] += float64(getGenoIntValue(geno, altAllele, eas[varid])) * wgts[varid]
				varCounts[sampleData[i]]++
			}
		}
	}
	for k, v := range genoScores {
		scoreLine := fmt.Sprintf("%s,%d,%.3f", k, varCounts[k], v)
		scoreList = append(scoreList, scoreLine)
	}
	return scoreList
}
func getGenoIntValue(geno string, alt string, ea string) int {
	gval := 0

	genoData := strings.Split(geno, ":")

	if genoData[0] == "1/0" || genoData[0] == "0/1" {
		gval = 1
	}
	if genoData[0] == "0/0" && alt == ea {
		gval = 0
	}
	if genoData[0] == "0/0" && alt != ea {
		gval = 2
	}
	if genoData[0] == "1/1" && alt == ea {
		gval = 2
	}
	if genoData[0] == "1/1" && alt != ea {
		gval = 0
	}
	return gval
}
