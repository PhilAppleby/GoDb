// Package grs ...
//
package grs

import (
	"fmt"
	"log"
	"strconv"
	"strings"
	"variant"
)

var firstSampleCol int

// GetScores ...
//---------------------------------------------------------------------
func GetScores(records []string, eas map[string]string, eafs map[string]float64, wgts map[string]float64) ([]string, []string) {
	firstSampleCol = 9

	genoScores := make(map[string]float64)
	genoScoresFlip := make(map[string]float64)
	varCounts := make(map[string]int)

	sampleData := strings.Split(records[0], "\t")[firstSampleCol:]
	scoreList := make([]string, 0, len(sampleData))
	scoreListFlip := make([]string, 0, len(sampleData))

	for _, record := range records[1:] {
		mult := 1.0
		recData := strings.Split(record, "\t")
		genoData := recData[firstSampleCol:]
		varid := variant.GetVarid(recData)
		refAllele, altAllele := variant.GetAlleles(recData)
		if (eas[varid] != refAllele) && (eas[varid] != altAllele) {
			log.Printf("REJect: %s [%s,%s,%s]\n", varid, refAllele, altAllele, eas[varid])
			continue
		}
		if eas[varid] == refAllele {
			log.Printf("Sign Change / FLIP for: %s [%s,%s,%s]\n", varid, refAllele, altAllele, eas[varid])
			mult = -1.0
		}
		for i, geno := range genoData {
			if _, ok := genoScores[sampleData[i]]; !ok {
				genoScores[sampleData[i]] = 0.0
				genoScoresFlip[sampleData[i]] = 0.0
				varCounts[sampleData[i]] = 0
			}
			if !strings.HasPrefix(geno, ".") {
				genoScores[sampleData[i]] += float64(getGenoIntValue(geno)) * wgts[varid] * mult
				genoScoresFlip[sampleData[i]] += float64(getGenoIntValueFlip(geno, altAllele, eas[varid])) * wgts[varid]
				varCounts[sampleData[i]]++
			}
		}
	}
	log.Printf("genoscore lengths %d, %d\n", len(genoScores), len(genoScoresFlip))
	for k, v := range genoScores {
		scoreLine := fmt.Sprintf("%s,%d,%.6f", k, varCounts[k], v)
		scoreList = append(scoreList, scoreLine)
	}

	for k, v := range genoScoresFlip {
		scoreLine := fmt.Sprintf("%s,%d,%.6f", k, varCounts[k], v)
		scoreListFlip = append(scoreListFlip, scoreLine)
	}
	return scoreList, scoreListFlip
}
func getGenoIntValueFlip(geno string, alt string, ea string) int {
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
func getGenoIntValue(geno string) int {
	gval := 0

	genoData := strings.Split(geno, ":")

	if genoData[0] == "1/0" || genoData[0] == "0/1" {
		gval = 1
	}
	if genoData[0] == "0/0" {
		gval = 0
	}
	if genoData[0] == "1/1" {
		gval = 2
	}
	return gval
}

// GetGrsMaps ...
// Break up GRS input into consituent parts
func GetGrsMaps(rsIDLines []string) ([]string, map[string]string, map[string]float64, map[string]float64) {
	rsIDList := make([]string, 0, len(rsIDLines))
	eaMap := make(map[string]string)
	eafMap := make(map[string]float64)
	wgtMap := make(map[string]float64)
	colMap := make(map[string]int)

	hdrStr := rsIDLines[0]
	hdrData := strings.Split(hdrStr, ",")
	for i, col := range hdrData {
		colMap[col] = i
	}

	for _, rLine := range rsIDLines[1:] {
		grsData := strings.Split(rLine, ",")
		rsIDList = append(rsIDList, grsData[colMap["varid"]])
		eaMap[grsData[colMap["varid"]]] = strings.ToUpper(grsData[colMap["ea"]])
		eafMap[grsData[colMap["varid"]]], _ = strconv.ParseFloat(grsData[colMap["eaf"]], 64)
		wgtMap[grsData[colMap["varid"]]], _ = strconv.ParseFloat(grsData[colMap["wgt"]], 64)
		log.Printf("%s : %f\n", grsData[colMap["varid"]], wgtMap[grsData[colMap["varid"]]])
	}
	return rsIDList, eaMap, eafMap, wgtMap
}
