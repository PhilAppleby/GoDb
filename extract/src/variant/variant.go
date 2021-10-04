// Package variant ...
// Everything VCF variant related
// Author: P Appleby, UoD
package variant

import (
	"fmt"
	"log"
	"strconv"
	"strings"
)

var sglDigitChrom map[string]int
var genoDelim = "/"
var genoStrings = []string{"0/0", "0/1", "1/1"}
var genoStringInts = []string{"0", "1", "2"}

const firstGenoIdx = 9
const chrIdx = 0
const posnIdx = 1
const varIdx = 2
const refIdx = 3
const altIdx = 4
const qcIdx = 5
const filtIdx = 6
const infoIdx = 7
const fmtIdx = 8

func init() {
	sglDigitChrom = make(map[string]int)
	sglDigitChrom["affy"] = 1
	sglDigitChrom["illumina"] = 1
	sglDigitChrom["broad"] = 1
	sglDigitChrom["exome"] = 1
	sglDigitChrom["metabo"] = 1
	sglDigitChrom["bigtest"] = 1
}

func check(e error) {
	//  fmt.Println("check")
	if e != nil {
		fmt.Println("err != nil")
		log.Fatal(e)
	}
}

// AppendToFmt ...
func AppendToFmt(prfx []string, addStr string) []string {
	prfx[fmtIdx] = prfx[fmtIdx] + ":" + addStr
	return prfx
}

// NormaliseChromosome ...
func NormaliseChromosome(prfx []string) []string {
	if prfx[chrIdx][0] == '0' {
		prfx[chrIdx] = prfx[chrIdx][1:]
	}
	return prfx
}

// SetInfoValue ...
func SetInfoValue(prfx []string, infoscore float64) []string {
	infomap := parseInfoStr(GetInfo(prfx))
	infomap["INFO"] = fmt.Sprintf("%.5f", infoscore)
	prfx[infoIdx] = buildinfoStr(infomap)
	return prfx
}

// GetGeno ...
// based on imputation probability threshold
//------------------------------------------------------------------------------
func GetGeno(geno string, threshold float64, probidx int) string {

	if geno == "." {
		return geno
	}

	mprob, maxProbIdx, genoarray := MaxProb(geno, probidx)
	//fmt.Printf("GGENO %s,%f,%f,%d\n", geno, mprob, threshold, probidx)
	if mprob < threshold {
		genoarray[0] = "./."
	} else {
		genoarray[0] = genoStrings[maxProbIdx]
	}
	return strings.Join(genoarray, ":")
}

// GetGenoAsIntStr ...
// based on imputation probability threshold
//------------------------------------------------------------------------------
func GetGenoAsIntStr(geno string, threshold float64, probidx int) string {

	if geno == "." {
		return geno
	}

	mprob, maxProbIdx, _ := MaxProb(geno, probidx)
	if mprob < threshold {
		return "-9"
	}
	return genoStringInts[maxProbIdx]
}

// MaxProb ...
// test for which slot contains the max probability for a genotype
//------------------------------------------------------------------------------
func MaxProb(geno string, probidx int) (float64, int, []string) {
	g := strings.Split(geno, ":")

	maxProb := 0.0
	maxProbIdx := -9
	probs := strings.Split(g[probidx], ",")

	for i, prob := range probs {
		probf, _ := strconv.ParseFloat(prob, 64)
		if probf > maxProb {
			maxProb = probf
			maxProbIdx = i
		}
	}
	return maxProb, maxProbIdx, g
}

// GetVCFPrfxSfx ...
func GetVCFPrfxSfx(recslice []string) ([]string, []string) {
	return recslice[:firstGenoIdx], recslice[firstGenoIdx:]
}

// GetChrom ...
func GetChrom(recslice []string) string {
	return recslice[chrIdx]
}

// GetVarid ...
func GetVarid(recslice []string) string {
	return recslice[varIdx]
}

// GetPosn ...
func GetPosn(recslice []string) int {
	value, _ := strconv.ParseInt(recslice[posnIdx], 0, 32)
	return int(value)
}

// GetPosnStr ...
func GetPosnStr(recslice []string) string {
	return recslice[posnIdx]
}

// GetAlleles ...
func GetAlleles(recslice []string) (string, string) {
	return recslice[refIdx], recslice[altIdx]
}

// GetA ...
func GetA(recslice []string) string {
	return recslice[refIdx]
}

// GetB ...
func GetB(recslice []string) string {
	return recslice[altIdx]
}

// GetProbIdx ...
func GetProbIdx(recslice []string) int {
	return getStrIdx(recslice[fmtIdx], "GP")
}

// HasFmt ...
func HasFmt(recslice []string, fmt string) bool {
	if getStrIdx(recslice[fmtIdx], fmt) == -9 {
		return false
	}
	return true
}

// GetRefPanelAF ...
func GetRefPanelAF(recslice []string) float64 {
	infomap := parseInfoStr(GetInfo(recslice))
	if maf, ok := infomap["RefPanelAF"]; ok {
		refpaf, _ := strconv.ParseFloat(maf, 64)
		return refpaf
	}
	return 0.0
}

// GetInfoScore ...
func GetInfoScore(recslice []string) float64 {
	infomap := parseInfoStr(GetInfo(recslice))
	if info, ok := infomap["INFO"]; ok {
		infoscore, _ := strconv.ParseFloat(info, 64)
		return infoscore
	}
	return 1.0
}

// GetInfo ...
func GetInfo(recslice []string) string {
	return recslice[infoIdx]
}

func parseInfoStr(infoStr string) map[string]string {
	infomap := make(map[string]string)
	infodata := strings.Split(infoStr, ";")
	for _, elem := range infodata {
		kv := strings.Split(elem, "=")
		if len(kv) == 2 {
			infomap[kv[0]] = kv[1]
		}
	}
	return infomap
}
func buildinfoStr(infomap map[string]string) string {
	infostr := ""
	return infostr
}

func getStrIdx(str string, matchStr string) int {
	strArr := strings.Split(str, ":")
	for i, mstr := range strArr {
		if mstr == matchStr {
			return i
		}
	}
	return -9
}
