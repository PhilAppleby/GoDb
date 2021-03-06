//---------------------------------------------------------
// File: variant.go
// Everything VCF variant related
// Author: P Appleby, University of Dundee
//---------------------------------------------------------
package variant

import (
	"fmt"
	"log"
	"strconv"
	"strings"
)

var sglDigitChrom map[string]int
var geno_strings = []string{"0/0", "0/1", "1/1"}

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
	sglDigitChrom["atest"] = 1
	sglDigitChrom["atest2"] = 1
	sglDigitChrom["atest3"] = 1
	sglDigitChrom["atest4"] = 1
	sglDigitChrom["atest5"] = 1
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

func AppendToFmt(prfx []string, add_str string) []string {
	prfx[fmtIdx] = prfx[fmtIdx] + ":" + add_str
	return prfx
}

func NormaliseChromosome(prfx []string) []string {
	if prfx[chrIdx][0] == '0' {
		prfx[chrIdx] = prfx[chrIdx][1:]
	}
	return prfx
}
// set INFO
func SetInfoValue(prfx []string, infoscore float64) []string {
	infomap := parseInfoStr(GetInfo(prfx))
  infomap["INFO"] = fmt.Sprintf("%.5f", infoscore)
  prfx[infoIdx] = buildinfoStr(infomap)
	return prfx
}

//------------------------------------------------------------------------------
// get geno based on threshold
//------------------------------------------------------------------------------
func Get_geno(geno string, threshold float64, probidx int) string {

  if geno == "." || geno == "./." {
    return geno
  }

	mprob, max_prob_idx, genoarray := MaxProb(geno, probidx)
	//fmt.Printf("GGENO %s,%f,%f,%d\n", geno, mprob, threshold, probidx)
	if mprob < threshold {
		genoarray[0] = "./."
	} else {
		genoarray[0] = geno_strings[max_prob_idx]
	}
	return strings.Join(genoarray, ":")
}

//------------------------------------------------------------------------------
// maxprob test for a genotype
//------------------------------------------------------------------------------
func MaxProb(geno string, probidx int) (float64, int, []string) {
	g := strings.Split(geno, ":")
	//if len(g) < (probidx + 1) {
	// fmt.Printf("GENO IDX PROBLEM %d, %s\n", probidx, geno)
	//	return 0.0, -9, g
	//}
	max_prob := 0.0
	max_prob_idx := -9
	probs := strings.Split(g[probidx], ",")

	for i, prob := range probs {
	 probf, _ := strconv.ParseFloat(prob, 64)
	 if probf > max_prob {
	  max_prob = probf
	  max_prob_idx = i
	 }
	}
	return max_prob, max_prob_idx, g
}

//------------------------------------------------------------------------------
// ---- a group of convenience fns to get array fields
//------------------------------------------------------------------------------
func GetVCFPrfx_Sfx(recslice []string) ([]string, []string) {
	return recslice[:firstGenoIdx], recslice[firstGenoIdx:]
}

func GetChrom(recslice []string) string {
	return recslice[chrIdx]
}

func GetVarid(recslice []string) string {
	return recslice[varIdx]
}

func GetPosn(recslice []string) int {
	value, _ := strconv.ParseInt(recslice[posnIdx], 0, 32)
	return int(value)
}

func GetPosnStr(recslice []string) string {
	return recslice[posnIdx]
}

func GetAlleles(recslice []string) (string, string) {
	return recslice[refIdx], recslice[altIdx]
}

func GetA(recslice []string) string {
	return recslice[refIdx]
}

func GetB(recslice []string) string {
	return recslice[altIdx]
}

func GetProbidx(recslice []string) int {
	return getStrIdx(recslice[fmtIdx], "GP")
}

func HasFmt(recslice []string, fmt string) bool {
  if getStrIdx(recslice[fmtIdx], fmt) == -9 {
    return false
  }
  return true
}

func GetRefPanelAF(recslice []string) float64 {
	infomap := parseInfoStr(GetInfo(recslice))
	if maf, ok := infomap["RefPanelAF"]; ok {
		refpaf, _ := strconv.ParseFloat(maf, 64)
		return refpaf
	}
	return 0.0
}

func GetInfoScore(recslice []string) float64 {
	infomap := parseInfoStr(GetInfo(recslice))
	if info, ok := infomap["INFO"]; ok {
		infoscore, _ := strconv.ParseFloat(info, 64)
		return infoscore
	}
	return 1.0
}

func GetInfo(recslice []string) string {
	return recslice[infoIdx]
}

func parseInfoStr(info_str string) map[string]string {
	infomap := make(map[string]string)
	infodata := strings.Split(info_str, ";")
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

func getStrIdx(str string, match_str string) int {
	str_arr := strings.Split(str, ":")
	for i, mstr := range str_arr {
		if mstr == match_str {
			return i
		}
	}
	return -9
}
