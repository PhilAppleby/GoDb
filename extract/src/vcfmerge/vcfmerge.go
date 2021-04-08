package vcfmerge

//---------------------------------------------------------
// File: vcfmerge.go
// Functions to merge (combine, merge is a misnomer) arrays
// of VCF genotype data
// Author: P Appleby, University of Dundee
//---------------------------------------------------------
import (
	"genometrics"
	"log"
	"sort"
	"strings"
	"variant"
)

// Vcfdata ...
// selected data extracted from a record
//-----------------------------------------------
type Vcfdata struct {
	Assaytype string
	Probidx   int
	Callrate  float64
	Infoscore float64
}

const hdrPrfx = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"

var assayTypeAbbreviation = map[string]string{
	"affy":        "A",
	"illumina":    "I",
	"affy1KG":     "A1",
	"illumina1KG": "I1",
	"broad":       "B",
	"metabo":      "M",
	"exome":       "E",
	"bigtest":     "T",
	"biggertest":  "G",
}

// GetCombinedColumnHeaders ...
//
func GetCombinedColumnHeaders(sampleNameMap map[string]int) (string, []string) {

	posns := make(map[int]string, len(sampleNameMap))
	samplePosns := make([]int, 0, len(sampleNameMap))
	sampleNames := make([]string, 0, len(sampleNameMap))

	for name, posn := range sampleNameMap {
		samplePosns = append(samplePosns, posn)
		posns[posn] = name
	}
	sort.Ints(samplePosns)

	for _, posn := range samplePosns {
		sampleNames = append(sampleNames, posns[posn])
	}
	return hdrPrfx + strings.Join(sampleNames, "\t"), sampleNames
}

// GetColumnHeaders ...
// Get single-string column header record
func GetColumnHeaders(samplePosnMap map[int]string) string {

	samplePosns := make([]int, 0, len(samplePosnMap))
	sampleNames := make([]string, 0, len(samplePosnMap))

	for posn := range samplePosnMap {
		samplePosns = append(samplePosns, posn)
	}
	sort.Ints(samplePosns)

	for _, posn := range samplePosns {
		sampleNames = append(sampleNames, samplePosnMap[posn])
	}
	return hdrPrfx + strings.Join(sampleNames, "\t")
}

// CombineOne ...
// single SNP version : wraps Combine
//------------------------------------------------------------------------------
func CombineOne(vcfset [][]string, vcfdataset []Vcfdata, rsid string,
	sampleNamesByPosn map[string]map[int]string, comboPosns map[string]int,
	comboNames []string, threshold float64, gmetrics *genometrics.AllMetrics) string {

	var rtnString string
	fileRecords := make(chan string, 1)
	Combine(vcfset, vcfdataset, rsid, sampleNamesByPosn, comboPosns, comboNames, threshold, gmetrics, fileRecords)

	close(fileRecords)
	for rec := range fileRecords {
		rtnString = rec
	}
	return rtnString
}

// Combine ...
// version III: build out full results arrays for each assay in the vcfset,
// then process in lockstep to allow comparision of all genotypes for the same
// sample at the same time
//------------------------------------------------------------------------------
func Combine(vcfset [][]string, vcfdataset []Vcfdata, rsid string,
	sampleNamesByPosn map[string]map[int]string, comboPosns map[string]int,
	comboNames []string, threshold float64, gmetrics *genometrics.AllMetrics, recs chan string) {

	var prfx []string
	var sfx []string
	probidx := 1

	comborec := make([]string, len(comboPosns))
	for i := range comborec {
		comborec[i] = "."
	}

	assayrecs := make([][]string, 0, len(vcfset))
	atypeList := make([]string, 0, len(vcfset))
	atypeMap := make(map[string][]string)

	savedVarid := ""
	savedRefAllele := ""
	savedAltAllele := ""

	if len(vcfset) > 0 {
		savedVarid = variant.GetVarid(vcfset[0][1:])
		savedRefAllele, savedAltAllele = variant.GetAlleles(vcfset[0][1:])
	}

	for _, rec := range vcfset {
		atype := rec[0]
		// currentRecord will have slots in the same order as the combined record
		currentRecord := make([]string, len(comboPosns))

		prfx, sfx = variant.GetVCFPrfxSfx(rec[1:])
		atypeList = append(atypeList, atype)
		atypeMap[atype] = prfx
		(*gmetrics).AllGenoCount += len(sfx)
		probidx = variant.GetProbIdx(prfx)
		varid := variant.GetVarid(prfx)
		refAllele, altAllele := variant.GetAlleles(prfx)
		hasAT := variant.HasFmt(prfx, "AT")
		if (varid == savedVarid) && (refAllele == savedRefAllele) && (altAllele == savedAltAllele) {
			for j, elem := range sfx {
				if !hasAT {
					elem = appendAssayAbbrev(elem, atype)
				}
				// this is the crux: map from sample_name at slot j in the assay type record to the position
				// aligned with the combination record
				currentRecord[comboPosns[sampleNamesByPosn[atype][j]]] = elem
			}
			assayrecs = append(assayrecs, currentRecord)
		} else {
			log.Printf("REJ: merge mismatch: %v (%s, %s, %s)\n", prfx, savedVarid, savedRefAllele, savedAltAllele)
		}
	}
	// At this point all "input" genotype data has been captured and is lined up with the comborec
	// Now look at each possible genotype for the comborec
	for i := range comborec {
		genoList := make([]string, 0, len(assayrecs))
		aList := make([]string, 0, len(assayrecs))
		for j, genos := range assayrecs {
			if genos[i] != "" {
				genoList = append(genoList, variant.GetGeno(genos[i], threshold, probidx))
				aList = append(aList, atypeList[j])
			}
		}
		(*gmetrics).UniqueGenoCount++
		if len(genoList) > 1 {
			(*gmetrics).OverlapTestCount++
			if len(genoList) == 2 {
				(*gmetrics).TwoOverlapCount++
			} else {
				(*gmetrics).GtTwoOverlapCount++
			}
			comborec[i] = getBestGeno(genoList, aList, atypeMap, probidx, rsid, gmetrics)
		} else {
			if len(genoList) == 1 {
				comborec[i] = genoList[0]
				if strings.HasPrefix(comborec[i], "./.") {
					(*gmetrics).MissingCount++
				}
			}
		}

	}
	if !variant.HasFmt(prfx, "AT") {
		prfx = variant.AppendToFmt(prfx, "AT")
	}
	// no leading chr zeros
	prfx = variant.NormaliseChromosome(prfx)
	comborec = append(prfx, comborec...)
	recStr := strings.Join(comborec, "\t")
	recs <- recStr
}

//------------------------------------------------------------------------------
// Equality test for genotypes
//------------------------------------------------------------------------------
func areEqual(geno1 string, geno2 string) bool {
	g1 := strings.Split(geno1, ":")
	g2 := strings.Split(geno2, ":")
	if g1[0] == g2[0] {
		return true
	}
	return false
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
func appendAssayAbbrev(geno string, assaytype string) string {
	if v, ok := assayTypeAbbreviation[assaytype]; ok {
		return geno + ":" + v
	} else {
		return geno
	}
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
func resolveGeno(geno1 string, geno2 string, probidx int) string {
	if areEqual(geno1, geno2) {
		return geno1
	}
	mp2, _, _ := variant.MaxProb(geno2, probidx)
	mp1, _, _ := variant.MaxProb(geno1, probidx)
	if mp2 > mp1 {
		return geno2
	}
	return geno1
}

//------------------------------------------------------------------------------
// Get the best genotype based on imputation probability
//------------------------------------------------------------------------------
func getBestGeno(genoList []string, assayList []string, assayTypeMap map[string][]string,
	probidx int, varid string, gmetrics *genometrics.AllMetrics) string {
	currentGeno := ""     // will contain just a genotype, for example "0/0" or "./."
	bestGenoString := "." // will be returned either as "." or as a gentoype string of the form "0/0:0:1,0,0"
	bestProb := 0.0
	bestInfoScore := 0.0

	gcount, mcount := countDiffGenos(genoList)
	(*gmetrics).MismatchCount += (gcount - 1)
	(*gmetrics).MissTestCount += mcount
	for i, geno := range genoList {
		if geno == "." {
			(*gmetrics).NoAssayCount++
			continue
		}
		genodata := strings.Split(geno, ":")
		if genodata[0] != currentGeno {
			prob, _, _ := variant.MaxProb(geno, probidx)
			if prob == bestProb {
				(*gmetrics).SameProbDiffs++
				infoScore := variant.GetInfoScore(assayTypeMap[assayList[i]])
				if infoScore > bestInfoScore {
					currentGeno = genodata[0]
					bestGenoString = geno
					bestProb = prob
					bestInfoScore = infoScore
				}
			} else {
				if bestProb != 0.0 {
					(*gmetrics).DiffProbDiffs++
				}
				if prob > bestProb {
					currentGeno = genodata[0]
					bestGenoString = geno
					bestProb = prob
					bestInfoScore = variant.GetInfoScore(assayTypeMap[assayList[i]])
				}
			}
		}
	}
	bestGenodata := strings.Split(bestGenoString, ":")
	if bestGenodata[0] == "./." {
		(*gmetrics).MissingCount++
	}
	return bestGenoString
}

//------------------------------------------------------------------------------
// Count the number of genotypes in a genotype list
//------------------------------------------------------------------------------
func countDiffGenos(genoList []string) (int, int) {
	var gcount, mcount int
	gmap := make(map[string]bool, len(genoList))
	for _, genoinfo := range genoList {
		geno := strings.Split(genoinfo, ":")
		if _, ok := gmap[geno[0]]; !ok {
			gmap[geno[0]] = true
			gcount++
		}
		if geno[0] == "./." {
			mcount++
		}
	}
	return gcount, mcount
}
