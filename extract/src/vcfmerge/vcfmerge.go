package vcfmerge

//---------------------------------------------------------
// File: vcfmerge.go
// Functions to merge (combine, merge is a misnomer) slices
// of VCF genotype data
// Author: P Appleby, University of Dundee
//---------------------------------------------------------
import (
	//"fmt"
	"genometrics"
	"log"
	"sort"
	"strings"
	"variant"
)

//-----------------------------------------------
// Vcfdata: data extracted from a record
//-----------------------------------------------
type Vcfdata struct {
	Assaytype string
	Probidx   int
	Callrate  float64
	Infoscore float64
}

const hdr_prfx = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"

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

func GetCombinedColumnHeaders(sample_name_map map[string]int) (string, []string) {

	posns := make(map[int]string, len(sample_name_map))
	sample_posns := make([]int, 0, len(sample_name_map))
	sample_names := make([]string, 0, len(sample_name_map))

	for name, posn := range sample_name_map {
		sample_posns = append(sample_posns, posn)
		posns[posn] = name
	}
	sort.Ints(sample_posns)

	for _, posn := range sample_posns {
		sample_names = append(sample_names, posns[posn])
	}
	return hdr_prfx + strings.Join(sample_names, "\t"), sample_names
}

func GetColumnHeaders(sample_posn_map map[int]string) string {

	sample_posns := make([]int, 0, len(sample_posn_map))
	sample_names := make([]string, 0, len(sample_posn_map))

	for posn, _ := range sample_posn_map {
		sample_posns = append(sample_posns, posn)
	}
	sort.Ints(sample_posns)

	for _, posn := range sample_posns {
		sample_names = append(sample_names, sample_posn_map[posn])
	}
	return hdr_prfx + strings.Join(sample_names, "\t")
}
//------------------------------------------------------------------------------
// Merge version II: build out full results arrays for each assay in the vcfset,
// then process in lockstep to allow comparision of all genotypes for the same
// sample at the same time
//------------------------------------------------------------------------------
func Mergeslices_full(vcfset [][]string, vcfdataset []Vcfdata, rsid string,
	sample_names_by_posn map[string]map[int]string, combo_posns map[string]int,
	combo_names []string, threshold float64, gmetrics *genometrics.AllMetrics) []string {

	var prfx []string
	var sfx []string
	probidx := 1

	comborec := make([]string, len(combo_posns))
	for i, _ := range comborec {
		comborec[i] = "."
	}

	assayrecs := make([][]string, 0, len(vcfset))
	atypes := make([]string, 0, len(vcfset))

	savedVarid := ""
	savedRefAllele := ""
	savedAltAllele := ""

	if len(vcfset) > 0 {
		savedVarid = variant.GetVarid(vcfset[0][1:])
		savedRefAllele, savedAltAllele = variant.GetAlleles(vcfset[0][1:])
	}

	for _, rec := range vcfset {
		atype := rec[0]
		currec := make([]string, len(combo_posns))
		atypes = append(atypes, atype)
		prfx, sfx = variant.GetVCFPrfx_Sfx(rec[1:])
		probidx = variant.GetProbidx(prfx)
		varid := variant.GetVarid(prfx)
		refAllele, altAllele := variant.GetAlleles(prfx)
		if (varid == savedVarid) && (refAllele == savedRefAllele) && (altAllele == savedAltAllele) {
			for j, elem := range sfx {
        if !variant.HasFmt(prfx, "AT") {
			    currec[combo_posns[sample_names_by_posn[atype][j]]] = appendAssayAbbrev(elem, atype)
			  }
      }
			assayrecs = append(assayrecs, currec)
		} else {
			log.Printf("REJ: merge mismatch: %v (%s, %s, %s)\n", prfx, savedVarid, savedRefAllele, savedAltAllele)
		}
	}
	// At this point all "input" genotype data has been captured - now
	// Look at each possible genotype for the comborec
	for i, _ := range comborec {
		geno_list := make([]string, 0, len(assayrecs))
		for _, genos := range assayrecs {
			if genos[i] != "" {
				geno_list = append(geno_list, variant.Get_geno(genos[i], threshold, probidx))
			}
		}
		//fmt.Printf("%s: %s\n", combo_names[i], geno_list)
		(*gmetrics).AllGenoCount += len(geno_list)
		(*gmetrics).UniqueGenoCount += 1
		if len(geno_list) > 1 {
			(*gmetrics).OverlapTestCount++
			if len(geno_list) == 2 {
				//fmt.Printf("OVERLAP_TWO %s:%s - %s\n", rsid, combo_names[i], geno_list)
				(*gmetrics).TwoOverlapCount++
			} else {
				//fmt.Printf("OVERLAP_GT2 %s:%s - %s\n", rsid, combo_names[i], geno_list)
				(*gmetrics).GtTwoOverlapCount++
			}
			comborec[i] = get_best_geno(geno_list, probidx, rsid, gmetrics)
		} else {
			if len(geno_list) == 1 {
				comborec[i] = geno_list[0]
			}
		}
	}
	prfx = variant.AppendToFmt(prfx, "AT")
	// no leading chr zeros
	prfx = variant.NormaliseChromosome(prfx)
	return append(prfx, comborec...)
}
//------------------------------------------------------------------------------
// Merge version III: build out full results arrays for each assay in the vcfset,
// then process in lockstep to allow comparision of all genotypes for the same
// sample at the same time
//------------------------------------------------------------------------------
func Mergeslices(vcfset [][]string, vcfdataset []Vcfdata, rsid string,
	sample_names_by_posn map[string]map[int]string, combo_posns map[string]int,
	combo_names []string, threshold float64, gmetrics *genometrics.AllMetrics, recs chan string) {

	var prfx []string
	var sfx []string
	probidx := 1

	comborec := make([]string, len(combo_posns))
	for i, _ := range comborec {
		comborec[i] = "."
	}

	assayrecs := make([][]string, 0, len(vcfset))
	atypes := make([]string, 0, len(vcfset))

	savedVarid := ""
	savedRefAllele := ""
	savedAltAllele := ""

	if len(vcfset) > 0 {
		savedVarid = variant.GetVarid(vcfset[0][1:])
		savedRefAllele, savedAltAllele = variant.GetAlleles(vcfset[0][1:])
	}

	for _, rec := range vcfset {
		atype := rec[0]
		currec := make([]string, len(combo_posns))
		atypes = append(atypes, atype)
		prfx, sfx = variant.GetVCFPrfx_Sfx(rec[1:])
		probidx = variant.GetProbidx(prfx)
		varid := variant.GetVarid(prfx)
		refAllele, altAllele := variant.GetAlleles(prfx)
		if (varid == savedVarid) && (refAllele == savedRefAllele) && (altAllele == savedAltAllele) {
			for j, elem := range sfx {
        if !variant.HasFmt(prfx, "AT") {
			    currec[combo_posns[sample_names_by_posn[atype][j]]] = appendAssayAbbrev(elem, atype)
			  }
      }
			assayrecs = append(assayrecs, currec)
		} else {
			log.Printf("REJ: merge mismatch: %v (%s, %s, %s)\n", prfx, savedVarid, savedRefAllele, savedAltAllele)
		}
	}
	// At this point all "input" genotype data has been captured - now
	// Look at each possible genotype for the comborec
	for i, _ := range comborec {
		geno_list := make([]string, 0, len(assayrecs))
		for _, genos := range assayrecs {
			if genos[i] != "" {
				geno_list = append(geno_list, variant.Get_geno(genos[i], threshold, probidx))
			}
		}
		//fmt.Printf("%s: %s\n", combo_names[i], geno_list)
		(*gmetrics).AllGenoCount += len(geno_list)
		(*gmetrics).UniqueGenoCount += 1
		if len(geno_list) > 1 {
			(*gmetrics).OverlapTestCount++
			if len(geno_list) == 2 {
				//fmt.Printf("OVERLAP_TWO %s:%s - %s\n", rsid, combo_names[i], geno_list)
				(*gmetrics).TwoOverlapCount++
			} else {
				//fmt.Printf("OVERLAP_GT2 %s:%s - %s\n", rsid, combo_names[i], geno_list)
				(*gmetrics).GtTwoOverlapCount++
			}
			comborec[i] = get_best_geno(geno_list, probidx, rsid, gmetrics)
		} else {
			if len(geno_list) == 1 {
				comborec[i] = geno_list[0]
			}
		}
	}
  if !variant.HasFmt(prfx, "AT") {
	  prfx = variant.AppendToFmt(prfx, "AT")
  }
	// no leading chr zeros
	prfx = variant.NormaliseChromosome(prfx)
  comborec = append(prfx, comborec...)
  rec_str := strings.Join(comborec, "\t")
	recs <- rec_str
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
//------------------------------------------------------------------------------
func get_best_geno(geno_list []string, probidx int, varid string, gmetrics *genometrics.AllMetrics) string {
	rgeno := ""
	bgeno := "."
	best_prob := 0.0

	for _, geno := range geno_list {
    if geno == "." {
		  (*gmetrics).NoAssayCount += 1
      continue
    }
		genodata := strings.Split(geno, ":")
		if genodata[0] != rgeno {
			if rgeno != "" {
				(*gmetrics).MismatchCount += 1
			}
			if genodata[0] == "./." {
				(*gmetrics).MissTestCount += 1
			}
			prob, _, _ := variant.MaxProb(geno, probidx)
			if prob > best_prob {
				rgeno = genodata[0]
				bgeno = geno
				best_prob = prob
			}
		}
	}
	bgenodata := strings.Split(bgeno, ":")
	if (bgeno == ".") || (bgenodata[0] == "./.") {
		(*gmetrics).MissingCount += 1
	}
	return bgeno
}
