package genometrics

//
// Calculate selected metrics for a VCF record-as-slice
// These include HWEP, CR, MAF
// Note that the use of runParameters is not implemented as yet (TODO)
//
//
import (
	"log"
	"strconv"
	"strings"
	"variant"
)

type AllMetrics struct {
	AllGenoCount      int
	UniqueGenoCount   int
	OverlapTestCount  int
	OverlapSampCount  int
	TwoOverlapCount   int
	GtTwoOverlapCount int
	MismatchCount     int
	DiffProbDiffs     int
	SameProbDiffs     int
	MissTestCount     int
	MissingCount      int
	NoAssayCount      int
}

// RunParameters ...
type RunParameters struct {
	TestNum   int
	MafDelta  float64
	CallRate  float64
	InfoScore float64
}

// Increment ...
// add the values from src to tgt
//-----------------------------------------------------------------------------
func Increment(tgt *AllMetrics, src *AllMetrics) {
	(*tgt).AllGenoCount += (*src).AllGenoCount
	(*tgt).UniqueGenoCount += (*src).UniqueGenoCount
	(*tgt).OverlapTestCount += (*src).OverlapTestCount
	(*tgt).OverlapSampCount += (*src).OverlapSampCount
	(*tgt).TwoOverlapCount += (*src).TwoOverlapCount
	(*tgt).GtTwoOverlapCount += (*src).GtTwoOverlapCount
	(*tgt).MismatchCount += (*src).MismatchCount
	(*tgt).DiffProbDiffs += (*src).DiffProbDiffs
	(*tgt).SameProbDiffs += (*src).SameProbDiffs
	(*tgt).MissTestCount += (*src).MissTestCount
	(*tgt).MissTestCount += (*src).MissTestCount
	(*tgt).MissingCount += (*src).MissingCount
	(*tgt).NoAssayCount += (*src).NoAssayCount
}

// LogMetrics ...
// log metrics output, detail depends on level
func LogMetrics(logLevel int, varid string, snpcount int, msg string, genomet *AllMetrics) {
	errorPct := (float64((*genomet).MismatchCount) / float64((*genomet).OverlapTestCount)) * 100
	if logLevel > 0 {
		log.Printf("%s,%s,SNPCount=%d\n", msg, varid, snpcount)
		log.Printf("%s,%s,ErrPct=%.3f\n", msg, varid, errorPct)
	}
	if logLevel > 1 {
		log.Printf("%s,%s,AllGenos=%d\n", msg, varid, (*genomet).AllGenoCount)
		log.Printf("%s,%s,UniqueGenos=%d\n", msg, varid, (*genomet).UniqueGenoCount)
		log.Printf("%s,%s,Alloverlap=%d\n", msg, varid, (*genomet).OverlapTestCount)
		log.Printf("%s,%s,Two=%d\n", msg, varid, (*genomet).TwoOverlapCount)
		log.Printf("%s,%s,GTTwo=%d\n", msg, varid, (*genomet).GtTwoOverlapCount)
		log.Printf("%s,%s,OverlapGenoDiffs=%d\n", msg, varid, (*genomet).MismatchCount)
		log.Printf("%s,%s,DiffProbDiffs=%d\n", msg, varid, (*genomet).DiffProbDiffs)
		log.Printf("%s,%s,SameProbDiffs=%d\n", msg, varid, (*genomet).SameProbDiffs)
		log.Printf("%s,%s,MissingGenoTested=%d\n", msg, varid, (*genomet).MissTestCount)
		log.Printf("%s,%s,MissingUnresolved=%d\n", msg, varid, (*genomet).MissingCount)
		log.Printf("%s,%s,NoAssay=%d\n", msg, varid, (*genomet).NoAssayCount)
	}
}

// HweExactForRecord ...
// caller passes a string array representing a whole VCF
// record, including prefix
// Return the results for the SNPHWE fn.
func HweExactForRecord(rec []string, threshold float64) float64 {
	homref, homalt, het, _, _, _, _ := getGenotypeCounts(rec, threshold)
	return SNPHWE(het, homref, homalt)
}

// MetricsForRecord ...
// return all SNP metrics
// CR, RAF, AAF, MAF, HWE_P
func MetricsForRecord(rec []string, threshold float64) (float64, float64,
	float64, float64, float64, int, int, int, int, int, int, float64) {
	homref, homalt, het, alln, miss, dot, refPAF := getGenotypeCounts(rec, threshold)
	n := alln - miss
	cr := float64(homref+het+homalt) / float64(alln)
	raf := float64(2*homref+het) / float64(2*n)
	aaf := float64(2*homalt+het) / float64(2*n)
	maf := aaf
	if raf < aaf {
		maf = raf
	}
	obsHomc := homref
	obsHomr := homalt
	if homalt > homref {
		obsHomc = homalt
		obsHomr = homref
	}
	return cr, raf, aaf, maf, SNPHWE(het, homref, homalt), het, obsHomc, obsHomr, n, miss, dot, refPAF
}

// GetRunParams ...
func GetRunParams(testnum string, mafdelta string, callrate string, infoscore string) RunParameters {
	var runParams RunParameters

	runParams.TestNum = 0
	runParams.MafDelta = 0.0
	runParams.CallRate = 0.0
	runParams.InfoScore = 0.0

	if testnum != "" {
		runParams.TestNum, _ = strconv.Atoi(testnum)
	}
	if mafdelta != "" {
		runParams.MafDelta, _ = strconv.ParseFloat(mafdelta, 32)
	}
	if callrate != "" {
		runParams.CallRate, _ = strconv.ParseFloat(callrate, 32)
	}
	if infoscore != "" {
		runParams.InfoScore, _ = strconv.ParseFloat(infoscore, 32)
	}

	return runParams
}

// SNPHWE ...
// Following is translated from the 'C' routine from UMICH
// Original comments:
//
// "This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
// Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
// Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000
//
// Written by Jan Wigginton"
//
func SNPHWE(obsHets int, obsHom1 int, obsHom2 int) float64 {
	obsHomc := obsHom1
	obsHomr := obsHom2

	if obsHom2 > obsHom1 {
		obsHomc = obsHom2
		obsHomr = obsHom1
	}

	rareCopies := 2*obsHomr + obsHets
	genotypes := obsHets + obsHomc + obsHomr

	hetProbs := make([]float64, rareCopies+1, rareCopies+1)

	// start at the midpoint
	mid := rareCopies * (2*genotypes - rareCopies) / (2 * genotypes)
	//fmt.Printf("mid %d, rare %d\n", mid, rareCopies)

	if ((rareCopies & 1) ^ (mid & 1)) != 0 {
		mid++
		//fmt.Printf("add 1 to mid %d\n", mid)
	}

	currHomr := (rareCopies - mid) / 2
	currHomc := genotypes - mid - currHomr

	hetProbs[mid] = 1.0
	sum := 1.0

	for currHets := mid; currHets > 1; currHets -= 2 {
		hetProbs[currHets-2] = hetProbs[currHets] * float64(currHets) * float64(currHets-1.0) / (4.0 * float64(currHomr+1.0) * float64(currHomc+1.0))
		sum += hetProbs[currHets-2]
		// 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
		currHomr++
		currHomc++
	}

	currHomr = (rareCopies - mid) / 2
	currHomc = genotypes - mid - currHomr

	for currHets := mid; currHets <= rareCopies-2; currHets += 2 {
		hetProbs[currHets+2] = hetProbs[currHets] * 4.0 * float64(currHomr) * float64(currHomc) / (float64(currHets+2.0) * float64(currHets+1.0))
		sum += hetProbs[currHets+2]
		// add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
		currHomr -= 1
		currHomc -= 1
	}

	for i, _ := range hetProbs {
		hetProbs[i] /= sum
	}

	pHwe := 0.0

	for i := 0; i <= rareCopies; i++ {
		if hetProbs[i] <= hetProbs[obsHets] {
			pHwe += hetProbs[i]
		}
	}
	if pHwe > 1.0 {
		pHwe = 1.0
	}
	return pHwe

}

func getGenotypeCounts(rec []string, threshold float64) (int, int, int, int, int, int, float64) {
	homr := 0
	homa := 0
	het := 0
	n := 0
	miss := 0
	dot := 0

	prfx, sfx := variant.GetVCFPrfxSfx(rec)
	probidx := variant.GetProbIdx(prfx)
	refPAF := variant.GetRefPanelAF(prfx)

	for _, geno := range sfx {
		if geno != "." {
			n++
			geno = variant.GetGeno(geno, threshold, probidx)
			genoA := strings.Split(geno, ":")
			if genoA[0] == "0/0" {
				homr++
			}
			if genoA[0] == "0/1" {
				het++
			}
			if genoA[0] == "1/0" {
				het++
			}
			if genoA[0] == "1/1" {
				homa++
			}
			if genoA[0] == "./." {
				miss++
			}
		} else {
			dot++
		}
	}

	return homr, homa, het, n, miss, dot, refPAF
}
