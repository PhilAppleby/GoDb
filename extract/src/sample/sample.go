package sample

//---------------------------------------------------
// sample.go:
// Methods for resolving sample lists both by assaytype
// and as a combined string to int map
//---------------------------------------------------

import (
	"fmt"
	"sort"
)

//---------------------------------------------------
//---------------------------------------------------
func MakeSamplesByAssaytype(sample_header_map map[string][]string) (map[string]map[string]int, map[string]map[int]string) {
	sampleNamePosn := make(map[string]map[string]int)
	samplePosnName := make(map[string]map[int]string)

	for at, sample_hdr := range sample_header_map {
		sampleNamePosn[at] = make(map[string]int)
		samplePosnName[at] = make(map[int]string)
		for i, sample := range sample_hdr {
			sampleNamePosn[at][sample] = i
			samplePosnName[at][i] = sample
		}
	}
	return sampleNamePosn, samplePosnName
}
//------------------------------------------------------------------------------
// Condense all the lists by assayttpe into a single map (sample_name to int)
// called by the filemerge code
//------------------------------------------------------------------------------
func GetCombinedSampleMap(samplesByAssayType map[string]map[string]int) map[string]int {
	sample_list := make([]string, 0, len(samplesByAssayType))
	sample_index := make(map[string]int, len(samplesByAssayType))

	for _, posns := range samplesByAssayType {
		for sampleId, _ := range posns {
			sample_list = append(sample_list, sampleId)
		}
	}
	sort.Strings(sample_list)

	i := 0

	for _, samp := range sample_list {
		if _, ok := sample_index[samp]; !ok {
			sample_index[samp] = i
			i++
		}
	}
	return sample_index
}
//------------------------------------------------------------------------------
// Condense all the lists by assayttpe into a single map (sample_name to int)
// goes through each required assaytype and places samples in a master combined
// list of the sample_name is not already there, whilc incrementing the idx
//------------------------------------------------------------------------------
func GetCombinedSampleMapByAssaytypes(samplesByAssayType map[string]map[string]int, assayTypeList []string) map[string]int {
	sample_index := make(map[string]int, len(samplesByAssayType))

	idx := 0
	for _, atype := range assayTypeList {
		if posns, ok := samplesByAssayType[atype]; ok {
			for samp, _ := range posns {
				if _, ok := sample_index[samp]; !ok {
					sample_index[samp] = idx
					idx++
				}
			}
		} else {
			fmt.Printf("##Horrible error %s (assaytype) not found in samples\n", atype)
		}
	}
	return sample_index
}
