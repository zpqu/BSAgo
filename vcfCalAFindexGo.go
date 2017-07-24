package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"regexp"
	"strconv"
	"strings"
)

//declare command arguments
var (
	vcfFile  = flag.String("in", "", "Input vcf file")
	failFile = flag.String("fail", "", "Output file with failed index")
	passFile = flag.String("pass", "", "Output file with passed index")
	skipFile = flag.String("skip", "", "Output file with skiped index")
)

// getSNVlong function reads lines from vcf file and return a slice type
func getSNVlong(path string) ([]string, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	header, _ := regexp.Compile("^[^#]")
	SNVlong := []string{}
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if header.MatchString(line) {
			SNVlong = append(SNVlong, line)
		}
	}
	return SNVlong, scanner.Err()
}

// writeSNVlong function will write vcf slices in to out file
func writeSNVlong(SNVlong []string, path string) error {
	file, err := os.Create(path)
	if err != nil {
		return err
	}
	defer file.Close()

	w := bufio.NewWriter(file)
	for _, line := range SNVlong {
		fmt.Fprintln(w, line)
	}
	return w.Flush()
}

//function of calculate index
func calculateSnvIndex(valLine string) float64 {
	var sampleIndex float64 = 0.0
	if valLine != "." {
		genoSlice := strings.Split(valLine, ":")
		dpSlice := strings.Split(genoSlice[2], ",")
		sampleRef, err := strconv.ParseFloat(dpSlice[0], 64)
		if err == nil {
			sampleAlt, err := strconv.ParseFloat(dpSlice[1], 64)
			if err == nil {
				sampleIndex = sampleAlt / (sampleAlt + sampleRef)
			}
		}
	}
	return sampleIndex
}

func getSnvIndex(lineSlice []string) []float64 {
	result := []float64{}

	var (
		v3Index = 0.0
		ceIndex = 0.0
		wtIndex = 0.0
		mtIndex = 0.0
	)

	for idxLine, valLine := range lineSlice {
		switch idxLine {
		case 9:
			v3Index = calculateSnvIndex(valLine)
		case 10:
			ceIndex = calculateSnvIndex(valLine)
		case 11:
			wtIndex = calculateSnvIndex(valLine)
		case 12:
			mtIndex = calculateSnvIndex(valLine)
		}
	}
	result = append(result, v3Index, ceIndex, wtIndex, mtIndex)
	return result
}

// main function
func main() {
	flag.Parse()

	vcfLines, err := getSNVlong(*vcfFile)
	if err != nil {
		log.Fatalf("read input vcf file: %s", err)
	}

	snvSkipSlice := []string{}
	snvPassSlice := []string{}
	snvFailSlice := []string{}
	indexHeader := []string{
		"#ID", "geno_v3", "geno_CE", "geno_WT", "geno_MT",
		"idx_v3", "idx_CE", "idx_WT", "idx_MT", "deltaIdx_Parent", "deltaIdx_F2",
	}
	snvPassSlice = append(snvPassSlice, strings.Join(indexHeader, "\t"))
	snvFailSlice = append(snvFailSlice, strings.Join(indexHeader, "\t"))

	for _, valSNV := range vcfLines {
		snvSlice := strings.Fields(valSNV)
		if snvSlice[9] != "./." && snvSlice[10] != "./." && snvSlice[11] != "./." && snvSlice[12] != "./." {
			newIdSlice := []string{}
			newSnvSlice := []string{}
			newIdSlice = append(newIdSlice, snvSlice[0:2]...)
			newIdSlice = append(newIdSlice, snvSlice[3:5]...)
			newSnvSlice = append(newSnvSlice, strings.Join(newIdSlice, "_"))
			newSnvSlice = append(newSnvSlice, snvSlice[9], snvSlice[10], snvSlice[11], snvSlice[12])

			vcfIndex := getSnvIndex(snvSlice)
			vcfIndexParent := vcfIndex[0] - vcfIndex[1]
			vcfIndexF2 := 0.0
			if vcfIndexParent < 0.0 {
				vcfIndexF2 = (vcfIndex[2] - vcfIndex[3]) * -1.0
			} else {
				vcfIndexF2 = (vcfIndex[2] - vcfIndex[3]) * 1.0
			}
			vcfIndex = append(vcfIndex, vcfIndexParent, vcfIndexF2)
			vcfIndexStrSlice := []string{}
			for _, valVcfIndex := range vcfIndex {
				vcfIndexStrSlice = append(vcfIndexStrSlice, strconv.FormatFloat(valVcfIndex, 'f', 2, 64))
			}
			newSnvSlice = append(newSnvSlice, vcfIndexStrSlice...)

			if math.Abs(vcfIndexParent) > 0.9 {
				if vcfIndex[2] < 0.3 && vcfIndex[3] < 0.3 {
					snvFailSlice = append(snvFailSlice, strings.Join(newSnvSlice, "\t"))
				} else {
					snvPassSlice = append(snvPassSlice, strings.Join(newSnvSlice, "\t"))
				}
			} else {
				snvFailSlice = append(snvFailSlice, strings.Join(newSnvSlice, "\t"))
			}
		} else {
			snvSkipSlice = append(snvSkipSlice, valSNV)
		}
	}

	//write skip snv file
	if err := writeSNVlong(snvSkipSlice, *skipFile); err != nil {
		log.Fatalf("write skip vcf: $s", err)
	}

	// Write pass snv file
	if err := writeSNVlong(snvPassSlice, *passFile); err != nil {
		log.Fatalf("write pass vcf: $s", err)
	}

	// Write fail snv file
	if err := writeSNVlong(snvFailSlice, *failFile); err != nil {
		log.Fatalf("write fail vcf: $s", err)
	}

}
