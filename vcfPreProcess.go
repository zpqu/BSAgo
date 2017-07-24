// Author: Zhipeng
// Date: 02/06/2017
// This tool is used to split SNVs with multiple alt alleles
// in a vcf file into single alt allele in multiple lines, with each
//line represents a single alt allele.
package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"log"
	"os"
	"regexp"
	"strconv"
	"strings"
	"time"
)

var vcfIn = flag.String("in", "", "Input vcf file")
var vcfOut = flag.String("out", "", "Output vcf file")

// getSNVheader will read header information from vcf file and return a slice type
func getSNVheader(path string) ([]string, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	header, _ := regexp.Compile("^##")
	scanner := bufio.NewScanner(file)
	headerStr := []string{}
	for scanner.Scan() {
		line := scanner.Text()
		if header.MatchString(line) {
			headerStr = append(headerStr, line)
		}
	}
	return headerStr, scanner.Err()
}

// getSNVlong will read SNVs from vcf file and return a string type
func getSNVlong(path string) ([]string, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	header, _ := regexp.Compile("^[^#]")
	scanner := bufio.NewScanner(file)
	SNVlong := []string{}
	for scanner.Scan() {
		line := scanner.Text()
		if header.MatchString(line) {
			SNVlong = append(SNVlong, line)
		}
	}
	return SNVlong, scanner.Err()
}

// getNewInfo function will modify INFO fileld of SNVs in vcf file and return
// new INFO field into a string
func getNewInfo(infoField string, infoIdx int) string {
	infoElts := strings.Split(infoField, "=")
	infoSlice := strings.Split(infoElts[1], ",")

	var afBuffer bytes.Buffer
	afBuffer.WriteString(infoElts[0])
	afBuffer.WriteString("=")
	afBuffer.WriteString(infoSlice[infoIdx])
	return afBuffer.String()
}

// getNewGeno function will split geno information into single alleles
func getNewGeno(genoField string, infoIdx int) string {
	var newGeno string
	if genoField == "." {
		newGeno = "."
	} else {
		genoElts := strings.Split(genoField, ":")
		if genoElts[2] == "." {
			newGeno = "."
		} else {
			dpSlice := strings.Split(genoElts[2], ",")
			refDp := dpSlice[0]
			altDp := dpSlice[infoIdx+1]
			newDpSlice := []string{refDp, altDp}

			aoSlice := strings.Split(genoElts[5], ",")
			newAO := aoSlice[infoIdx]

			qaSlice := strings.Split(genoElts[6], ",")
			newQA := qaSlice[infoIdx]

			var newGTstr string
			if refDp == "0" && altDp == "0" {
				newGeno = "."
			} else {
				switch {
				case refDp == "0":
					newGTstr = "1/1"
				case altDp == "0":
					newGTstr = "0/0"
				default:
					newGTstr = "0/1"
				}

				// Modify GL field
				//glSlice := strings.Split(genoElts[7], ",")
				//newGlStr := strings.Join(glSlice[(3*infoIdx):(3*infoIdx+3)], ",")
				newGlStr := "0,0,0"

				//start loop to modify geno field
				newGenoSlice := []string{}
				for idxGeno, valGeno := range genoElts {
					switch idxGeno {
					case 0:
						newGenoSlice = append(newGenoSlice, newGTstr)
					case 2:
						newGenoSlice = append(newGenoSlice, strings.Join(newDpSlice, ","))
					case 5:
						newGenoSlice = append(newGenoSlice, newAO)
					case 6:
						newGenoSlice = append(newGenoSlice, newQA)
					case 7:
						newGenoSlice = append(newGenoSlice, newGlStr)
					default:
						newGenoSlice = append(newGenoSlice, valGeno)
					}
				} //end loop of modify geno field
				newGeno = strings.Join(newGenoSlice, ":")
			}
		}
	}
	return newGeno
}

// splitAllele function will split vcf line with multiple alt alleles and
// return a slice type
func splitAllele(line string) []string {
	var splitAlleleSlice []string
	lineSlice := strings.Fields(line)
	smSlice := []string{lineSlice[9], lineSlice[10], lineSlice[11], lineSlice[12]}

	var newLineSlice []string
	altSlice := strings.Split(lineSlice[4], ",")
	infoSlice := strings.Split(lineSlice[7], ";")
	//start loop of all alter alleles
	for idxAlt, valAlt := range altSlice {
		newLineSlice = newLineSlice[:0] //empty newLineSlice for each altAllele
		newAltStr := valAlt             //get new alter allele

		newInfoSlice := []string{}
		newInfoSlice = newInfoSlice[:0]
		for _, valInfo := range infoSlice { //start loop to modify info field
			infoTest := regexp.MustCompile(",")
			if infoTest.MatchString(valInfo) {
				newInfoSlice = append(newInfoSlice, getNewInfo(valInfo, idxAlt))
			} else {
				newInfoSlice = append(newInfoSlice, valInfo)
			}
		} //end loop of modify info field
		newInfoStr := strings.Join(newInfoSlice, ";")

		//modify DP of GENO filed
		newGenoSlice := []string{}
		newGenoSlice = newGenoSlice[:0]
		for _, valGeno := range smSlice {
			newGenoSlice = append(newGenoSlice, getNewGeno(valGeno, idxAlt))
		}
		newGenoStr := strings.Join(newGenoSlice, "\t")

		//start make new line slice
		newLineSlice := lineSlice[0:4]
		newLineSlice = append(newLineSlice, newAltStr)
		newLineSlice = append(newLineSlice, lineSlice[5:7]...)
		newLineSlice = append(newLineSlice, newInfoStr)
		newLineSlice = append(newLineSlice, lineSlice[8])
		newLineSlice = append(newLineSlice, newGenoStr)

		splitAlleleSlice = append(splitAlleleSlice, strings.Join(newLineSlice, "\t"))
	} //end loop of all alter alleles
	return splitAlleleSlice
}

// getDPslice function will get all depth information for each alleles
func getDPslice(line string) []int {
	dpSlice := []int{}
	lineField := strings.Fields(line)
	smSlice := []string{lineField[9], lineField[10], lineField[11], lineField[12]}
	for _, valGeno := range smSlice {
		if valGeno == "." {
			dpSlice = append(dpSlice, 0)
		} else {
			genoSlice := strings.Split(valGeno, ":")
			if genoSlice[2] == "." {
				dpSlice = append(dpSlice, 0)
			} else {
				dpElts := strings.Split(genoSlice[2], ",")
				refDp, _ := strconv.Atoi(dpElts[0])
				altDp, _ := strconv.Atoi(dpElts[1])
				dpSlice = append(dpSlice, refDp+altDp)
			}
		}
	}
	return dpSlice
}

// sortSMslice function will get all depth information for each alleles
func sortSMslice(line string) string {
	lineField := strings.Fields(line)
	smSlice := []string{lineField[9], lineField[12], lineField[11], lineField[10]}
	newLineSlice := []string{}
	newLineSlice = append(newLineSlice, lineField[0:9]...)
	newLineSlice = append(newLineSlice, smSlice...)
	return strings.Join(newLineSlice, "\t")
}

//function of write files
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

//main function
func main() {
	flag.Parse()
	fmt.Println("[", time.Now(), "] ", "Program start ...")

	//read header of vcf file
	vcfHeader, err := getSNVheader(*vcfIn)
	if err != nil {
		log.Fatalf("read vcf header: $s", err)
	}

	//read vcf lines from vcf file
	vcfLines, err := getSNVlong(*vcfIn)
	if err != nil {
		log.Fatalf("read vcf file: %s", err)
	}

	fmt.Println("Total number of vcf is: ", len(vcfLines))

	//put header into slice outVcf
	var outVcf []string
	for _, sHeader := range vcfHeader {
		outVcf = append(outVcf, sHeader)
	}

	// append new title for vcf file
	newVcfTitleSlice := []string{
		"#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
		"FILTER", "INFO", "FORMAT", "v3", "CE2_1", "WT", "MT",
	}

	outVcf = append(outVcf, strings.Join(newVcfTitleSlice, "\t"))

	//put vcf lines (split if multiple alt alleles exist) into outVcf
	newVcfLines := []string{}
	for _, sLine := range vcfLines {
		lineField := strings.Fields(sLine)
		mnp, _ := regexp.Compile(",")
		if mnp.MatchString(lineField[4]) {
			newVcfLines = append(newVcfLines, splitAllele(sLine)...)
		} else {
			newVcfLines = append(newVcfLines, sLine)
		}
	}

	fmt.Println("Total number of split VCF is: ", len(newVcfLines))

	i := 0
	for _, line := range newVcfLines {
		dpIntSlice := getDPslice(line)
		if dpIntSlice[0] > 9 && dpIntSlice[1] > 9 && dpIntSlice[2] > 9 && dpIntSlice[3] > 9 {
			newLine := sortSMslice(line)
			outVcf = append(outVcf, newLine)
			i++
		}
	}

	fmt.Println("Total number of passed split VCF is: ", i)

	//write new vcf file
	if err := writeSNVlong(outVcf, *vcfOut); err != nil {
		log.Fatalf("writeLines: $s", err)
	}

	fmt.Println("[", time.Now(), "] ", "Program end ...")
}
