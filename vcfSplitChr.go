// Author: Zhipeng
// Date: 02/06/2017
// This tool is used to split SNVs with multiple alt alleles
// in a vcf file into single alt allele in multiple lines, with each
//line represents a single alt allele.
package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"regexp"
	"strings"
	"time"
)

var vcfIn = flag.String("in", "", "Input vcf file")
var chr1 = flag.String("chr1", "", "Output vcf file in chr1")
var chr2 = flag.String("chr2", "", "Output vcf file in chr2")
var chr3 = flag.String("chr3", "", "Output vcf file in chr3")
var chr4 = flag.String("chr4", "", "Output vcf file in chr4")
var chr5 = flag.String("chr5", "", "Output vcf file in chr5")
var chr6 = flag.String("chr6", "", "Output vcf file in chr6")
var chr7 = flag.String("chr7", "", "Output vcf file in chr7")
var chr8 = flag.String("chr8", "", "Output vcf file in chr8")
var chr9 = flag.String("chr9", "", "Output vcf file in chr9")
var chr10 = flag.String("chr10", "", "Output vcf file in chr10")
var chrUn = flag.String("chrUn", "", "Output vcf file in chrUn")

// getSNVheader will read header information from vcf file and return a slice type
func getSNVheader(path string) ([]string, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	header, _ := regexp.Compile("^#")
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
	vcfChr1 := []string{}
	vcfChr2 := []string{}
	vcfChr3 := []string{}
	vcfChr4 := []string{}
	vcfChr5 := []string{}
	vcfChr6 := []string{}
	vcfChr7 := []string{}
	vcfChr8 := []string{}
	vcfChr9 := []string{}
	vcfChr10 := []string{}
	vcfChrUn := []string{}
	vcfChr1 = append(vcfChr1, vcfHeader...)
	vcfChr2 = append(vcfChr2, vcfHeader...)
	vcfChr3 = append(vcfChr3, vcfHeader...)
	vcfChr4 = append(vcfChr4, vcfHeader...)
	vcfChr5 = append(vcfChr5, vcfHeader...)
	vcfChr6 = append(vcfChr6, vcfHeader...)
	vcfChr7 = append(vcfChr7, vcfHeader...)
	vcfChr8 = append(vcfChr8, vcfHeader...)
	vcfChr9 = append(vcfChr9, vcfHeader...)
	vcfChr10 = append(vcfChr10, vcfHeader...)

	//put vcf lines (split if multiple alt alleles exist) into outVcf
	for _, sLine := range vcfLines {
		lineField := strings.Fields(sLine)
		switch lineField[0] {
		case "1":
			vcfChr1 = append(vcfChr1, sLine)
		case "2":
			vcfChr2 = append(vcfChr2, sLine)
		case "3":
			vcfChr3 = append(vcfChr3, sLine)
		case "4":
			vcfChr4 = append(vcfChr4, sLine)
		case "5":
			vcfChr5 = append(vcfChr5, sLine)
		case "6":
			vcfChr6 = append(vcfChr6, sLine)
		case "7":
			vcfChr7 = append(vcfChr7, sLine)
		case "8":
			vcfChr8 = append(vcfChr8, sLine)
		case "9":
			vcfChr9 = append(vcfChr9, sLine)
		case "10":
			vcfChr10 = append(vcfChr10, sLine)
		default:
			vcfChrUn = append(vcfChrUn, sLine)
		}
	}

	fmt.Println("Total number of VCF in Chr1 is: ", len(vcfChr1))
	fmt.Println("Total number of VCF in Chr2 is: ", len(vcfChr2))
	fmt.Println("Total number of VCF in Chr3 is: ", len(vcfChr3))
	fmt.Println("Total number of VCF in Chr4 is: ", len(vcfChr4))
	fmt.Println("Total number of VCF in Chr5 is: ", len(vcfChr5))
	fmt.Println("Total number of VCF in Chr6 is: ", len(vcfChr6))
	fmt.Println("Total number of VCF in Chr7 is: ", len(vcfChr7))
	fmt.Println("Total number of VCF in Chr8 is: ", len(vcfChr8))
	fmt.Println("Total number of VCF in Chr9 is: ", len(vcfChr9))
	fmt.Println("Total number of VCF in Chr10 is: ", len(vcfChr10))
	fmt.Println("Total number of VCF in ChrUn is: ", len(vcfChrUn))

	//write new vcf file
	if err := writeSNVlong(vcfChr1, *chr1); err != nil {
		log.Fatalf("write vcf in chr1: $s", err)
	}

	if err := writeSNVlong(vcfChr2, *chr2); err != nil {
		log.Fatalf("write vcf in chr2: $s", err)
	}

	if err := writeSNVlong(vcfChr3, *chr3); err != nil {
		log.Fatalf("write vcf in chr3: $s", err)
	}

	if err := writeSNVlong(vcfChr4, *chr4); err != nil {
		log.Fatalf("write vcf in chr4: $s", err)
	}

	if err := writeSNVlong(vcfChr5, *chr5); err != nil {
		log.Fatalf("write vcf in chr5: $s", err)
	}

	if err := writeSNVlong(vcfChr6, *chr6); err != nil {
		log.Fatalf("write vcf in chr6: $s", err)
	}

	if err := writeSNVlong(vcfChr7, *chr7); err != nil {
		log.Fatalf("write vcf in chr7: $s", err)
	}

	if err := writeSNVlong(vcfChr8, *chr8); err != nil {
		log.Fatalf("write vcf in chr8: $s", err)
	}

	if err := writeSNVlong(vcfChr9, *chr9); err != nil {
		log.Fatalf("write vcf in chr9: $s", err)
	}

	if err := writeSNVlong(vcfChr10, *chr10); err != nil {
		log.Fatalf("write vcf in chr10: $s", err)
	}

	if err := writeSNVlong(vcfChrUn, *chrUn); err != nil {
		log.Fatalf("write vcf in chrUn: $s", err)
	}

	fmt.Println("[", time.Now(), "] ", "Program end ...")
}
