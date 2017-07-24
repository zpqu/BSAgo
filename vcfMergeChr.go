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
	"time"
)

var vcfOut = flag.String("out", "", "output vcf file")
var chr1 = flag.String("chr1", "", "input vcf file in chr1")
var chr2 = flag.String("chr2", "", "input vcf file in chr2")
var chr3 = flag.String("chr3", "", "input vcf file in chr3")
var chr4 = flag.String("chr4", "", "input vcf file in chr4")
var chr5 = flag.String("chr5", "", "input vcf file in chr5")
var chr6 = flag.String("chr6", "", "input vcf file in chr6")
var chr7 = flag.String("chr7", "", "input vcf file in chr7")
var chr8 = flag.String("chr8", "", "input vcf file in chr8")
var chr9 = flag.String("chr9", "", "input vcf file in chr9")
var chr10 = flag.String("chr10", "", "input vcf file in chr10")

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
	vcfHeader, err := getSNVheader(*chr1)
	if err != nil {
		log.Fatalf("read vcf header: $s", err)
	}

	//read vcf lines from vcf file
	vcfChr1, err := getSNVlong(*chr1)
	if err != nil {
		log.Fatalf("read vcf chr1 file: %s", err)
	}

	vcfChr2, err := getSNVlong(*chr2)
	if err != nil {
		log.Fatalf("read vcf chr2 file: %s", err)
	}

	vcfChr3, err := getSNVlong(*chr3)
	if err != nil {
		log.Fatalf("read vcf chr3 file: %s", err)
	}

	vcfChr4, err := getSNVlong(*chr4)
	if err != nil {
		log.Fatalf("read vcf chr4 file: %s", err)
	}

	vcfChr5, err := getSNVlong(*chr5)
	if err != nil {
		log.Fatalf("read vcf chr5 file: %s", err)
	}

	vcfChr6, err := getSNVlong(*chr6)
	if err != nil {
		log.Fatalf("read vcf chr6 file: %s", err)
	}

	vcfChr7, err := getSNVlong(*chr7)
	if err != nil {
		log.Fatalf("read vcf chr7 file: %s", err)
	}

	vcfChr8, err := getSNVlong(*chr8)
	if err != nil {
		log.Fatalf("read vcf chr8 file: %s", err)
	}

	vcfChr9, err := getSNVlong(*chr9)
	if err != nil {
		log.Fatalf("read vcf chr9 file: %s", err)
	}

	vcfChr10, err := getSNVlong(*chr10)
	if err != nil {
		log.Fatalf("read vcf chr10 file: %s", err)
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

	newVcfLine := []string{}
	newVcfLine = append(newVcfLine, vcfHeader...)
	newVcfLine = append(newVcfLine, vcfChr1...)
	newVcfLine = append(newVcfLine, vcfChr2...)
	newVcfLine = append(newVcfLine, vcfChr3...)
	newVcfLine = append(newVcfLine, vcfChr4...)
	newVcfLine = append(newVcfLine, vcfChr5...)
	newVcfLine = append(newVcfLine, vcfChr6...)
	newVcfLine = append(newVcfLine, vcfChr7...)
	newVcfLine = append(newVcfLine, vcfChr8...)
	newVcfLine = append(newVcfLine, vcfChr9...)
	newVcfLine = append(newVcfLine, vcfChr10...)

	fmt.Println("Total number of vcf is: ", len(newVcfLine))

	//write new vcf file
	if err := writeSNVlong(newVcfLine, *vcfOut); err != nil {
		log.Fatalf("write vcf: $s", err)
	}

	fmt.Println("[", time.Now(), "] ", "Program end ...")
}
