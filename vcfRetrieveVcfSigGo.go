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

var (
	vcfFile  = flag.String("vcf", "", "Input vcf file of sorted merged vcfs")
	passFile = flag.String("pass", "", "Input text file including passed vcf info")
	binFile  = flag.String("bin", "", "Input sliding window file with ave idx")
	passOut  = flag.String("passOut", "", "Output vcf file including passed vcf with idx")
	sigOut   = flag.String("sigOut", "", "Output vcf file including only vcfs in sig bins")
	sigBin   = flag.String("sigBin", "", "Output bin file including sig bins")
)

// getSNVheader function will read lines from vcf file and return all header
// lines into a slice type
func getSNVheader(path string) (map[string]string, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	fileformat := regexp.MustCompile("^##fileformat=")
	filterField := regexp.MustCompile("^##FILTER=")
	formatField := regexp.MustCompile("^##FORMAT=")
	infoField := regexp.MustCompile("^##INFO=")
	contigField := regexp.MustCompile("^##contig=")
	referenceField := regexp.MustCompile("^##reference=")

	scanner := bufio.NewScanner(file)
	headerMap := map[string]string{}
	fileformatSlice := []string{}
	filterSlice := []string{}
	formatSlice := []string{}
	infoSlice := []string{}
	contigSlice := []string{}
	referenceSlice := []string{}
	otherSlice := []string{}
	for scanner.Scan() {
		line := scanner.Text()
		switch {
		case fileformat.MatchString(line):
			fileformatSlice = append(fileformatSlice, line)
		case filterField.MatchString(line):
			filterSlice = append(filterSlice, line)
		case formatField.MatchString(line):
			formatSlice = append(formatSlice, line)
		case infoField.MatchString(line):
			infoSlice = append(infoSlice, line)
		case contigField.MatchString(line):
			contigSlice = append(contigSlice, line)
		case referenceField.MatchString(line):
			referenceSlice = append(referenceSlice, line)
		default:
			otherSlice = append(otherSlice, line)
		}
	}
	headerMap["fileformat"] = strings.Join(fileformatSlice, "\n")
	headerMap["filter"] = strings.Join(filterSlice, "\n")
	headerMap["format"] = strings.Join(formatSlice, "\n")
	headerMap["info"] = strings.Join(infoSlice, "\n")
	headerMap["contig"] = strings.Join(contigSlice, "\n")
	headerMap["reference"] = strings.Join(referenceSlice, "\n")
	headerMap["other"] = strings.Join(otherSlice, "\n")

	return headerMap, scanner.Err()
}

// getSNVtitle function will read title Line from vcf and reformat it
func getSNVtitle(path string) ([]string, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	header, _ := regexp.Compile("^#CHROM")
	scanner := bufio.NewScanner(file)
	titleSlice := []string{}
	for scanner.Scan() {
		line := scanner.Text()
		if header.MatchString(line) {
			titleField := strings.Fields(line)[:9]
			titleField = append(titleField, "v3", "CE2_1", "WT", "MT")
			titleSlice = append(titleSlice, strings.Join(titleField, "\t"))
		}
	}
	return titleSlice, scanner.Err()
}

//getSNVmap function will read lines from vcf file and return a map type
func getSNVmap(path string) (map[string]string, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	header, _ := regexp.Compile("^[^#]")
	SNVlong := map[string]string{}
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if header.MatchString(line) {
			lineSlice := strings.Fields(line)

			newInfoSlice := []string{}
			//infoSlice := strings.Split(lineSlice[7], ";")
			//matchInfo := regexp.MustCompile("^(AC=|AF=|AN=|DP=|MQ=)")
			//for _, valSlice := range infoSlice {
			//	if matchInfo.MatchString(valSlice) {
			//		newInfoSlice = append(newInfoSlice, valSlice)
			//	}
			//}

			newLineSlice := []string{}
			newLineSlice = append(newLineSlice, lineSlice[:7]...)
			newLineSlice = append(newLineSlice, strings.Join(newInfoSlice, ";"))
			newLineSlice = append(newLineSlice, lineSlice[8], lineSlice[9], lineSlice[10])
			newLineSlice = append(newLineSlice, lineSlice[11], lineSlice[12])
			keySNPlongSlice := []string{newLineSlice[0], newLineSlice[1], newLineSlice[3], newLineSlice[4]}
			SNVlong[strings.Join(keySNPlongSlice, "_")] = strings.Join(newLineSlice, "\t")
		}
	}
	return SNVlong, scanner.Err()
}

//getSNVkey function will read lines from vcf file and return the keys in order
func getSNVkey(path string) ([]string, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	header, _ := regexp.Compile("^[^#]")
	SNVkey := []string{}
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if header.MatchString(line) {
			lineSlice := strings.Fields(line)

			newLineSlice := lineSlice[:2]
			newLineSlice = append(newLineSlice, lineSlice[3], lineSlice[4])
			SNVkey = append(SNVkey, strings.Join(newLineSlice, "_"))
		}
	}
	return SNVkey, scanner.Err()
}

// getPASSmap function will get pass infomation and put into a map
func getPASSmap(path string) (map[string]string, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	header, _ := regexp.Compile("^[^#]")
	passMap := map[string]string{}
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if header.MatchString(line) {
			lineSlice := strings.Fields(line)
			var buffer bytes.Buffer
			buffer.WriteString("IdxA=")
			buffer.WriteString(lineSlice[5])
			buffer.WriteString(";IdxB=")
			buffer.WriteString(lineSlice[6])
			buffer.WriteString(";IdxF2A=")
			buffer.WriteString(lineSlice[7])
			buffer.WriteString(";IdxF2B=")
			buffer.WriteString(lineSlice[8])
			buffer.WriteString(";DeltaIdxPa=")
			buffer.WriteString(lineSlice[9])
			buffer.WriteString(";DeltaIdxF2=")
			buffer.WriteString(lineSlice[10])
			buffer.WriteString(";P90L=")
			buffer.WriteString(lineSlice[13])
			buffer.WriteString(";P90H=")
			buffer.WriteString(lineSlice[14])
			buffer.WriteString(";P95L=")
			buffer.WriteString(lineSlice[15])
			buffer.WriteString(";P95H=")
			buffer.WriteString(lineSlice[16])
			buffer.WriteString(";P99L=")
			buffer.WriteString(lineSlice[17])
			buffer.WriteString(";P99H=")
			buffer.WriteString(lineSlice[18])

			passMap[lineSlice[0]] = buffer.String()
		}
	}
	return passMap, scanner.Err()
}

// getSigPASSmap function will get pass infomation and put into a map
func getSigPASSmap(path string) (map[string]string, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	header, _ := regexp.Compile("^[^#]")
	sigPassMap := map[string]string{}
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if header.MatchString(line) {
			lineSlice := strings.Fields(line)
			lineFloat := []float64{}
			for _, valLine := range lineSlice[9:] {
				valLineFloat, _ := strconv.ParseFloat(valLine, 64)
				lineFloat = append(lineFloat, valLineFloat)
			}
			//if lineFloat[1] < lineFloat[4] || lineFloat[1] > lineFloat[5] {
			if lineFloat[1] > lineFloat[5] {
				var buffer bytes.Buffer
				buffer.WriteString("IdxA=")
				buffer.WriteString(lineSlice[5])
				buffer.WriteString(";IdxB=")
				buffer.WriteString(lineSlice[6])
				buffer.WriteString(";IdxF2A=")
				buffer.WriteString(lineSlice[7])
				buffer.WriteString(";IdxF2B=")
				buffer.WriteString(lineSlice[8])
				buffer.WriteString(";DeltaIdxPa=")
				buffer.WriteString(lineSlice[9])
				buffer.WriteString(";DeltaIdxF2=")
				buffer.WriteString(lineSlice[10])
				buffer.WriteString(";P90L=")
				buffer.WriteString(lineSlice[13])
				buffer.WriteString(";P90H=")
				buffer.WriteString(lineSlice[14])
				buffer.WriteString(";P95L=")
				buffer.WriteString(lineSlice[15])
				buffer.WriteString(";P95H=")
				buffer.WriteString(lineSlice[16])
				buffer.WriteString(";P99L=")
				buffer.WriteString(lineSlice[17])
				buffer.WriteString(";P99H=")
				buffer.WriteString(lineSlice[18])

				sigPassMap[lineSlice[0]] = buffer.String()
			}
		}
	}
	return sigPassMap, scanner.Err()
}

// getSigBINmap function will get significant slidng windows infomation and put into a slice
func getSigBIN(path string) ([]string, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	header, _ := regexp.Compile("^[^#]")
	sigBIN := []string{}
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if header.MatchString(line) {
			lineSlice := strings.Fields(line)
			lineFloatSlice := []float64{}
			for _, valLine := range lineSlice[4:] {
				valLineFloat, err := strconv.ParseFloat(valLine, 64)
				if err != nil {
					fmt.Println("Error for string converstion in Bin file: ", line)
				}

				lineFloatSlice = append(lineFloatSlice, valLineFloat)
			}

			//sigBIN > 0.1 confidence interval
			if lineFloatSlice[6] > lineFloatSlice[10] {
				sigBIN = append(sigBIN, line)
			}
		}
	}
	return sigBIN, scanner.Err()
}

// getSNVsigBINmap function will check whether vcf is inside sig bins
func getSNVsigBINmap(SNVkey []string, sigBIN []string) map[string]string {
	mapSNVsigBIN := map[string]string{}
	for _, valSNV := range SNVkey {
		valSNVSlice := strings.Split(valSNV, "_")
		valSNVchr := valSNVSlice[0]
		valSNVpos, _ := strconv.Atoi(valSNVSlice[1])
		for _, valBin := range sigBIN {
			valBinSlice := strings.Fields(valBin)
			valBinChr := valBinSlice[0]
			valBinStart, _ := strconv.Atoi(valBinSlice[1])
			valBinEnd, _ := strconv.Atoi(valBinSlice[2])
			valBinF2, _ := strconv.ParseFloat(valBinSlice[10], 64)
			valP90H, _ := strconv.ParseFloat(valBinSlice[14], 64)
			valP95H, _ := strconv.ParseFloat(valBinSlice[16], 64)
			valP99H, _ := strconv.ParseFloat(valBinSlice[18], 64)
			if valSNVchr == valBinChr {
				if valSNVpos >= valBinStart && valSNVpos <= valBinEnd {
					if valBinF2 > valP90H && valBinF2 <= valP95H {
						var buffer bytes.Buffer
						buffer.WriteString("SigBinP90H=")
						buffer.WriteString(valBinSlice[3])
						mapSNVsigBIN[valSNV] = buffer.String()
					} else if valBinF2 > valP95H && valBinF2 <= valP99H {
						var buffer bytes.Buffer
						buffer.WriteString("SigBinP95H=")
						buffer.WriteString(valBinSlice[3])
						mapSNVsigBIN[valSNV] = buffer.String()
					} else if valBinF2 > valP99H {
						var buffer bytes.Buffer
						buffer.WriteString("SigBinP99H=")
						buffer.WriteString(valBinSlice[3])
						mapSNVsigBIN[valSNV] = buffer.String()
					}
				}
			}
		}
	}
	return mapSNVsigBIN
}

//function of write files
func writeLines(lines []string, path string) error {
	file, err := os.Create(path)
	if err != nil {
		return err
	}
	defer file.Close()

	w := bufio.NewWriter(file)
	for _, line := range lines {
		fmt.Fprintln(w, line)
	}
	return w.Flush()
}

//main function
func main() {
	flag.Parse()

	fmt.Println("[", time.Now(), "] ", "Program start ...")

	// get header of vcf files
	vcfHeaderMap, err := getSNVheader(*vcfFile)
	if err != nil {
		log.Fatalf("read vcf header: %s", err)
	}

	// get vcf title
	vcfTitleLines, err := getSNVtitle(*vcfFile)
	if err != nil {
		log.Fatalf("read vcf title: $s", err)
	}

	// get vcf map
	vcfMap, err := getSNVmap(*vcfFile)
	if err != nil {
		log.Fatalf("read vcf map: $s", err)
	}

	// get vcf key
	vcfKeyLines, err := getSNVkey(*vcfFile)
	if err != nil {
		log.Fatalf("read vcf keys: $s", err)
	}
	fmt.Println("[", time.Now(), "] ", "The total number of VCFs is: ", len(vcfKeyLines))

	// get pass map
	passMap, err := getPASSmap(*passFile)
	if err != nil {
		log.Fatalf("read pass map: %s", err)
	}
	fmt.Println("[", time.Now(), "] ", "The total number of pass VCFs is: ", len(passMap))

	//get sig pass map
	sigPassMap, err := getSigPASSmap(*passFile)
	if err != nil {
		log.Fatalf("read sig pass map: %s", err)
	}
	fmt.Println("[", time.Now(), "] ", "The total number of sig pass VCFs is: ", len(sigPassMap))

	//get sig bin lines
	sigBinLines, err := getSigBIN(*binFile)
	if err != nil {
		log.Fatalf("read sliding windows: $s", err)
	}
	fmt.Println("[", time.Now(), "] ", "The total number of sig bins is: ", len(sigBinLines))

	// start merge informaiton
	keySigPassSlice := make([]string, len(sigPassMap))
	i := 0
	for keySigPass := range sigPassMap {
		keySigPassSlice[i] = keySigPass
		i++
	}
	vcfSigBinMap := getSNVsigBINmap(keySigPassSlice, sigBinLines)

	passVcfLines := []string{}
	sigVcfLines := []string{}
	for _, valVcfKey := range vcfKeyLines {
		newVcfSlice := strings.Split(vcfMap[valVcfKey], "\t")
		newInfoSlice := []string{}
		if valPass, ok := passMap[valVcfKey]; ok {
			newInfoSlice = append(newInfoSlice, valPass)
			if _, ok := sigPassMap[valVcfKey]; ok {
				if valSigBin, ok := vcfSigBinMap[valVcfKey]; ok {
					newInfoSlice = append(newInfoSlice, valSigBin)
					sigVcfSlice := newVcfSlice[:7]
					sigVcfSlice = append(sigVcfSlice, strings.Join(newInfoSlice, ";"))
					sigVcfSlice = append(sigVcfSlice, newVcfSlice[8:]...)
					sigVcfLines = append(sigVcfLines, strings.Join(sigVcfSlice, "\t"))
				}
			}
			passVcfSlice := newVcfSlice[:7]
			passVcfSlice = append(passVcfSlice, strings.Join(newInfoSlice, ";"))
			passVcfSlice = append(passVcfSlice, newVcfSlice[8:]...)
			passVcfLines = append(passVcfLines, strings.Join(passVcfSlice, "\t"))
		}
	}

	//write sig bins
	headerBinSlice := []string{"#CHR", "START", "END", "mid_pos", "num_SNVs",
		"AveIdx_v3", "AveIdx_ce", "AveIdx_wt", "AveIdx_mt",
		"AveIdx_parent", "AveIdx_F2", "AveDp_wt", "AveDp_mt",
		"Ave_p90L", "Ave_p90H", "Ave_p95L", "Ave_p95H", "Ave_p99L", "Ave_p99H",
	}

	newSigBinLines := []string{strings.Join(headerBinSlice, "\t")}
	newSigBinLines = append(newSigBinLines, sigBinLines...)
	if err := writeLines(newSigBinLines, *sigBin); err != nil {
		log.Fatalf("Write sig bin file: $s", err)
	}

	//get new header Slice
	newVcfHeader := []string{}
	newVcfHeaderInfo := []string{"##INFO=<ID=IdxA,Number=1,Type=Float,Description=\"SNP index of parent A\">",
		"##INFO=<ID=IdxB,Number=1,Type=Float,Description=\"SNP index of parent B\">",
		"##INFO=<ID=IdxF2A,Number=1,Type=Float,Description=\"SNP index of F2 A\">",
		"##INFO=<ID=IdxF2B,Number=1,Type=Float,Description=\"SNP index of F2 B\">",
		"##INFO=<ID=DeltaIdxPa,Number=1,Type=Float,Description=\"Delta SNP index of parents, A - B\">",
		"##INFO=<ID=DeltaIdxF2,Number=1,Type=Float,Description=\"Delta SNP index of F2, A - B\">",
		"##INFO=<ID=P90L,Number=1,Type=Float,Description=\"95% confidence interval of low delta SNP index of parents, A - B\">",
		"##INFO=<ID=P90H,Number=1,Type=Float,Description=\"95% confidence interval of high delta SNP index of parents, A - B\">",
		"##INFO=<ID=P95L,Number=1,Type=Float,Description=\"99% confidence interval of low delta SNP index of parents, A - B\">",
		"##INFO=<ID=P95H,Number=1,Type=Float,Description=\"99% confidence interval of high delta SNP index of parents, A - B\">",
		"##INFO=<ID=P99L,Number=1,Type=Float,Description=\"99% confidence interval of low delta SNP index of parents, A - B\">",
		"##INFO=<ID=P99H,Number=1,Type=Float,Description=\"99% confidence interval of high delta SNP index of parents, A - B\">",
		"##INFO=<ID=SigBinP90H,Number=1,Type=String,Description=\"F2 index of vcf is in  90% < CI <= 95%\">",
		"##INFO=<ID=SigBinP95H,Number=1,Type=String,Description=\"F2 index of vcf is in  95% < CI <= 99%\">",
		"##INFO=<ID=SigBinP99H,Number=1,Type=String,Description=\"F2 index of vcf is in  CI > 99%\">",
	}
	newVcfHeader = append(newVcfHeader, vcfHeaderMap["fileformat"])
	//newVcfHeader = append(newVcfHeader, vcfHeaderMap["filter"])
	newVcfHeader = append(newVcfHeader, vcfHeaderMap["format"])
	newVcfHeader = append(newVcfHeader, newVcfHeaderInfo...)
	newVcfHeader = append(newVcfHeader, vcfHeaderMap["contig"])
	newVcfHeader = append(newVcfHeader, vcfHeaderMap["reference"])

	//write PASS vcfs
	newPassLines := []string{}
	newPassLines = append(newPassLines, newVcfHeader...)
	newPassLines = append(newPassLines, vcfTitleLines...)
	newPassLines = append(newPassLines, passVcfLines...)

	if err := writeLines(newPassLines, *passOut); err != nil {
		log.Fatalf("Write pass vcf file: $s", err)
	}

	//write sig PASS vcfs
	newSigPassLines := []string{}
	newSigPassLines = append(newSigPassLines, newVcfHeader...)
	newSigPassLines = append(newSigPassLines, vcfTitleLines...)
	newSigPassLines = append(newSigPassLines, sigVcfLines...)
	fmt.Println("[", time.Now(), "] ", "The total number of sig pass in bin VCFs is: ", len(sigVcfLines))

	if err := writeLines(newSigPassLines, *sigOut); err != nil {
		log.Fatalf("Write sig pass vcf file: $s", err)
	}

	fmt.Println("[", time.Now(), "] ", "Program end ...")
}
