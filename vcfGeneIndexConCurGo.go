package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"log"
	"os"
	"regexp"
	"runtime"
	"strconv"
	"strings"
	"time"
)

// declare command arguments
var (
	bedFile  = flag.String("bed", "", "Input chromosome file")
	vcfFile  = flag.String("vcf", "", "Input vcf file")
	outFile  = flag.String("out", "", "Output file with index in sliding windows")
	cpus     = flag.Int("c", 1, "Number of working CPUs")
	upSize   = flag.Int("up", 0, "Up stream flanking region, default 0")
	downSize = flag.Int("down", 0, "Down stream flanking region, default 0")
)

// getLines function reads lines from vcf file and return a slice,
// with each element representing each line
func getLines(path string) ([]string, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	header, _ := regexp.Compile("^[^#]")
	lines := []string{}
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if header.MatchString(line) {
			lines = append(lines, line)
		}
	}
	return lines, scanner.Err()
}

// getSNVlong function will read from vcf file and return a map type,
// with chromosome and position as key, and index fields as value
func getSNVlong(path string) (map[string]string, error) {
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
			idSlice := strings.Split(lineSlice[0], "_")
			idKey := strings.Join(idSlice[:2], "_")
			idVal := strings.Join(lineSlice[5:], "\t")
			SNVlong[idKey] = idVal
		}
	}
	return SNVlong, scanner.Err()
}

// writeLines function will write lines in a slice in to out file
func writeLines(SNVlong []string, path string) error {
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

// calculateAveIndex function will calculate avergege index
func calculateAveIndex(binLine []int, vcfLine map[string]string) string {
	indexSlice := []string{}
	v3IndexSum := 0.0
	ceIndexSum := 0.0
	wtIndexSum := 0.0
	mtIndexSum := 0.0
	paIndexSum := 0.0
	f2IndexSum := 0.0
	wtDpSum := 0.0
	mtDpSum := 0.0
	p90LSum := 0.0
	p90HSum := 0.0
	p95LSum := 0.0
	p95HSum := 0.0
	p99LSum := 0.0
	p99HSum := 0.0

	numVcf := 0
	for j := binLine[1]; j <= binLine[2]; j++ {
		var posBuffer bytes.Buffer
		posBuffer.WriteString(strconv.Itoa(binLine[0]))
		posBuffer.WriteString("_")
		posBuffer.WriteString(strconv.Itoa(j))
		posKey := posBuffer.String()
		if posVal, ok := vcfLine[posKey]; ok {
			posValSlice := strings.Fields(posVal)
			v3VcfIndex, _ := strconv.ParseFloat(posValSlice[0], 64)
			ceVcfIndex, _ := strconv.ParseFloat(posValSlice[1], 64)
			wtVcfIndex, _ := strconv.ParseFloat(posValSlice[2], 64)
			mtVcfIndex, _ := strconv.ParseFloat(posValSlice[3], 64)
			paVcfIndex, _ := strconv.ParseFloat(posValSlice[4], 64)
			f2VcfIndex, _ := strconv.ParseFloat(posValSlice[5], 64)
			wtVcfDp, _ := strconv.ParseFloat(posValSlice[6], 64)
			mtVcfDp, _ := strconv.ParseFloat(posValSlice[7], 64)
			p90LVcfDp, _ := strconv.ParseFloat(posValSlice[8], 64)
			p90HVcfDp, _ := strconv.ParseFloat(posValSlice[9], 64)
			p95LVcfDp, _ := strconv.ParseFloat(posValSlice[10], 64)
			p95HVcfDp, _ := strconv.ParseFloat(posValSlice[11], 64)
			p99LVcfDp, _ := strconv.ParseFloat(posValSlice[12], 64)
			p99HVcfDp, _ := strconv.ParseFloat(posValSlice[13], 64)

			v3IndexSum += v3VcfIndex
			ceIndexSum += ceVcfIndex
			wtIndexSum += wtVcfIndex
			mtIndexSum += mtVcfIndex
			paIndexSum += paVcfIndex
			f2IndexSum += f2VcfIndex
			wtDpSum += wtVcfDp
			mtDpSum += mtVcfDp
			p90LSum += p90LVcfDp
			p90HSum += p90HVcfDp
			p95LSum += p95LVcfDp
			p95HSum += p95HVcfDp
			p99LSum += p99LVcfDp
			p99HSum += p99HVcfDp

			numVcf++
		}
	}
	v3IndexAve := 0.0
	ceIndexAve := 0.0
	wtIndexAve := 0.0
	mtIndexAve := 0.0
	paIndexAve := 0.0
	f2IndexAve := 0.0
	wtDpAve := 0.0
	mtDpAVe := 0.0
	p90LAve := 0.0
	p90HAve := 0.0
	p95LAve := 0.0
	p95HAve := 0.0
	p99LAve := 0.0
	p99HAve := 0.0

	if numVcf > 3 {
		v3IndexAve = v3IndexSum / float64(numVcf)
		ceIndexAve = ceIndexSum / float64(numVcf)
		wtIndexAve = wtIndexSum / float64(numVcf)
		mtIndexAve = mtIndexSum / float64(numVcf)
		wtIndexAve = wtIndexSum / float64(numVcf)
		paIndexAve = paIndexSum / float64(numVcf)
		f2IndexAve = f2IndexSum / float64(numVcf)
		wtDpAve = wtDpSum / float64(numVcf)
		mtDpAVe = mtDpSum / float64(numVcf)
		p90LAve = p90LSum / float64(numVcf)
		p90HAve = p90HSum / float64(numVcf)
		p95LAve = p95LSum / float64(numVcf)
		p95HAve = p95HSum / float64(numVcf)
		p99LAve = p99LSum / float64(numVcf)
		p99HAve = p99HSum / float64(numVcf)
	}

	indexSlice = append(indexSlice, strconv.Itoa(binLine[0]))
	indexSlice = append(indexSlice, strconv.Itoa(binLine[1]))
	indexSlice = append(indexSlice, strconv.Itoa(binLine[2]))
	indexSlice = append(indexSlice, strconv.Itoa(numVcf))
	indexSlice = append(indexSlice, strconv.FormatFloat(v3IndexAve, 'f', 2, 64))
	indexSlice = append(indexSlice, strconv.FormatFloat(ceIndexAve, 'f', 2, 64))
	indexSlice = append(indexSlice, strconv.FormatFloat(wtIndexAve, 'f', 2, 64))
	indexSlice = append(indexSlice, strconv.FormatFloat(mtIndexAve, 'f', 2, 64))
	indexSlice = append(indexSlice, strconv.FormatFloat(paIndexAve, 'f', 2, 64))
	indexSlice = append(indexSlice, strconv.FormatFloat(f2IndexAve, 'f', 2, 64))
	indexSlice = append(indexSlice, strconv.FormatFloat(wtDpAve, 'f', 2, 64))
	indexSlice = append(indexSlice, strconv.FormatFloat(mtDpAVe, 'f', 2, 64))
	indexSlice = append(indexSlice, strconv.FormatFloat(p90LAve, 'f', 2, 64))
	indexSlice = append(indexSlice, strconv.FormatFloat(p90HAve, 'f', 2, 64))
	indexSlice = append(indexSlice, strconv.FormatFloat(p95LAve, 'f', 2, 64))
	indexSlice = append(indexSlice, strconv.FormatFloat(p95HAve, 'f', 2, 64))
	indexSlice = append(indexSlice, strconv.FormatFloat(p99LAve, 'f', 2, 64))
	indexSlice = append(indexSlice, strconv.FormatFloat(p99HAve, 'f', 2, 64))

	return strings.Join(indexSlice, "\t")
}

//get flanking regions
func getGeneBin(binLine []string, upSize int, downSize int) []string {
	newGeneFlank := []string{}
	for _, valGeneBin := range binLine {
		binSlice := strings.Fields(valGeneBin)
		chrPos, _ := strconv.Atoi(binSlice[0])
		startPos, _ := strconv.Atoi(binSlice[1])
		newStartPos := startPos - upSize
		if newStartPos < 1 {
			newStartPos = 1
		}
		endPos, _ := strconv.Atoi(binSlice[2])
		newEndPos := endPos + downSize
		binPosSlice := []int{chrPos, newStartPos, newEndPos}
		binFlankSlice := []string{}
		for _, valBinPosSlice := range binPosSlice {
			binFlankSlice = append(binFlankSlice, strconv.Itoa(valBinPosSlice))
		}
		binFlankSlice = append(binFlankSlice, binSlice[3])
		binFlankStr := strings.Join(binFlankSlice, "\t")
		newGeneFlank = append(newGeneFlank, binFlankStr)
	}
	return newGeneFlank
}

//worker function for making worker pools
func worker(vcfLine map[string]string, jobs <-chan string, results chan<- string) {
	for j := range jobs {
		binSlice := strings.Fields(j)
		chrPos, _ := strconv.Atoi(binSlice[0])
		startPos, _ := strconv.Atoi(binSlice[1])
		endPos, _ := strconv.Atoi(binSlice[2])
		binPosSlice := []int{chrPos, startPos, endPos}
		geneAveIndex := []string{binSlice[3]}
		geneAveIndex = append(geneAveIndex, calculateAveIndex(binPosSlice, vcfLine))
		geneAveIndexStr := strings.Join(geneAveIndex, "\t")

		results <- geneAveIndexStr
	}
}

// maxParallelism function will return number of CPUs
func maxParallelism(cpus int) int {
	maxProcs := runtime.GOMAXPROCS(cpus)
	numCPU := runtime.NumCPU()
	if maxProcs < numCPU {
		return maxProcs
	}
	return numCPU
}

// main function
func main() {
	flag.Parse()

	fmt.Println("[", time.Now(), "] ", "Program start ...")
	numThreads := maxParallelism(*cpus)
	fmt.Println("Total available CPU number is: ", runtime.NumCPU())
	fmt.Println("Working CPU number is: ", numThreads)

	vcfLines, err := getSNVlong(*vcfFile)
	if err != nil {
		log.Fatalf("read input vcf file: %s", err)
	}

	bedLines, err := getLines(*bedFile)
	if err != nil {
		log.Fatalf("read input bed file: %s", err)
	}

	binSlice := getGeneBin(bedLines, *upSize, *downSize)
	fmt.Println("The total number of genes is: ", len(binSlice))

	//start worker
	jobs := make(chan string, len(binSlice))
	results := make(chan string, len(binSlice))
	for w := 1; w <= numThreads; w++ {
		go worker(vcfLines, jobs, results)
	}

	for _, valBinSlice := range binSlice {
		jobs <- valBinSlice
	}
	close(jobs)

	binLines := []string{}
	for a := 1; a <= len(binSlice); a++ {
		binLines = append(binLines, <-results)
	}
	//fmt.Println("The test is: ", len(binLines))

	mapBinLines := map[string]string{}
	for _, valBinLine := range binLines {
		valBinLineSlice := strings.Fields(valBinLine)
		keyMapBinLines := []string{valBinLineSlice[1], valBinLineSlice[2], valBinLineSlice[3], valBinLineSlice[0]}
		valMapBinlines := append(keyMapBinLines, valBinLineSlice[4:]...)
		mapBinLines[strings.Join(keyMapBinLines, "\t")] = strings.Join(valMapBinlines, "\t")
	}

	binHeader := []string{
		"#CHR", "START", "END", "mid_pos", "num_SNVs",
		"AveIdx_v3", "AveIdx_ce", "AveIdx_wt", "AveIdx_mt",
		"AveIdx_parent", "AveIdx_F2", "AveDp_wt", "AveDp_mt",
		"Ave_p90L", "Ave_p90H", "Ave_p95L", "Ave_p95H", "Ave_p99L", "Ave_p99H",
	}

	newBinLines := []string{}
	newBinLines = append(newBinLines, strings.Join(binHeader, "\t"))
	for _, valBinSlice := range binSlice {
		newBinLines = append(newBinLines, mapBinLines[valBinSlice])
	}

	//write skip snv file
	if err := writeLines(newBinLines, *outFile); err != nil {
		log.Fatalf("write gene intervals: $s", err)
	}

	fmt.Println("[", time.Now(), "] ", "Program end ...")
}
