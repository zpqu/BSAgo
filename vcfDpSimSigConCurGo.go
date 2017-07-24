package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"regexp"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/leesper/go_rng"
)

// declare command arguments
var (
	passFile  = flag.String("in", "", "Input pass vcf file")
	dpFile    = flag.String("dp", "", "Output depth simulation file")
	outFile   = flag.String("out", "", "Output file with index and sig merged")
	popStruct = flag.String("p", "", "Population struction: RIL or F2")
	numIndvl  = flag.Int("n", 0, "Number of individuals in each bulk")
	rep       = flag.Int("r", 0, "Number of replication in simulation")
	filterVal = flag.Float64("f", 0.3, "Filter value")
	cpus      = flag.Int("c", 1, "Number of workding CPUs")
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

// getDP function will retrieve depth information from VCF files and store
// them into a string, this function is used to keep the order of all DPs
func getDP(vcfLine string) string {
	lineSlice := strings.Fields(vcfLine)

	wtDp := "0"
	if lineSlice[3] == "." {
		wtDp = "0"
	} else {
		wtDpSlice := strings.Split(lineSlice[3], ":")
		wtDpAF := strings.Split(wtDpSlice[2], ",")
		wtRef, _ := strconv.Atoi(wtDpAF[0])
		wtAlt, _ := strconv.Atoi(wtDpAF[1])
		wtDp = strconv.Itoa(wtRef + wtAlt)
	}

	mtDp := "0"
	if lineSlice[4] == "." {
		mtDp = "0"
	} else {
		mtDpSlice := strings.Split(lineSlice[4], ":")
		mtDpAF := strings.Split(mtDpSlice[2], ",")
		mtRef, _ := strconv.Atoi(mtDpAF[0])
		mtAlt, _ := strconv.Atoi(mtDpAF[1])
		mtDp = strconv.Itoa(mtRef + mtAlt)
	}
	dpMapKeySlice := []string{}
	dpMapKey := ""
	if wtDp == "0" && mtDp == "0" {
	} else {
		if wtDp == "0" {
			dpMapKeySlice = append(dpMapKeySlice, mtDp, mtDp)
		} else if mtDp == "0" {
			dpMapKeySlice = append(dpMapKeySlice, wtDp, wtDp)
		} else {
			dpMapKeySlice = append(dpMapKeySlice, wtDp, mtDp)
		}
		dpMapKey = strings.Join(dpMapKeySlice, "_")
	}
	return dpMapKey
}

// getDPmap function will retrieve depth information for two F2s from VCF files and
// store them in a map
func getDPmap(vcfLines []string) map[string]string {
	dpMap := map[string]string{}
	for _, valLine := range vcfLines {
		lineSlice := strings.Fields(valLine)

		wtDp := "0"
		if lineSlice[3] == "." {
			wtDp = "0"
		} else {
			wtDpSlice := strings.Split(lineSlice[3], ":")
			wtDpAF := strings.Split(wtDpSlice[2], ",")
			wtRef, _ := strconv.Atoi(wtDpAF[0])
			wtAlt, _ := strconv.Atoi(wtDpAF[1])
			wtDp = strconv.Itoa(wtRef + wtAlt)
		}

		mtDp := "0"
		if lineSlice[4] == "." {
			mtDp = "0"
		} else {
			mtDpSlice := strings.Split(lineSlice[4], ":")
			mtDpAF := strings.Split(mtDpSlice[2], ",")
			mtRef, _ := strconv.Atoi(mtDpAF[0])
			mtAlt, _ := strconv.Atoi(mtDpAF[1])
			mtDp = strconv.Itoa(mtRef + mtAlt)
		}
		dpMapKeySlice := []string{}
		if wtDp == "0" && mtDp == "0" {
			fmt.Println("Skip depth: ", valLine)
		} else {
			if wtDp == "0" {
				dpMapKeySlice = append(dpMapKeySlice, mtDp, mtDp)
			} else if mtDp == "0" {
				dpMapKeySlice = append(dpMapKeySlice, wtDp, wtDp)
			} else {
				dpMapKeySlice = append(dpMapKeySlice, wtDp, mtDp)
			}
			dpMapKey := strings.Join(dpMapKeySlice, "_")
			dpMapValSlice := []string{wtDp, mtDp}
			dpMapVal := strings.Join(dpMapValSlice, "_")
			dpMap[dpMapKey] = dpMapVal
		}
	}
	fmt.Println("Total number of depth is: ", len(dpMap))
	return dpMap
}

// genotype function will randomly get genotype given population struction
func genotype(popStrut string) float64 {
	count := 0.0
	if popStrut == "RIL" {
		randge := rng.NewUniformGenerator(time.Now().UnixNano())
		frq := randge.Float64()

		if frq <= 0.5 {
			count = 1.0
		} else {
			count = 0.0
		}
	} else {
		for i := 1; i <= 2; i++ {
			randge := rng.NewUniformGenerator(time.Now().UnixNano())
			frq := randge.Float64()
			if frq <= 0.5 {
				count += 0.5
			}
		}
	}
	return count
}

// calIndvlGeno function will calculate individual genotype using randomly
// generated genotype
func calIndvlGeno(numIndvl int, popStrut string) float64 {
	genoTotal := 0.0
	for j := 1; j <= numIndvl; j++ {
		genoTotal += genotype(popStrut)
	}
	indvlGeno := genoTotal / float64(int64(numIndvl))
	return indvlGeno
}

// calIndvlIndex function will calculate indeividual index using binomial
// distribution with given depth and expected genotype ratio (default 0.3)
func calIndvlIndex(dp int, ratioGeno float64) float64 {
	binomalge := rng.NewBinomialGenerator(time.Now().UnixNano())
	indvlIndex := binomalge.Binomial(int64(dp), ratioGeno)
	indvlIndexAve := float64(indvlIndex) / float64(dp)
	return indvlIndexAve
}

// simIndex function will do QTL simulation with given times of replication
func simIndex(numIndvl int, dpSlice []int, rep int, filterVal float64, popStrut string) []float64 {
	p90L := 0.0
	p90H := 0.0
	p95L := 0.0
	p95H := 0.0
	p99L := 0.0
	p99H := 0.0
	delIndvlIndexSlice := []float64{}
	for k := 1; k <= rep; k++ {
		wtRatioGeno := calIndvlGeno(numIndvl, popStrut)
		wtIndvlIndex := calIndvlIndex(dpSlice[0], wtRatioGeno)
		mtRatioGeno := calIndvlGeno(numIndvl, popStrut)
		mtIndvlIndex := calIndvlIndex(dpSlice[1], mtRatioGeno)

		if wtIndvlIndex >= filterVal || mtIndvlIndex >= filterVal {
			delIndvlIndex := wtIndvlIndex - mtIndvlIndex
			delIndvlIndexSlice = append(delIndvlIndexSlice, delIndvlIndex)
		}
	}
	sort.Float64s(delIndvlIndexSlice)
	delIndvlIndexLen := float64(len(delIndvlIndexSlice))

	//delta index probability 0.1
	if math.Floor(0.05*delIndvlIndexLen) > 0.0 {
		p90L = delIndvlIndexSlice[int(math.Floor(0.05*delIndvlIndexLen))]
	} else {
		p90L = delIndvlIndexSlice[0]
	}

	if math.Ceil(0.95*delIndvlIndexLen) < delIndvlIndexLen {
		p90H = delIndvlIndexSlice[int(math.Ceil(0.95*delIndvlIndexLen))]
	} else {
		p90H = delIndvlIndexSlice[int(len(delIndvlIndexSlice))-1]
	}

	//delta index probability 0.05
	if math.Floor(0.025*delIndvlIndexLen) > 0.0 {
		p95L = delIndvlIndexSlice[int(math.Floor(0.025*delIndvlIndexLen))]
	} else {
		p95L = delIndvlIndexSlice[0]
	}

	if math.Ceil(0.975*delIndvlIndexLen) < delIndvlIndexLen {
		p95H = delIndvlIndexSlice[int(math.Ceil(0.975*delIndvlIndexLen))]
	} else {
		p95H = delIndvlIndexSlice[int(len(delIndvlIndexSlice))-1]
	}

	//delta index probability 0.01
	if math.Floor(0.005*delIndvlIndexLen) > 0.0 {
		p99L = delIndvlIndexSlice[int(math.Floor(0.005*delIndvlIndexLen))]
	} else {
		p99L = delIndvlIndexSlice[0]
	}

	if math.Ceil(0.995*delIndvlIndexLen) < delIndvlIndexLen {
		p99H = delIndvlIndexSlice[int(math.Ceil(0.995*delIndvlIndexLen))]
	} else {
		p99H = delIndvlIndexSlice[int(len(delIndvlIndexSlice))-1]
	}
	sigSlice := []float64{p90L, p90H, p95L, p95H, p99L, p99H}
	return sigSlice
}

// worker function makes working pools to do QTL simulation and return 4 confidence
// intervals: 95% low, 95% high, 99% low, 99% high
func worker(numIndvl int, rep int, filterVal float64, popStrut string, jobs <-chan string, results chan<- string) {
	for j := range jobs {
		keyDpSlice := strings.Split(j, "_")
		wtDp, _ := strconv.Atoi(keyDpSlice[0])
		mtDp, _ := strconv.Atoi(keyDpSlice[1])

		dpIntSlice := []int{wtDp, mtDp}
		dpSimIndex := simIndex(numIndvl, dpIntSlice, rep, filterVal, popStrut)
		vcfDpIndexSlice := []string{keyDpSlice[0], keyDpSlice[1]}
		for _, valDpSim := range dpSimIndex {
			vcfDpIndexSlice = append(vcfDpIndexSlice, strconv.FormatFloat(valDpSim, 'f', 2, 64))
		}
		results <- strings.Join(vcfDpIndexSlice, "\t")
	}
}

// maxParallelism function will return working number of CPUs, if requested number
// of CPUs is greater than available CPUs, use available CPUs (machine CPUs)
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

	vcfLines, err := getSNVlong(*passFile)
	if err != nil {
		log.Fatalf("read input pass vcf file: %s", err)
	}

	vcfDpMap := getDPmap(vcfLines)

	indexHeader := []string{
		"#ID", "DP_wt", "DP_mt", "p90L", "p90H",
		"p95L", "p95H", "p99L", "p99H",
	}

	numThreads := maxParallelism(*cpus)
	fmt.Println("Total available CPU number is: ", runtime.NumCPU())
	fmt.Println("Working CPU number is: ", numThreads)
	jobs := make(chan string, len(vcfDpMap))
	results := make(chan string, len(vcfDpMap))

	for w := 1; w <= numThreads; w++ {
		go worker(*numIndvl, *rep, *filterVal, *popStruct, jobs, results)
	}

	for keyDpMap, _ := range vcfDpMap {
		jobs <- keyDpMap
	}
	close(jobs)

	dpLines := []string{}
	dpLines = append(dpLines, strings.Join(indexHeader, "\t"))
	for a := 1; a <= len(vcfDpMap); a++ {
		dpLines = append(dpLines, <-results)
	}

	//write skip snv file
	if err := writeSNVlong(dpLines, *dpFile); err != nil {
		log.Fatalf("write out dp: $s", err)
	}

	//make new vcfDpMap
	dpLineMap := map[string]string{}
	for _, valDpLine := range dpLines {
		dpLineSlice := strings.Fields(valDpLine)
		dpLineMap[strings.Join(dpLineSlice[:2], "_")] = valDpLine
	}

	//merge index with significance
	outHeader := []string{
		"#ID", "geno_v3", "geno_WT", "geno_WT", "gene_MT",
		"idx_v3", "idx_CE", "idx_WT", "idx_MT", "deltaIdx_parent",
		"deltaIdx_F2", "dp_WT", "dp_MT", "p90L", "p90H",
		"p95L", "p95H", "p99L", "p99H",
	}

	mergedVcfLines := []string{}
	mergedVcfLines = append(mergedVcfLines, strings.Join(outHeader, "\t"))
	for _, valVcfLine := range vcfLines {
		vcfDpKey := getDP(valVcfLine)
		if valDp, ok := dpLineMap[vcfDpKey]; ok {
			newVcfSlice := []string{valVcfLine, valDp}
			mergedVcfLines = append(mergedVcfLines, strings.Join(newVcfSlice, "\t"))
		} else {
			fmt.Println("Skiped vcf: ", valVcfLine)
		}
	}

	//write merged file
	if err := writeSNVlong(mergedVcfLines, *outFile); err != nil {
		log.Fatalf("write out merged index and sig: $s", err)
	}

	fmt.Println("[", time.Now(), "] ", "Program end ...")

}
