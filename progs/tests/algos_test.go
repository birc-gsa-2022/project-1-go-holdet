package test

import (
	"encoding/csv"
	"fmt"
	"log"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"testing"
	"time"

	gsa "birc.au.dk/gsa/helpers"
)

func TestOutput(t *testing.T) {
	genome, reads := BuildSomeFastaAndFastq(100, 20, 1500, AB, 3)
	resNaive := runNaive(genome, reads)
	resLin := runLin(genome, reads)
	fmt.Println(len(resNaive))
	fmt.Println(len(resLin))

	for pos, s1 := range resNaive {
		/*
			fmt.Println(s1)
			fmt.Println(resLin[pos])
			fmt.Println("")
		*/
		if s1 != resLin[pos] {
			t.Error("error at ", pos, "naive: "+s1+"   lin: "+resLin[pos])
		}
	}
}

func TestVaryingAlphabets(t *testing.T) {
	time_Naive := 0
	time_Lin := 0

	Alphabets := []Alphabet{
		English, DNA, AB, A}

	for _, v := range Alphabets {
		genome, reads := BuildSomeFastaAndFastq(300, 200, 500, v, 3)

		parsedGenomes := gsa.GeneralParserStub(genome, gsa.Fasta, len(genome)+1)

		parsedReads := gsa.GeneralParserStub(reads, gsa.Fastq, len(reads)+1)

		for _, read := range parsedReads {
			for _, gen := range parsedGenomes {
				time_Start := time.Now()

				Naive(gen.Rec, read.Rec)
				time_Naive += int(time.Since(time_Start))

				time_Start = time.Now()
				lin(gen.Rec, read.Rec)
				time_Lin += int(time.Since(time_Start))

			}
		}

	}
	fmt.Println(time_Naive / 1000)
	fmt.Println(time_Lin / 1000)

}

func TestMakeDataFixN(t *testing.T) {
	csvFile, err := os.Create("fixed_n_data.csv")
	if err != nil {
		log.Fatalf("failed creating file: %s", err)
	}
	csvwriter := csv.NewWriter(csvFile)
	_ = csvwriter.Write([]string{"x_size", "p_size", "naive", "kmp"})

	num_of_n := 0
	time_Naive := 0
	time_Lin := 0

	for i := 1; i < 15; i++ {

		num_of_n *= 2
		num_of_m := 200
		genome, reads := BuildARepetitiveFastaAndFastq(20, num_of_n, 78)
		parsedGenomes := gsa.GeneralParserStub(genome, gsa.Fasta, num_of_n*num_of_m+1)
		parsedReads := gsa.GeneralParserStub(reads, gsa.Fastq, num_of_n*num_of_m+1)

		for i := 0; i < 10; i++ {
			for _, read := range parsedReads {
				for _, gen := range parsedGenomes {
					time_StartN := time.Now()

					//resultN := Naive(gen.Rec, read.Rec)
					time_endN := int(time.Since(time_StartN))
					time_Naive += time_endN

					time_StartL := time.Now()
					lin(gen.Rec, read.Rec)
					time_endL := int(time.Since(time_StartL))
					time_Lin += time_endL

					//fmt.Println(resultN)
					//for i, v := range resultN {

					//	if v != resultL[i] {
					//		t.Error("not same result")
					//	}
					//}

				}
			}
			fmt.Println("NAIVE", int((time_Naive)))
			fmt.Println("LINEA", int(time_Lin))
			_ = csvwriter.Write([]string{strconv.Itoa(num_of_n), strconv.Itoa(num_of_m), strconv.Itoa(time_Naive), strconv.Itoa(time_Lin)})
			time_Naive, time_Lin = 0, 0
		}

		csvwriter.Flush()

	}

}

func TestDoesMakeDataWork(t *testing.T) {
	size := 100
	fasta, fastq := BuildSomeFastaAndFastq(size, size, 1, DNA, 0)
	x := gsa.GeneralParserStub(fasta, gsa.Fasta, size+1)
	p := gsa.GeneralParserStub(fastq, gsa.Fastq, size+1)
	if x[0].Rec != p[0].Rec {
		t.Error("not identical genomes")
	}

	expected_len := 266000

	fastaa, _ := BuildSomeFastaAndFastq(expected_len, 100, 1, DNA, 0)
	xx := gsa.GeneralParserStub(fastaa, gsa.Fasta, expected_len+1)
	//pp := gsa.GeneralParserStub(fastqq, gsa.Fastq)
	genomexx := xx[0].Rec
	if len(genomexx) != expected_len {
		fmt.Println(len(genomexx), expected_len)
		t.Error("sus")
	}
}

func runNaive(genomeString string, readsString string) []string {

	parsedGenomes := gsa.GeneralParserStub(genomeString, gsa.Fasta, len(genomeString))

	parsedReads := gsa.GeneralParserStub(readsString, gsa.Fastq, len(genomeString))

	sams := make([]string, 0)
	for _, read := range parsedReads {
		for _, gen := range parsedGenomes {
			matches := Naive(gen.Rec, read.Rec)
			for _, match := range matches {
				sams = append(sams, gsa.SamStub(read.Name, gen.Name, match, read.Rec))
			}
		}
	}
	return sams
}

func Naive(x string, p string) (matches []int) {
outer_loop:
	for i := 0; i < len(x)-len(p)+1; i++ {
		for j, char := range []byte(p) {
			if char != x[i+j] {
				continue outer_loop
			}
		}
		matches = append(matches, i)
	}

	return matches
}

func runLin(genomeString string, readsString string) []string {

	parsedGenomes := gsa.GeneralParserStub(genomeString, gsa.Fasta, len(genomeString))

	parsedReads := gsa.GeneralParserStub(readsString, gsa.Fastq, len(genomeString))
	sams := make([]string, 0)

	for _, read := range parsedReads {
		for _, gen := range parsedGenomes {
			matches := lin(gen.Rec, read.Rec)
			for _, match := range matches {
				sams = append(sams, gsa.SamStub(read.Name, gen.Name, match, read.Rec))
			}
		}
	}
	return sams

}

func lin(x string, p string) (matches []int) {
	b := Borderarray(p)
	i := 0
	j := 0

	for len(x)-i >= (len(p) - j) {
		if x[i] == p[j] {
			i += 1
			j += 1
		}

		if j == len(p) {
			matches = append(matches, i-j)
			j = b[j-1]
		} else if i < len(x) && x[i] != p[j] {
			if j != 0 {
				j = b[j-1]
			} else {
				i += 1
			}
		}
	}
	return matches
}

func Borderarray(x string) []int {
	ba := make([]int, len(x))

	b := 0
	i := 1
	for i < len(x) {
		if x[b] == x[i] {
			b++
			ba[i] = b
			i++
			continue
		}
		if b > 0 {
			b = ba[b-1]
			continue
		}
		i++
	}

	return ba
}

type Alphabet string

const (
	English Alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
	DNA     Alphabet = "ACTG"
	AB      Alphabet = "ab"
	A       Alphabet = "a"
)

func randString(number int, alphabet Alphabet) string {
	var letters = alphabet
	b := make([]byte, number)
	for i := range b {
		b[i] = letters[rand.Intn(len(letters))]
	}
	return string(b)
}

func BuildSomeFastaAndFastq(len_Fasta int, len_Fastq int, amount int, alphabet Alphabet, seed int64) (string, string) {

	rand.Seed(seed)

	var sb strings.Builder
	var sc strings.Builder

	for i := 0; i < amount; i++ {
		sb.WriteString("> chr" + strconv.Itoa(i) + "\n")

		b_str := randString(len_Fasta, alphabet)
		sb.WriteString(b_str + "\n")

		sc.WriteString("@read" + strconv.Itoa(i) + "\n")

		dif := len_Fasta - len_Fastq
		idx := 0
		//allows for fasta and fastq to be same len. rand.intn panics for n=0 input.
		if dif > 0 {
			idx = rand.Intn(len_Fasta - len_Fastq)

		}
		c_str := b_str[idx:(idx + len_Fastq)]

		sc.WriteString(c_str + "\n")
	}

	return sb.String(), sc.String()
}

func BuildARepetitiveFastaAndFastq(repetitions int, len_Fastq int, seed int64) (string, string) {
	rand.Seed(seed)

	var sb strings.Builder
	var sc strings.Builder

	single_rep := strings.Repeat("a", len_Fastq-1)

	sb.WriteString("> chr" + strconv.Itoa(0) + "\n")
	sb.WriteString(strings.Repeat(single_rep+"b", repetitions))

	sb.WriteString("@ read" + strconv.Itoa(0) + "\n")
	sc.WriteString(single_rep + "c")

	return sb.String(), sc.String()

}

func BuildARepetitiveFastaAndFastqNotQuadratic(len_Fasta int, len_Fastq int, seed int64) (string, string) {
	//crashes if len fastq cant be divided by len fasta
	rand.Seed(seed)

	var sb strings.Builder
	var sc strings.Builder

	single_rep := strings.Repeat("a", len_Fastq-1)

	sb.WriteString("> chr" + strconv.Itoa(0) + "\n")
	sb.WriteString(strings.Repeat(single_rep+"b", len_Fasta/len_Fastq))

	sb.WriteString("@ read" + strconv.Itoa(0) + "\n")
	sc.WriteString(single_rep + "a")

	return sb.String(), sc.String()

}
