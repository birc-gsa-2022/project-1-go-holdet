package test

import (
	"fmt"
	"math/rand"
	"strconv"
	"strings"
	"testing"

	gsa "birc.au.dk/gsa/helpers"
)

func TestOutput(t *testing.T) {
	genome, reads := BuildSomeFastaAndFastq(1000, 100, Repetitive, 3)
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

func runNaive(genomeString string, readsString string) []string {

	parsedGenomes := gsa.GeneralParserStub(genomeString, gsa.Fasta)

	parsedReads := gsa.GeneralParserStub(readsString, gsa.Fastq)

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

	parsedGenomes := gsa.GeneralParserStub(genomeString, gsa.Fasta)

	parsedReads := gsa.GeneralParserStub(readsString, gsa.Fastq)
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
	English    Alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
	DNA        Alphabet = "ACTG"
	Repetitive Alphabet = "ab"
	A          Alphabet = "a"
)

func randString(number int, alphabet Alphabet) string {
	var letters = alphabet
	b := make([]byte, number)
	for i := range b {
		b[i] = letters[rand.Intn(len(letters))]
	}
	return string(b)
}

func BuildSomeFastaAndFastq(longest int, chr_amount int, alphabet Alphabet, seed int64) (string, string) {

	rand.Seed(seed)

	var sb strings.Builder
	var sc strings.Builder

	for i := 0; i < chr_amount; i++ {
		sb.WriteString("> chr" + strconv.Itoa(i) + "\n")

		b_str := randString(longest, alphabet)
		sb.WriteString(b_str + "\n")

		sc.WriteString("@read" + strconv.Itoa(i) + "\n")

		idx := rand.Intn(len(b_str))
		c_str := ""
		if len(b_str) < idx+3 {
			c_str = b_str[idx:]
		} else {
			c_str = b_str[idx : idx+3]
		}
		sc.WriteString(c_str + "\n")
	}

	return sb.String(), sc.String()
}
