package main

import (
	"fmt"
	"os"

	// Directories in the root of the repo can be imported
	// as long as we pretend that they sit relative to the
	// url birc.au.dk/gsa, like this for the example 'shared':
	gsa "birc.au.dk/gsa/helpers"
)

func main() {
	if len(os.Args) != 3 {
		fmt.Fprintf(os.Stderr, "Usage: genome-file reads-file\n")
		os.Exit(1)
	}
	genome := os.Args[1]
	reads := os.Args[2]

	parsedGenomes := gsa.GeneralParser(genome, gsa.Fasta)

	parsedReads := gsa.GeneralParser(reads, gsa.Fastq)

	for _, read := range parsedReads {
		for _, gen := range parsedGenomes {
			matches := lin(gen.Rec, read.Rec)
			for _, match := range matches {
				gsa.Sam(read.Name, gen.Name, match, read.Rec)
			}
		}
	}

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
