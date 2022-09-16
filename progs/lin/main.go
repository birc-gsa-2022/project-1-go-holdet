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
outer_loop:
	for i < len(x)-len(p)+1 {
		for j, char := range []byte(p) {
			if char != x[i+j] {
				if j == 0 {
					i += 1
					continue outer_loop
				}
				i += b[j-1] + 1
				continue outer_loop
			}
		}
		matches = append(matches, i)
		if b[len(b)-1] == 0 {
			i += len(b)
		} else {
			i += b[len(b)-1]
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
