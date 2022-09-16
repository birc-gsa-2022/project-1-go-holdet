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
			matches := Naive(gen.Rec, read.Rec)
			for _, match := range matches {
				gsa.Sam(read.Name, gen.Name, match, read.Rec)
			}
		}
	}
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
