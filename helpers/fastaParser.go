package gsa

import (
	"bufio"
	"fmt"
	"os"
	"strings"
)

type Recs struct {
	Name string
	Rec  string
}

type Format string

const (
	Fasta Format = ">"
	Fastq Format = "@"
)

func GeneralParser(file string, format Format) []Recs {
	f, err := os.Open(file)

	if err != nil {
		fmt.Fprintf(os.Stderr, "%s", err.Error())
		os.Exit(1)
	}
	defer f.Close()

	//###########################################

	output := ""
	fileScanner := bufio.NewScanner(f)
	activeRec := new(Recs)

	recs := make([]Recs, 0)
	//scan file line by line
	for fileScanner.Scan() {
		line := fileScanner.Text()

		if len(line) == 0 {
			continue
		}

		//handle 'name of sequence' cases
		if line[0:1] == string(format) {
			if len(output) != 0 {
				activeRec.Rec = output
				recs = append(recs, *activeRec)
			}
			output = ""
			activeRec = new(Recs)
			activeRec.Name = strings.TrimSpace(line[1:])
			//handle 'sequence' cases
		} else {
			output = output + line
		}
	}
	activeRec.Rec = output
	recs = append(recs, *activeRec)

	return recs
}

func GeneralParserStub(file string, format Format) []Recs {
	output := ""
	fileScanner := bufio.NewScanner(strings.NewReader(file))
	activeRec := new(Recs)

	recs := make([]Recs, 0)
	//scan file line by line
	for fileScanner.Scan() {
		line := fileScanner.Text()

		if len(line) == 0 {
			continue
		}

		//handle 'name of sequence' cases
		if line[0:1] == string(format) {
			if len(output) != 0 {
				activeRec.Rec = output
				recs = append(recs, *activeRec)
			}
			output = ""
			activeRec = new(Recs)
			activeRec.Name = strings.TrimSpace(line[1:])
			//handle 'sequence' cases
		} else {
			output = output + line
		}
	}
	activeRec.Rec = output
	recs = append(recs, *activeRec)

	return recs
}
