package gsa

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

func FastaSubseqs(file string) {
	fastaFile, err := os.Open(file)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s", err.Error())
		os.Exit(1)
	}
	defer fastaFile.Close()
	//########################################################
	sequences := make([]string, 0)
	fileScanner := bufio.NewScanner(fastaFile)

	//scan file line by line
	for fileScanner.Scan() {

		line := fileScanner.Text()

		if len(line) == 0 {
			continue
		}

		//handle 'name of sequence' cases
		if line[0:1] == ">" {
			sequences = append(sequences, "")
			sequences[len(sequences)-1] = strings.TrimSpace(line[1:])

			//handle 'sequence' cases
		} else {
			sequences[len(sequences)-1] = sequences[len(sequences)-1] + strings.TrimSpace(line)
		}
	}
	fmt.Println(sequences)

	// #######################################

	var coordFile = os.Stdin
	if len(os.Args) == 3 && os.Args[2] != "-" {
		coordFile, err = os.Open(os.Args[2])
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s", err.Error())
			os.Exit(1)
		}
		defer coordFile.Close()
	}

	//######################

	output := ""
	fileScannerCoord := bufio.NewScanner(coordFile)

	for fileScannerCoord.Scan() {
		lineCoord := fileScannerCoord.Text()
		if len(lineCoord) == 0 {
			continue
		}
		split_str := strings.Fields(lineCoord)
		if len(split_str) != 3 {
			fmt.Println("Error, length of coordinates")
		}
		for i, v := range split_str {
			split_str[i] = strings.TrimSpace(v)
		}

		for _, v := range sequences {

			if strings.HasPrefix(v, split_str[0]) {

				basepair := strings.TrimSpace(strings.TrimPrefix(v, split_str[0]))
				low, er1 := strconv.Atoi(split_str[1])
				high, er2 := strconv.Atoi(split_str[2])
				if er1 != nil || er2 != nil {
					fmt.Printf("Error")

				}
				output = output + basepair[low-1:high-1] + "\n"
				break
			}

		}
	}
	fmt.Print(output)

	//#########################################
}

// Borderarray computes the border array over the string x. The border
// array ba will at index i have the length of the longest proper border
// of the string x[:i+1], i.e. the longest non-empty string that is both
// a prefix and a suffix of x[:i+1].
func Borderarray(x string) []int {
	ba := make([]int, len(x))
	fmt.Println("ASDF")
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

// StrictBorderarray computes the strict border array over the string x.
// This is almost the same as the border array, but ba[i] will be the
// longest proper border of the string x[:i+1] such that x[ba[i]] != x[i].
func StrictBorderarray(x string) []int {
	ba := make([]int, len(x))

	b := 0
	i := 1
	for i < len(x) {
		if x[b] == x[i] {
			b++
			if i == len(x)-1 {
				ba[i] = b
			}
			i++
			continue
		}
		if ba[i-1] == 0 {
			ba[i-1] = b
		}
		if b > 0 {
			b = ba[b-1]
			continue
		}
		i++
	}
	return ba

	/* **simpler but 2*n solution***
	ba = Borderarray(x)
	for i := 1; i < len(ba); i++ {
		if ba[i-1] == ba[i]-1 {
			ba[i-1] = 0
		}
	}
	return ba
	*/
}

func Align(p, q, edits string) (pRow, qRow string) {

	var sb_p, sb_q strings.Builder
	var idx_p, idx_q int

	mid := []rune("MID")
	// Align p and q based on edits
	for _, v := range edits {
		if v == mid[0] {
			sb_p.WriteString(string(p[idx_p]))
			sb_q.WriteString(string(q[idx_q]))
			idx_p++
			idx_q++
			continue
		}
		if v == mid[1] {
			sb_q.WriteByte(q[idx_q])
			sb_p.WriteString("-")
			idx_q++
			continue
		}
		if v == mid[2] {
			sb_p.WriteByte(p[idx_p])
			sb_q.WriteString("-")
			idx_p++

		}
	}

	pRow = sb_p.String()
	qRow = sb_q.String()

	return pRow, qRow
}

// Align two sequences from a sequence of edits.
//
//  Args:
//      p: The first sequence to align
//      x: The second sequence to align; we only align locally
//      i: Start position of the alignment in x
//      edits: The list of edits to apply, given as a string
//
//  Returns:
//      The two rows in the pairwise alignment
func LocalAlign(p, x string, i int, edits string) (pRow, xRow string) {
	pRow, xRow = "", ""
	// Align p and q based on edits
	pRow, xRow = Align(p, x[i:], edits)
	return pRow, xRow
}

func CigarToEdits(cigar string) (edits string) {
	var sb strings.Builder
	i := 0
	for i < len(cigar) {
		numb, _ := strconv.Atoi(string(cigar[i]))
		op := string(cigar[i+1])
		i += 2

		sb.WriteString(strings.Repeat(op, numb))
	}
	edits = sb.String()
	return edits
}

// Encode a sequence of edits as a CIGAR.
//
//  Args:
//      edits: A sequence of edit operations
//
//  Returns:
//      The CIGAR encoding of edits.
func EditsToCigar(edits string) (cigar string) {
	var sb strings.Builder
	var cur byte
	app := 0
	for i, v := range []byte(edits) {
		if i == 0 {
			cur = v
			app = 1
			continue
		}
		if cur == v {
			app++
		}
		if cur != v {
			sb.WriteString(fmt.Sprint(app))
			sb.WriteByte(cur)
			cur = v
			app = 1
		}
		if i == len(edits)-1 {
			sb.WriteString(fmt.Sprint(app))
			sb.WriteByte(cur)
		}
	}
	cigar = sb.String()
	return cigar
}

func GetEdits(p, q string) (gapFreeP, gapFreeQ, edits string) {
	var sb_p, sb_q, sb_edits strings.Builder

	for i := 0; i < len(p); i++ {

		if p[i] == []byte("-")[0] {
			sb_edits.WriteString("I")
			sb_q.WriteByte(q[i])
			continue
		}
		if q[i] == []byte("-")[0] {
			sb_edits.WriteString("D")
			sb_p.WriteByte(p[i])
			continue
		}
		sb_edits.WriteString("M")
		sb_p.WriteByte(p[i])
		sb_q.WriteByte(q[i])

	}
	edits = sb_edits.String()
	gapFreeP = sb_p.String()
	gapFreeQ = sb_q.String()

	return gapFreeP, gapFreeQ, edits
}

//  Get the distance between p and the string that starts at x[i:]
//  using the edits.
//
//  Args:
//      p: The read string we have mapped against x
//      x: The longer string we have mapped against
//      i: The location where we have an approximative match
//      edits: The list of edits to apply, given as a string
//
//  Returns:
//      The distance from p to x[i:?] described by edits
func EditDist(p, x string, i int, edits string) int {
	pRow, xRow := LocalAlign(p, x, i, edits)
	result := 0
	for i := 0; i < len(pRow); i++ {
		if pRow[i] != xRow[i] {
			result++
		}
	}

	return result
}
