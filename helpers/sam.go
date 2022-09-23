package gsa

import (
	"fmt"
	"strconv"
)

func Sam(readName string, chrom string, pos int, readString string) {
	output := readName + "	" + chrom + "	" + fmt.Sprint(pos+1) + "	" + strconv.Itoa(len(readString)) + "M" + "	" + readString + "\n"

	fmt.Print(output)
}

func SamStub(readName string, chrom string, pos int, readString string) string {
	output := readName + "	" + chrom + "	" + fmt.Sprint(pos+1) + "	" + strconv.Itoa(len(readString)) + "M" + "	" + readString + "\n"

	return output
}
