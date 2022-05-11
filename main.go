package lanthanideTraj

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
)

var shellDist = 2.8

func main() {
	atoms := readFile("tail.arc")
	writeFile("tail.xyz", atoms)
}

func writeFile(outPath string, atoms []atom) {
	thisFile, err := os.Create(outPath)
	if err != nil {
		fmt.Println("Failed to create new fragment file: " + outPath)
		log.Fatal(err)
	}
	_, _ = thisFile.WriteString(strconv.Itoa(len(atoms)) + "\n\n")
	for _, thisAtom := range atoms {
		outString := strconv.Itoa(thisAtom.newID) + " " + thisAtom.element + " " + fmt.Sprintf("%.6f", thisAtom.x) + " " +
			fmt.Sprintf("%.6f", thisAtom.y) + " " + fmt.Sprintf("%.6f", thisAtom.z) + " " + thisAtom.atomType + " "
		for _, bondID := range thisAtom.bonds {
			outString += strconv.Itoa(bondID) + " "
		}
		_, _ = thisFile.WriteString(outString + "\n")
	}
}

func readFile(filePath string) []atom {
	// open file
	file, err := os.Open(filePath)
	if err != nil {
		fmt.Println("Failed to open molecule file: " + filePath)
		log.Fatal(err)
	}

	scanner := bufio.NewScanner(file)
	var atoms []atom
	structureCenter := []float64{1e10, 0.0, 0.0}
	for scanner.Scan() {
		// get next line
		line := scanner.Text()
		fields := strings.Fields(line)
		if len(fields) > 5 {
			newAtom := line2atom(fields)
			if fields[0] == "1" {
				structureCenter = []float64{newAtom.x, newAtom.y, newAtom.z}
			}
			if dist2center(structureCenter, newAtom) < shellDist {
				atoms = append(atoms, newAtom)
			}
		}
	}
	atoms = renumber(atoms)
	return atoms
}

func renumber(atoms []atom) []atom {
	intMap := make(map[int]int)
	for i := 0; i < len(atoms); i++ {
		intMap[atoms[i].id] = i + 1
	}
	for _, thisAtom := range atoms {
		for i := 0; i < len(thisAtom.bonds); i++ {
			thisAtom.bonds[i] = intMap[thisAtom.bonds[i]]
		}
	}
	return atoms
}

func dist2center(c []float64, a atom) float64 {
	dist := math.Sqrt(math.Pow(c[0]-a.x, 2) + math.Pow(c[1]-a.y, 2) + math.Pow(c[2]-a.z, 2))
	return dist
}

func line2atom(fields []string) atom {
	var newAtom atom
	newAtom.id, _ = strconv.Atoi(fields[0])
	newAtom.element = fields[1]
	newAtom.x, _ = strconv.ParseFloat(fields[2], 64)
	newAtom.y, _ = strconv.ParseFloat(fields[3], 64)
	newAtom.z, _ = strconv.ParseFloat(fields[4], 64)
	newAtom.atomType = fields[5]
	newAtom.bonds = []int{}
	for i := 6; i < len(fields); i++ {
		bondInt, _ := strconv.Atoi(fields[i])
		newAtom.bonds = append(newAtom.bonds, bondInt)
	}
	return newAtom
}

type atom struct {
	element  string
	x        float64
	y        float64
	z        float64
	bonds    []int
	atomType string
	id       int
	newID    int
}
