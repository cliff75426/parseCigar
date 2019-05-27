package countInDel

import(
	"fmt"
	"github.com/biogo/hts/sam"
	"os"
	"log"
	"flag"
//	"strconv"
)


func main(){

	var querySAM string
	flag.StringVar(&querySAM, "querySAM", "query.sam","querySAM")
	flag.Parse()



	/* Open SAM file */
	quSam , err := os.Open(querySAM)
	if err != nil {
		log.Println(err)
	}
	defer quSam.Close()



	/* Get SAM Header */
	var samRecord *sam.Record
	samLine, err := sam.NewReader(quSam)

	/* read Record */
	count := 0
	for samRecord, err = samLine.Read();samRecord != nil;samRecord, err = samLine.Read() {
		for _, cigar := range samRecord.Cigar{
			switch cigar.Type().String() {
				case "I":
					count += cigar.Len()
				case "D":
					count += cigar.Len()
				default:
			}

		}
	}
	fmt.Println("InDel COUNT:")
	fmt.Println(count)
}
