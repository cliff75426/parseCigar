package countNM

import(
	"fmt"
	"github.com/biogo/hts/sam"
	"os"
	"log"
	"flag"
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
	NMtag := sam.NewTag("NM")
	count := 0
	for samRecord, err = samLine.Read();samRecord != nil;samRecord, err = samLine.Read() {
		if samRecord.AuxFields.Get(NMtag) != nil{
			if intcount,ok := samRecord.AuxFields.Get(NMtag).Value().(uint16); ok{
				count += int(intcount)
			}else if intcount,ok := samRecord.AuxFields.Get(NMtag).Value().(uint8); ok{
				count += int(intcount)
			}
			
		}
	}
	fmt.Println(count)
}
