package main

import(
	"fmt"
	"github.com/biogo/hts/sam"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/biogo/alphabet"
	"os"
	"log"
	"flag"
	"strconv"
)

const (
    CaseSensitive = true
	RefType = true
	QueryType = false
)
/*
- This function get correction's sequence
- Input - (samRecord, refSeq, querySeq)
-			samRecord -- sam record
-			refSeq -- one reference sequence
-			querySeq -- one query sequence
- Output - correction's sequence
*/
func correction_seq(samRecord *sam.Record,
				    refSeq alphabet.Slice,
				    querySeq alphabet.Slice,
					correctType string) *linear.Seq {

	/* Detect strand */
	corSeq := linear.NewSeq("", []alphabet.Letter(""), alphabet.DNA)
	tempSeq := linear.NewSeq("", []alphabet.Letter(""), alphabet.DNA)
	if samRecord.Strand() ==  -1 {
		if querySeq, ok := querySeq.Slice(0,querySeq.Len()).(alphabet.Letters); ok{
			tempSeq.AppendLetters(querySeq...)
			tempSeq.RevComp()
		}
		querySeq = tempSeq.Seq.Slice(0,tempSeq.Seq.Len())
	}

	var tempSlice alphabet.Slice
	refPos, queryPos, numI, numD,numS, numH := samRecord.Start(), 0, 0, 0, 0, 0
	//FirstFlag - first cigar = true , otherwise = false
	FirstFlag := true
	TypeFlag := RefType
	if correctType == "query" {
		TypeFlag = QueryType
	}
	/* Parser Cigar */
	for _, cigar := range samRecord.Cigar {
		tempSlice, refPos, queryPos, numI, numD, numS, numH = fix_cigar(cigar, refSeq, querySeq, refPos, queryPos, numI, numD, numS, numH, FirstFlag, TypeFlag)
		if tempSlice != nil {
			if allSlice, ok := tempSlice.Slice(0,tempSlice.Len()).(alphabet.Letters); ok{
				corSeq.AppendLetters(allSlice...)
			}
		}
		FirstFlag = false
	}

	/* Get Sequence Name */
	queryRecordStart := strconv.Itoa(1+numS+numH)
	queryRecordEnd := strconv.Itoa(samRecord.End()+numI+numS+numH)
	refRecordStart := strconv.Itoa(samRecord.Start())
	refRecordEnd := strconv.Itoa(samRecord.End())
	corSeq.Annotation.ID = samRecord.Name+":"+queryRecordStart+"-"+queryRecordEnd+"|"+samRecord.Ref.Name()+":"+refRecordStart+"-"+refRecordEnd

	return corSeq
}


func fix_cigar(cigar sam.CigarOp,
			   refSeq alphabet.Slice,
			   querySeq alphabet.Slice,
			   refPos int, queryPos int,
			   numI int, numD int, numS int, numH int,
			   FirstFlag bool,TypeFlag bool) (alphabet.Slice, int, int, int, int, int, int) {

	var appSeq alphabet.Letters
	var seq alphabet.Slice
	switch cigar.Type().String() {
		case "M":
			if(TypeFlag == true){
				seq = appSeq.Append(refSeq.Slice(refPos, refPos+cigar.Len()))
			}else{
				seq = appSeq.Append(querySeq.Slice(queryPos, queryPos+cigar.Len()))
			}
				refPos += cigar.Len()
				queryPos += cigar.Len()
		case "I":
			seq = appSeq.Append(querySeq.Slice(queryPos, queryPos+cigar.Len()))
			queryPos += cigar.Len()
			numI += cigar.Len()
		case "D":
			seq = appSeq.Append(refSeq.Slice(refPos, refPos+cigar.Len()))
			refPos += cigar.Len()
			numD += cigar.Len()
		case "N":
		case "S":
			if(FirstFlag == true){
				queryPos += cigar.Len()
				numS += cigar.Len()
			}
		case "H":
			if(FirstFlag == true){
				queryPos += cigar.Len()
				numH += cigar.Len()
			}
		case "P":
		case "=":
		case "X":
		default:
			log.Printf("NO cigar: %s",cigar.Type().String())
	}
	return  seq, refPos, queryPos, numI, numD, numS, numH
}

func main(){

	/*
	-input file name
	-	-refFile="ref.fa"
	-	-querySAM="query.sam"
	*/
	var refFile string
	var querySAM string
	var queryFile string
	var correctType string
	flag.StringVar(&refFile, "refFile", "ref.fa","refFile")
	flag.StringVar(&queryFile, "queryFile", "query.fa","queryFile")
	flag.StringVar(&querySAM, "querySAM", "query.sam","querySAM")
	flag.StringVar(&correctType, "correctType", "ref or query --> use {ref} correct Sequence || use {query} correct Sequence(Default REFERENCE)","correctType")
	flag.Parse()


	/* Open Ref(FASTA) file */
	refFa , err := os.Open(refFile)
	if err != nil {
		log.Println(err)
	}
	defer refFa.Close()

	/* Open Query(FASTA) file */
	queryFa , err := os.Open(queryFile)
	if err != nil {
		log.Println(err)
	}
	defer queryFa.Close()

	/* Open SAM file */
	quSam , err := os.Open(querySAM)
	if err != nil {
		log.Println(err)
	}
	defer quSam.Close()

	/* Get the refFASTA content*/
	seqIndex := make(map[string]alphabet.Slice)
	var faReader *fasta.Reader
	faReader = fasta.NewReader(refFa, linear.NewSeq("", nil, alphabet.DNA))
	seqIO := seqio.NewScanner(faReader)
	for seqIO.Next() {
		seqFA := seqIO.Seq()
		seqIndex[seqFA.Name()] = seqFA.Slice()
	}

	/* Get the queryFASTA content*/
	querySeqIndex := make(map[string]alphabet.Slice)
	faReader = fasta.NewReader(queryFa, linear.NewSeq("", nil, alphabet.DNA))
	seqIO = seqio.NewScanner(faReader)
	for seqIO.Next() {
		seqFA := seqIO.Seq()
		querySeqIndex[seqFA.Name()] = seqFA.Slice()
	}

	/* Get SAM Header */
	var samRecord *sam.Record
	samLine, err := sam.NewReader(quSam)
	if err != nil{
		log.Println(err)
	}

	/* read Record */
	for samRecord, err = samLine.Read();samRecord != nil;samRecord, err = samLine.Read() {
		if err != nil && sam.IsValidRecord(samRecord) != true{
			log.Println(err)
		}else{
			corSeq := correction_seq(samRecord,seqIndex[samRecord.Ref.Name()],querySeqIndex[samRecord.Name],correctType)
			fmt.Printf(">%s\n",corSeq.Annotation.ID)
			fmt.Printf("%s\n",corSeq.Seq)
		}
	}
}
