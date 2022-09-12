## A script to extract the lengths of the sequences after each process
## Usage: bash seq.lengths.sh <path/to/output>

OUTPUTFOLDER=$1

cd "$OUTPUTFOLDER"/demultiplexed

for file in */*_round1.fastq; do
awk -v VARIABLE=$file 'NR%4==2 {print f,length, VARIABLE } {f=$1}' $file >> "$OUTPUTFOLDER"/seq.lengths.demult.by.plate.txt
 done

cd ../noprimers

 for file in */*_Well_*.fastq; do awk -v VARIABLE=$file 'NR%4==2 {print f,length, VARIABLE } {f=$1}' $file >> "$OUTPUTFOLDER"/seq.lengths.demult.by.well.txt ; done

# cd "$OUTPUTFOLDER"/noprimers
#


# for file in */*/*.fastq; do awk -v VARIABLE=$file 'NR%4==2 {print f,length, VARIABLE } {f=$1}' $file >> "$OUTPUTFOLDER"/seq.lengths.demult.primer.removed.txt ; done

cd ../demultiplexed

for file in */*_adap.fastq; do
awk -v VARIABLE=$file 'NR%4==2  {print f,length, VARIABLE } {f=$1}' $file >> "$OUTPUTFOLDER"/seq.lengths.adapters.txt

awk -v VARIABLE=$file 'NR%4==2 {print f,substr($1, 1,10), length, VARIABLE } {f=$1}' $file >> "$OUTPUTFOLDER"/seq.adapters.txt
 done

for file in */*.anchored.rc.fastq; do
awk -v VARIABLE=$file 'NR%4==2 {print f,length, VARIABLE } {f=$1}' $file >> "$OUTPUTFOLDER"/seq.lengths.after.rev.primer.txt
 done