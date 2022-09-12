#!/bin/bash
# Start by activating decona
# Usage bash demultiplex_both_fastqs.sh banzai_params.sh ...
#This script is built using banzai (github.com/jimmyodonnell/banzai) as template

##### LOCATING SCRIPTS AND SOURCING

MAIN_DIR="$(dirname "$0")"
Current_dir=$(pwd)
SCRIPT_DIR="${MAIN_DIR}"/scripts
for file in "${SCRIPT_DIR}"/*.sh ; do
	echo "sourcing script ${file}"
	source "${file}"
done

###### LOCATING PARAMETERS FILE

param_file=${1}

echo "Reading analysis parameters from:"
echo "${param_file}"
source "${param_file}"

###### LOCATING METADATA

# Check if the metadata file exists
if [[ -s "${SEQUENCING_METADATA}" ]]; then
	echo "Reading metadata from:"
	echo "${SEQUENCING_METADATA}"
else
	echo 'ERROR! Could not find metadata file. You specified the file path:'
	echo
	echo "${SEQUENCING_METADATA}"
	echo
	echo 'That file is empty or does not exist. Aborting script.'
	exit
fi


# Now fix line ends if needed

if [[ $( file "${SEQUENCING_METADATA}" ) == *"CRLF"* ]]; then

  echo "The file has CRLF endings. Let me fix that for you..."

  BASE="${SEQUENCING_METADATA%.*}"

  EXT="${SEQUENCING_METADATA##*.}"

  NEWLINES_FIXED="${BASE}"_fix."${EXT}"

  tr -d '\r' < "${SEQUENCING_METADATA}" > "${NEWLINES_FIXED}"

  echo "the old file was: ${SEQUENCING_METADATA}"

  echo "The new file is here:"

  echo "${NEWLINES_FIXED}"

else

  echo "The file passes test for CRLF. Everybody dance!"
  echo

fi

if [[ -s "${NEWLINES_FIXED}" ]]; then
	SEQUENCING_METADATA="${NEWLINES_FIXED}"
fi

###### CREATING OUTPUT DIRECTORY

START_TIME=$(date +%Y%m%d_%H%M)
OUTPUT_DIR="${OUTPUT_DIRECTORY}"/demultiplexed_"${START_TIME}"

mkdir "${OUTPUT_DIR}"
echo "Output directory is ${OUTPUT_DIR}"
echo

# copy metadata and parameters file to output directory
cp "${SEQUENCING_METADATA}" "${OUTPUT_DIR}"/metadata.csv
SEQUENCING_METADATA="${OUTPUT_DIR}"/metadata.csv
cp "${param_file}" "${OUTPUT_DIR}"/banzai_params.sh

# Write a log file
LOGFILE="${OUTPUT_DIR}"/logfile.txt
exec > >(tee "${LOGFILE}") 2>&1

FINAL_DIR="${OUTPUT_DIR}"/noprimers
mkdir "${FINAL_DIR}"

DEMULT_DIR="${OUTPUT_DIR}"/demultiplexed
mkdir "${DEMULT_DIR}"

BARCODES_DIR="${OUTPUT_DIR}"/barcodes
mkdir "${BARCODES_DIR}"

################################################################################
# READ METADATA
################################################################################
# report metadata dimensions
METADATA_DIM=($( awk -F, 'END{print NR, NF}' "${SEQUENCING_METADATA}" ))

echo " #### Reading metadata file"
echo ""

echo "Metadata has" "${METADATA_DIM[0]}" "rows and" "${METADATA_DIM[1]}" "columns including header."
N_SAMPLES=$( echo "${METADATA_DIM[0]}" - 1 | bc )
echo "Expecting" "${N_SAMPLES}" "samples total."
echo
## NOW WE HAVE LOADED THE SEQUENCING_METADATA - WE NEED to find the columns specified
## in the params file. We should set up an alert & quit if a critical column is not found

# Filenames
COLNUM_FILE1=$( get_colnum "${COLNAME_FILE1}" "${SEQUENCING_METADATA}")



# PLATE names and SEQS
COLNUM_ID1=$( get_colnum "${COLNAME_ID1_NAME}" "${SEQUENCING_METADATA}")

COLNUM_ID1_SEQ=$( get_colnum "${COLNAME_ID1_SEQ}" "${SEQUENCING_METADATA}")

# WELL Names and SEQS
COLNUM_ID2_SEQ=$( get_colnum "${COLNAME_ID2_SEQ}" "${SEQUENCING_METADATA}")
COLNUM_ID2_WELL=$(get_colnum "${COLNAME_ID2_WELL}" "${SEQUENCING_METADATA}")

# Sample names
COLNUM_SAMPLE=$( get_colnum "${COLNAME_SAMPLE_ID}" "${SEQUENCING_METADATA}")

# Primers
COLNUM_PRIMER1=$( get_colnum "${COLNAME_PRIMER1}" "${SEQUENCING_METADATA}")
COLNUM_PRIMER2=$( get_colnum "${COLNAME_PRIMER2}" "${SEQUENCING_METADATA}")
COLNUM_LOCUS=$( get_colnum "${COLNAME_LOCUS}" "${SEQUENCING_METADATA}")

# Run away from the script if any of the previous columns was not found

all_columns=( COLNUM_FILE1 COLNUM_ID1 COLNUM_ID1_SEQ COLNUM_ID2_SEQ COLNUM_ID2_WELL \
COLNUM_SAMPLE COLNUM_PRIMER1 COLNUM_PRIMER2 COLNUM_LOCUS)

echo "Checking that all columns in metadata are there"

for column in "${all_columns[@]}" ; do

 if [ "${!column}"  > 0 ]; then
	 echo "looking good, ${column}"
 else
  echo "Something went wrong with column name ${column}"
	echo "exiting script"
	exit
fi
done
echo "All columns passed test"

#############
# ADD Reverse complements of Well indices and PCR primers
#############

	 ## First rc the i7, then join it to the metadata

awk -F',' -v COLNUM="${COLNUM_ID2_SEQ}"    \
'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" | \
while read line ; do revcom $line ; done |\
sed '1s/^/RC_i7\n/'  > "${OUTPUT_DIR}"/rci7.txt

paste "${SEQUENCING_METADATA}" "${OUTPUT_DIR}"/rci7.txt -d ',' > "${OUTPUT_DIR}"/metadata.new.csv
 rm "${OUTPUT_DIR}"/rci7.txt
   ## First rc the PCR PRIMER, then join it to the metadata


awk -F',' -v COLNUM="${COLNUM_PRIMER2}"    \
'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" | \
while read line ; do revcom $line ; done |\
sed '1s/^/RC_PCR_REV\n/'  > "${OUTPUT_DIR}"/rc_pcr.txt

paste "${OUTPUT_DIR}"/metadata.new.csv "${OUTPUT_DIR}"/rc_pcr.txt -d ',' > "${OUTPUT_DIR}"/metadata.new.new.csv

rm "${OUTPUT_DIR}"/rc_pcr.txt

mv "${OUTPUT_DIR}"/metadata.new.new.csv "${SEQUENCING_METADATA}"

rm "${OUTPUT_DIR}"/metadata.new.csv


##LOCATE the COLUMN OF the RCs
COLNUM_PRIMER2_RC=$( get_colnum "RC_PCR_REV" "${SEQUENCING_METADATA}")
COLNUM_ID2_RCSEQ=$( get_colnum "RC_i7" "${SEQUENCING_METADATA}")
## Update metadata DIM

METADATA_DIM=($( awk -F, 'END{print NR, NF}' "${SEQUENCING_METADATA}" ))

echo "Metadata has" "${METADATA_DIM[0]}" "rows and" "${METADATA_DIM[1]}" "columns including header."


################################################################################
# CHECK FILES
################################################################################


	FILE1=($(awk -F',' -v COLNUM=$COLNUM_FILE1 \
	  'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA |\
	  sort | uniq))



	NFILE1="${#FILE1[@]}"


	  echo 'Files read from metadata columns' "${COLNUM_FILE1}"
	  echo 'File names:'
		for (( i=0; i < "${NFILE1}"; ++i)); do
			printf '%s\t%s\n' "${FILE1[i]}"
		done
		echo

    echo "Primary Indices read from column ${COLNAME_ID1_SEQ}"

	ID1_ALL=($(awk -F',' -v COLNAME=$COLNUM_ID1 -v COLSEQ=$COLNUM_ID1_SEQ \
	  'NR>1 { print $COLNAME, $COLSEQ }' "${SEQUENCING_METADATA}" |\
			sort | uniq))
	NID1="${#ID1_ALL[@]}"

	echo "There are ${NID1} primary indices"

	awk -F',' -v COLNAME=$COLNUM_ID1 -v COLSEQ=$COLNUM_ID1_SEQ \
	  'NR>1 { print $COLNAME, $COLSEQ }' $SEQUENCING_METADATA |sort | uniq |\
		awk '{printf ">%s\n%s\n", $1, $2}' >> "${BARCODES_DIR}"/barcodes_P5.fasta

cat "${OUTPUT_DIR}"/barcodes_P5.fasta

	# awk -F',' -v COLNAME=$COLNUM_ID1 -v COLSEQ=$COLNUM_ID1_SEQ \
	#   'NR>1 { print $COLNAME, $COLSEQ }' $SEQUENCING_METADATA |sort | uniq |\
	# 	awk '{printf ">%s\n^%s\n", $1, $2}' >> "${BARCODES_DIR}"/barcodes_P5_anchored.fasta


################################################################################
# Read in primers
################################################################################
# Modify this so in case there are multiple primer pairs, they are not artificially paired
# alphabetically

	PRIMER1=($(awk -F',' -v COLNUM=$COLNUM_PRIMER1 \
	  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA ))

	PRIMER2=($(awk -F',' -v COLNUM=$COLNUM_PRIMER2 \
	  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA ))

	LOCI=($(awk -F',' -v COLNUM=$COLNUM_LOCUS \
	  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA ))

	PRIMER2_RC=($(awk -F',' -v COLNUM=$COLNUM_PRIMER2_RC \
	  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA ))

	if [[ -n "${PRIMER1}" && -n "${PRIMER2}" ]]; then
	  echo 'Primers read from metadata columns' "${COLNUM_PRIMER1}" 'and' "${COLNUM_PRIMER2}"
	  echo 'Primer sequences:' "${PRIMER1}" "${PRIMER2}"
		echo
		echo 'Reverse Complemented ' "${PRIMER2}" 'into' "${PRIMER2_RC}"
		echo
	else
	  echo 'ERROR:' 'At least one primer is not valid'
	  echo 'Looked in metadata columns' "${COLNUM_PRIMER1}" 'and' "${COLNUM_PRIMER2}"
	  echo 'Aborting script'
	  exit
	fi

##############################################
# Create a fasta with all Fwds to get all sequences in the same direction
############################################

awk -F','  -v FWD=$COLNUM_PRIMER1 -v LOCI=$COLNUM_LOCUS \
	  ' NR > 1 { print $LOCI , $FWD }' $SEQUENCING_METADATA | sort | uniq  |\
	  awk '{printf ">%s\n%s\n", $1, $2}' >> "${BARCODES_DIR}"/direction.fasta

		awk -F','  -v FWD=$COLNUM_PRIMER2 -v LOCI=$COLNUM_LOCUS \
			  ' NR > 1 { print $LOCI , $FWD }' $SEQUENCING_METADATA | sort | uniq  |\
			  awk '{printf ">%s\n%s\n", $1, $2}' >> "${BARCODES_DIR}"/revprimer.fasta
#######
#Unique samples are given by combining the primary and secondary indexes
######
	ID_COMBO=$( awk -F',' -v COLNUM1=$COLNUM_ID1 -v COLNUM2=$COLNUM_ID2_WELL \
	'NR>1 {
	  print "Plate=" $COLNUM1 ";Well=" $COLNUM2
	}' "${SEQUENCING_METADATA}" )

	SAMPLE_NAMES=($(awk -F',' -v COLNUM=$COLNUM_SAMPLE \
	  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" ))



	OUTPUT_SUMMARY="${OUTPUT_DIR}/summary.csv"
	printf "Sample,locus,plated,demultiplexed,noprimers\n" \
	> "${OUTPUT_SUMMARY}"

################################################################################
# BEGIN LOOP TO PERFORM LIBRARY-LEVEL ACTIONS
################################################################################

	for (( i=0; i < "${#FILE1[@]}"; i++ )); do
	  # Identify the fastq files. It is usually one, but just in case

	  READ1="${PARENT_DIR}/${FILE1[i]}"

	  BASE1="${FILE1[i]%.*}"


	  mkdir "${DEMULT_DIR}"/"${FILE1[i]}"
		mkdir "${FINAL_DIR}"/"${FILE1[i]}"



		echo "Working on input file $[i+1] out of ${#FILE1[@]}"

	##First cutdapt:


	# Look for the fwd primer and if found on rc, reverse the output. New version: Do it
	# by using the anchor sequence

	echo ""
	echo "Using cutadapt to select by sequence length, reorient sequences and trim them to the anchor"
	echo "this might take a while"
	# Did it using a pipe cutadapt as it is a notch quicker to filter first and find and reverse later
	cutadapt -m "${MIN_LENGTH}" -M "${MAX_LENGTH}" "${READ1}" --cores="${N_CORES}" -o "${READ1}".goodlength.fastq

	 cutadapt -g "file:"${BARCODES_DIR}"/direction.fasta;min_overlap=5" -o "${READ1}".anchored.fastq \
	--quiet --discard-untrimmed -e 0.2 --rc --cores="${N_CORES}" --rename='{id} locus_{adapter_name} {comment}' "${READ1}".goodlength.fastq
# Secind cutadapt to keep only tentative adapters
 cutadapt -a "file:"${BARCODES_DIR}"/direction.fasta;min_overlap=5" -o "${READ1}".tentative_i5.fastq \
--quiet --discard-untrimmed -e 0.2 "${READ1}".goodlength.fastq --rc --cores="${N_CORES}"

  echo "Number of sequences to begin with"

  awk 'NR % 4 == 2' "${READ1}" | wc -l

  echo "Number of sequences after size selection and anchoring"

  awk 'NR % 4 == 2' "${READ1}".anchored.fastq | wc -l

	echo "Number of sequences reversed"

	grep rc "${READ1}".anchored.fastq | wc -l

  echo ""


# Cutadapt on the short file
    cutadapt -g file:"${BARCODES_DIR}"/barcodes_P5.fasta -o "${DEMULT_DIR}"/"${FILE1[i]}"/{name}_adap.fastq \
   --quiet --untrimmed-output "${OUTPUT_DIR}"/${BASE1}_nop5.fastq -e 0.1 "${READ1}".tentative_i5.fastq --cores="${N_CORES}"

	# cutadapt -g file:"${OUTPUT_DIR}"/barcodes_P5.fasta -o "${DEMULT_DIR}"/"${FILE1[i]}"/{name}_round1.fastq \
	#  --quiet --untrimmed-output "${OUTPUT_DIR}"/${BASE1}_nop5.fastq -e 0.2 "${READ1}".new.fastq




		n_files=("${DEMULT_DIR}"/"${FILE1[i]}"/*_adap.fastq)





	 ls "${DEMULT_DIR}"/"${FILE1[i]}"/

	 for file in "${n_files[@]}"; do

		 echo "working on plate " "${file}"



	 BASE_OUTPUT=$(basename "${file}" |  sed 's/_adap.fastq//g')

	 echo "${BASE_OUTPUT}"

## get in awk the list of ids that I need
awk 'NR %4==1 {print substr($1,2)}' ${file} > "${OUTPUT_DIR}"/ids.txt
echo
seqkit grep -f "${OUTPUT_DIR}"/ids.txt "${READ1}".anchored.fastq > "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"_round1.fastq
	 ## subset in awk
echo
# finishing loops here

# done
# done

awk -F',' -v COLNAME="${COLNUM_ID1}" -v VALUE="${BASE_OUTPUT}" \
	 -v COLSEQ2="${COLNUM_ID2_SEQ}" -v COLID2="${COLNUM_ID2_WELL}" \
	 'NR>1 { if ($COLNAME == VALUE) {printf ">Well_%s\n^%s\n", $COLID2, $COLSEQ2} }'  "${SEQUENCING_METADATA}" > "${BARCODES_DIR}"/barcodes.p7.for.rc.fasta

	 awk -F',' -v COLNAME="${COLNUM_ID1}" -v VALUE="${BASE_OUTPUT}" \
	 	 -v COLSEQ2="${COLNUM_ID2_SEQ}" -v COLID2="${COLNUM_ID2_WELL}" \
	 	 'NR>1 { if ($COLNAME == VALUE) {printf ">Well_%s\n%s\n", $COLID2, $COLSEQ2} }'  "${SEQUENCING_METADATA}" > "${BARCODES_DIR}"/barcodes.p7.for.rc.fasta


#awk -F',' -v COLNAME="${COLNUM_ID1}" -v VALUE="${BASE_OUTPUT}" \
#-v COLSEQ2="${COLNUM_ID2_RCSEQ}" -v COLID2="${COLNUM_ID2_WELL}" \
#'NR>1 { if ($COLNAME == VALUE) {printf ">Well_%s\n%s\n", $COLID2, $COLSEQ2} }'  "${SEQUENCING_METADATA}" > "${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/barcodes.p7.fasta

#awk -F',' -v COLNAME="${COLNUM_ID1}" -v VALUE="${BASE_OUTPUT}" \
#-v COLSEQ2="${COLNUM_ID2_RCSEQ}" -v COLID2="${COLNUM_ID2_WELL}" \
#'NR>1 { if ($COLNAME == VALUE) {printf ">Well_%s\n%s$\n", $COLID2, $COLSEQ2} }'  "${SEQUENCING_METADATA}" > "${FINAL_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"/barcodes.p7.strict.fasta


## How does it work if we rc the files prior to finding the p7
## TODO: CAAGCAGAAGACGGCATACGAGAT can works as a way of anchoring the p7
	seqkit seq -r -p -t DNA "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"_round1.fastq \
	 -o "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"_round1.rc.fastq --quiet

	cutadapt -g "file:"${BARCODES_DIR}"/revprimer.fasta;min_overlap=5" \
         -o  "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}".anchored.rc.fastq \
         "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"_round1.rc.fastq \
				   --discard-untrimmed  -e 0.2  --cores="${N_CORES}" >> "${LOGFILE}"

					 cutadapt -a "file:"${BARCODES_DIR}"/revprimer.fasta;min_overlap=5" \
				          -o  "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}".adapters.fastq \
				          "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"_round1.rc.fastq \
				 				   --discard-untrimmed  -e 0.2 --cores="${N_CORES}" >> "${LOGFILE}"
				 # Change adaptor anchoring -e value



	cutadapt -g "file:"${BARCODES_DIR}"/barcodes.p7.for.rc.fasta;min_overlap=5" \
				 -o "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"_{name}_round2_rc.fastq --no-indels \
				 --quiet --untrimmed-output "${BASE_OUTPUT}"_nop7.fastq -e 0.1 "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}".adapters.fastq --cores="${N_CORES}" >> "${LOGFILE}"




### NOW FIND HOW MANY PRIMERS PER Plate-Well combo

n_files2=("${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}"_*_round2_rc.fastq)

	for file2 in "${n_files2[@]}"; do
FINAL_BASENAME=$(basename "${file2}" |  sed 's/_round2_rc.fastq//g')
	#	echo "${file2}"


	awk 'NR %4==1 {print substr($1,2)}' ${file2} > "${OUTPUT_DIR}"/ids.txt
	seqkit grep -f "${OUTPUT_DIR}"/ids.txt "${DEMULT_DIR}"/"${FILE1[i]}"/"${BASE_OUTPUT}".anchored.rc.fastq > "${DEMULT_DIR}"/"${FILE1[i]}"/"${FINAL_BASENAME}"_rc.fastq


		newp7name=$(echo "${file2}" | sed 's/_round2_rc.fastq$/.fastq/g')
		pcr_basename=$(basename "${newp7name}")



		seqkit seq -r -p -t DNA "${DEMULT_DIR}"/"${FILE1[i]}"/"${FINAL_BASENAME}"_rc.fastq \
		 -o "${FINAL_DIR}"/"${FILE1[i]}"/"${FINAL_BASENAME}".fastq --quiet

	  BASE_P7=$(echo "${pcr_basename}" |  sed 's/.fastq//g' | sed 's/.*Well_//g' )

   echo -ne "In plate ${BASE_OUTPUT} and Well  ${BASE_P7}"'\r'

if [[ "${DEMULT_BY_PRIMER}" == "YES" ]]; then

## This is fast, bc it only looks at the metadata, but misses unexpected
## matches - an alternative is to awk the second field

## awk 'NR % 4 ==1 && /locus_/{print $2}'"${FINAL_DIR}"/"${FILE1[i]}"/"${FINAL_BASENAME}".fastq

	awk -F',' -v COLNAME="${COLNUM_ID1}" -v VALUE="${BASE_OUTPUT}" \
	-v WELL_COL="${COLNUM_ID2_WELL}" -v WELL_NOW="${BASE_P7}" \
		 -v LOCI="${COLNUM_LOCUS}" \
		 'NR>1 { if ($COLNAME == VALUE && $WELL_COL == WELL_NOW) {printf "locus_%s\n", $LOCI} }'  "${SEQUENCING_METADATA}" | \
		 while read line ; do
		 # check if directory exists - if not, make it
if [[ ! -d "${FINAL_DIR}"/"$line" ]]; then
	mkdir "${FINAL_DIR}"/"$line"
fi

		 seqkit grep -n -r -p $line "${FINAL_DIR}"/"${FILE1[i]}"/"${FINAL_BASENAME}".fastq > \
       "${FINAL_DIR}"/"$line"/"${FILE1[i]}"_"${FINAL_BASENAME}"_"$line".fastq
		  done
#To avoid cluttering, remove the fastq before splitting by loci

rm  "${FINAL_DIR}"/"${FILE1[i]}"/"${FINAL_BASENAME}".fastq

fi


	nseq_deplated=$(cat ${file} | wc -l)

	nseq_demult=$(cat ${file2} | wc -l)

	nseq_noprimer=$(cat "${FINAL_DIR}"/"${FILE1[i]}"/"${FINAL_BASENAME}".fastq | wc -l )

	 printf "%s,%s,%s,%s,%s\n" \
	 "${BASE_OUTPUT}" "${pcr_basename}" \
	 "${nseq_deplated}" "${nseq_demult}" "${nseq_noprimer}" >> "${OUTPUT_SUMMARY}"



	done # End of loop across all wells within a plate




done # End of loop across all plates within a file
# AN extra thing to do is to create a csv file with all the sequence headers that ended on each sample - so we can run decona with all samples and
# demultiplex again
for file in "${FINAL_DIR}"/"${FILE1[i]}"/*.fastq; do
	 awk 'BEGIN {OFS=","} NR %4==1 {printf "%s,%s\n", FILENAME, $1}' $file >> "${OUTPUT_DIR}"/hash.list.csv


  done

	done # End of loop across all initial fastq files


if [[ "${HOARD}" != "YES" ]]; then
	rm "${OUTPUT_DIR}"/*nop5.fastq
	rm -rd "${DEMULT_DIR}"
fi
