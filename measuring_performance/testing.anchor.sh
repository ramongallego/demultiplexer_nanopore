OUTPUT_SUMMARY=/home/mk1b/Documents/test_demult/ouput.csv
printf "threshold,start,right.length,anchored,deplated,demultiplexed\n" \
> "${OUTPUT_SUMMARY}"
READ1=/data/minknow/data/Run003_20220520_longAmplicons/no_sample/20220520_1400_MN40189_FAT20945_29cfc02a/fast5_pass/fastq_sup/all.passing.fastq
DIRECTION=/home/mk1b/Documents/test_demult/demultiplexed_20220822_1524/
OUTPUT_DIR=/home/mk1b/Documents/test_demult/testing_anchors
START_TIME=$(date +%Y%m%d_%H%M)
OUTPUT_FOLDER="${OUTPUT_DIR}"/demultiplexed_"${START_TIME}"
mkdir "${OUTPUT_FOLDER}"
cd "${OUTPUT_FOLDER}"
cutadapt -m 1000 -M 2000 "${READ1}" -o "${READ1}".step1.fastq --quiet


echo "Number of sequences to begin with"

Start=$(grep ^@ "${READ1}" | wc -l)

echo "${Start}"
echo

echo "Number of sequences of the right length"

Good_length=$(grep ^@ "${READ1}".step1.fastq | wc -l)

echo "${Good_length}"
echo


awk 'NR % 4 == 2 {print length ($0)}' "${READ1}" > seq_lengths_start.txt
awk 'NR % 4 == 2 {print length ($0)}' "${READ1}".step1.fastq > seq_lengths_filter.txt

for threshold in 0.3; do # Keeping 0.3 from now onwards
  # Start by keeping only the reads of reasonable lengths



  echo "$threshold"
#    cutadapt -g file:"${DIRECTION}"/direction.fasta -o "${READ1}".new.fastq --quiet --action=none \
#   --rc -e "$threshold" "${READ1}".step1.fastq --discard-untrimmed



  echo "Number of sequences after"

#  Direction_check=$(grep ^@ "${READ1}".new.fastq | wc -l)

  echo "Number of sequences reversed"

#  Reversed=$(grep rc "${READ1}".new.fastq | wc -l)

  echo "Adding an anchor cutadapt so we are all on the same spot"

#  cutadapt -g AATGATACGGCGACCACCGAGATCTACAC -o "${READ1}".two.step.fastq \
#  --quiet --discard-untrimmed -e "$threshold" "${READ1}".new.fastq

  echo "Number of lines after anchoring"

#  Anchored=$(grep ^@ "${READ1}".two.step.fastq | wc -l )

  echo "What about anchoring and reversing at the same time"

  cutadapt -g AATGATACGGCGACCACCGAGATCTACAC -o "${READ1}".one.step.fastq \
  --quiet --discard-untrimmed -e "$threshold" "${READ1}".step1.fastq --rc

  echo "Number of lines after anchoring directly"

  Direct_Anchor=$(grep ^@ "${READ1}".one.step.fastq | wc -l )





# Demultiplex fwd and calculate success rate
echo "Deplating the two step at ${threshold} tolerance"
echo
#cutadapt -g file:"${DIRECTION}"/barcodes_P5_anchored.fasta -o two.step_"${threshold}"_{name}_round1.fastq \
# --quiet --untrimmed-output anchored_nop5.fastq -e "${threshold}" "${READ1}".two.step.fastq
echo "Deplating the one step at ${threshold} tolerance"
echo
 cutadapt -g file:"${DIRECTION}"/barcodes_P5_anchored.fasta -o one.step_"${threshold}"_{name}_round1.fastq \
  --quiet --untrimmed-output directly_nop5.fastq -e "${threshold}" "${READ1}".one.step.fastq
# A FOR LOOP TO GET THE STATS FOR EACH DEMULTIPLEXING
tot_demult_anchored=0
tot_demult_direct=0
tot_deplat_anchored=0
tot_deplat_direct=0
tot_demult_anchored_back=0

#  for file in  two.step_"${threshold}"_*_round1.fastq; do

#    basename=$(echo "${file}" | sed 's/_round1.fastq//g'  )

#     Seqs_deplat=$(awk 'NR %4 == 2' $file | wc -l)
#     tot_deplat_anchored=$(expr "${tot_deplat_anchored}" + "${Seqs_deplat}")
#     awk 'NR % 4 == 2 {print length ($0)}' "${file}" >> seq_lengths_deplated_anchored_direction_"${basename}".txt

#     cutadapt -a ATCTCGTATGCCGTCTTCTGCAAGCAGAAGACGGCATACGAGATCTTG \
#              -o  "${file}".short.fastq \
#     	 "${file}"  --discard-untrimmed --quiet -e "${threshold}"



#     cutadapt -a file:"${DIRECTION}"/noprimers/all.passing.fastq/A7/barcodes.p7.strict.fasta \
#      -o "${basename}"_{name}_round2.fastq \
#      --quiet --untrimmed-output directly_nop7.fastq -e "${threshold}" "${file}".short.fastq

#      Seqs_demult=$(awk 'NR %4 == 2' "${basename}"_*_round2.fastq | wc -l)
#      tot_demult_ancCAAGCAGAAGACGGCATACGAGAThored=$(expr "${tot_demult_anchored}" + "${Seqs_demult}")
#      awk 'NR % 4 == 2 {print length ($0)}' "${basename}"_*_round2.fastq >> seq_lengths_demulted_anchored_two.step_"${basename}".txt

#   done
   #Now direct
   for file in  one.step_"${threshold}"_*_round1.fastq; do

     basename=$(echo "${file}" | sed 's/_round1.fastq//g'  )

      Seqs_deplat=$(awk 'NR %4 == 2' $file | wc -l)
      tot_deplat_direct=$(expr "${tot_deplat_direct}" + "${Seqs_deplat}")
      awk 'NR % 4 == 2 {print length ($0)}' "${file}" >> seq_lengths_deplated_anchored_one.step_"${basename}".txt
#      cutadapt -a ATCTCGTATGCCGTCTTCTGCTTG \
#               -o  "${file}".short.fastq \
#        "${file}"  --discard-untrimmed --quiet -e "${threshold}"



#      cutadapt -a file:"${DIRECTION}"/noprimers/all.passing.fastq/A7/barcodes.p7.strict.fasta \
#       -o "${basename}"_{name}_round2.fastq \
#       --quiet --untrimmed-output directly_nop7.fastq -e "${threshold}" "${file}".short.fastq

#       Seqs_demult=$(awk 'NR %4 == 2' "${basename}"_*_round2.fastq | wc -l)
#       tot_demult_direct=$(expr "${tot_demult_direct}" + "${Seqs_demult}")
#       awk 'NR % 4 == 2 {print length ($0)}' "${basename}"_*_round2.fastq >> seq_lengths_demulted_anchored_one.step_"${basename}".txt

## Seems like the second demultiplexing step is the worst performer, with only 55-60% success. Lets do it with anchors as well

seqkit seq -r -p -t DNA "${file}" -o "${basename}".rc.fastq


  cutadapt -g CAAGCAGAAGACGGCATACGAGAT \
         -o  "${basename}".anchored.rc.fastq \
         "${basename}".rc.fastq  --discard-untrimmed --quiet -e 0.4 # Change adaptor anchoring -e value

  cutadapt -g file:"${DIRECTION}"/noprimers/all.passing.fastq/A7/barcodes.p7.rc.anchored.fasta \
  -o "${basename}"_{name}_round2_rc.fastq --no-indels \
  --quiet --untrimmed-output "${basename}"_nop7.fastq -e "${threshold}" "${basename}".anchored.rc.fastq

  Seqs_demult=$(awk 'NR %4 == 2' "${basename}"_*_round2_rc.fastq | wc -l)
  tot_demult_anchored_back=$(expr "${tot_demult_anchored_back}" + "${Seqs_demult}")
  awk 'NR % 4 == 2 {print length ($0)}' "${basename}"_*_round2_rc.fastq >> seq_lengths_demulted_anchored_one.step_rc_"${basename}".txt

    done

    printf "%s,%s,%s,%s,%s,%s\n" \
    "${threshold}" "${Start}" "${Good_length}" "${Direct_Anchor}" \
    "${tot_deplat_direct}" \
    "${tot_demult_anchored_back}">> "${OUTPUT_SUMMARY}"

    # Length distributions

  #  awk 'NR % 4 == 2 {print length ($0)}' "${READ1}".new.fastq > seq_lengths_direction_"${threshold}".txt
  #  awk 'NR % 4 == 2 {print length ($0)}' "${READ1}".anchored.fastq > seq_lengths_achored_"${threshold}".txt
    awk 'NR % 4 == 2 {print length ($0)}' "${READ1}".anchored.directly.fastq > seq_lengths_achored_directly_"${threshold}".txt

done
