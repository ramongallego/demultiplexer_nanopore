OUTPUT_SUMMARY=/home/mk1b/Documents/test_demult/ouput.csv
printf "threshold,start,right.length,direction,reversed,anchored,
direct_anchor,deplated_direction,deplated_direct, demult_twostep, demult_onestep\n" \
> "${OUTPUT_SUMMARY}"
READ1=/data/minknow/data/Run003_20220520_longAmplicons/no_sample/20220520_1400_MN40189_FAT20945_29cfc02a/fast5_pass/fastq_sup/all.passing.fastq
DIRECTION=/home/mk1b/Documents/test_demult/demultiplexed_20220820_1059/
OUTPUT_FOLDER=/home/mk1b/Documents/test_demult/testing_anchors
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

for threshold in 0.2 0.3 0.35; do
  # Start by keeping only the reads of reasonable lengths



  echo "$threshold"
    cutadapt -g file:"${DIRECTION}"/direction.fasta -o "${READ1}".new.fastq --quiet --action=none \
   --rc -e "$threshold" "${READ1}".step1.fastq --discard-untrimmed



  echo "Number of sequences after"

  Direction_check=$(grep ^@ "${READ1}".new.fastq | wc -l)

  echo "Number of sequences reversed"

  Reversed=$(grep rc "${READ1}".new.fastq | wc -l)

  echo "Adding an anchor cutadapt so we are all on the same spot"

  cutadapt -g AATGATACGGCGACCACCGAGATCTACAC -o "${READ1}".two.step.fastq \
  --quiet --discard-untrimmed -e "$threshold" "${READ1}".new.fastq

  echo "Number of lines after anchoring"

  Anchored=$(grep ^@ "${READ1}".two.step.fastq | wc -l )

  echo "What about anchoring and reversing at the same time"

  cutadapt -g AATGATACGGCGACCACCGAGATCTACAC -o "${READ1}".one.step.fastq \
  --quiet --discard-untrimmed -e "$threshold" "${READ1}".step1.fastq --rc

  echo "Number of lines after anchoring directly"

  Direct_Anchor=$(grep ^@ "${READ1}".one.step.fastq | wc -l )





# Demultiplex fwd and calculate success rate
echo "Deplating the two step at ${threshold} tolerance"
echo
cutadapt -g file:"${DIRECTION}"/barcodes_P5_anchored.fasta -o two.step_"${threshold}"_{name}_round1.fastq \
 --quiet --untrimmed-output anchored_nop5.fastq -e "${threshold}" "${READ1}".two.step.fastq
echo "Deplating the one step at ${threshold} tolerance"
echo
 cutadapt -g file:"${DIRECTION}"/barcodes_P5_anchored.fasta -o one.step_"${threshold}"_{name}_round1.fastq \
  --quiet --untrimmed-output directly_nop5.fastq -e "${threshold}" "${READ1}".one.step.fastq
# A FOR LOOP TO GET THE STATS FOR EACH DEMULTIPLEXING
tot_demult_anchored=0
tot_demult_direct=0
tot_deplat_anchored=0
tot_deplat_direct=0

  for file in  two.step_"${threshold}"_*_round1.fastq; do

     Seqs_deplat=$(awk 'NR %4 == 2' $file | wc -l)
     tot_deplat_anchored=$(expr "${tot_deplat_anchored}" + "${Seqs_deplat}")
     awk 'NR % 4 == 2 {print length ($0)}' "${file}" > seq_lengths_deplated_anchored_direction_"${file}"_"${threshold}".txt

     cutadapt -a ATCTCGTATGCCGTCTTCTGCTTG \
              -o  "${file}".short.fastq \
     	 "${file}"  --discard-untrimmed --quiet -e "${threshold}"



     cutadapt -a file:"${DIRECTION}"/noprimers/all.passing.fastq/A7/barcodes.p7.strict.fasta \
      -o "${file}"_{name}_round2.fastq \
      --quiet --untrimmed-output directly_nop7.fastq -e "${threshold}" "${file}".short.fastq

      Seqs_demult=$(awk 'NR %4 == 2' directly_"${file}"_*_round2.fastq | wc -l)
      tot_demult_anchored=$(expr "${tot_demult_anchored}" + "${Seqs_demult}")
      awk 'NR % 4 == 2 {print length ($0)}' directly_"${file}"_*_round2.fastq > seq_lengths_demulted_anchored_two.step_"${file}".txt

   done
   #Now direct
   for file in  one.step_"${threshold}"_*_round1.fastq; do

      Seqs_deplat=$(awk 'NR %4 == 2' $file | wc -l)
      tot_deplat_direct=$(expr "${tot_deplat_direct}" + "${Seqs_deplat}")
      awk 'NR % 4 == 2 {print length ($0)}' "${file}" > seq_lengths_deplated_anchored_one.step_"${file}"_"${threshold}".txt
      cutadapt -a ATCTCGTATGCCGTCTTCTGCTTG \
               -o  "${file}".short.fastq \
        "${file}"  --discard-untrimmed --quiet -e "${threshold}"



      cutadapt -a file:"${DIRECTION}"/noprimers/all.passing.fastq/A7/barcodes.p7.strict.fasta \
       -o "${file}"_{name}_round2.fastq \
       --quiet --untrimmed-output directly_nop7.fastq -e "${threshold}" "${file}".short.fastq

       Seqs_demult=$(awk 'NR %4 == 2' directly_"${file}"_*_round2.fastq | wc -l)
       tot_demult_direct=$(expr "${tot_demult_direct}" + "${Seqs_demult}")
       awk 'NR % 4 == 2 {print length ($0)}' "${file}"_*_round2.fastq > seq_lengths_demulted_anchored_one.step_"${file}".txt

    done

    printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
    "${threshold}" "${Start}" "${Good_length}" "${Direction_check}"\
    "${Reversed}" "${Anchored}" "${Direct_Anchor}" \
    "${tot_deplat_anchored}" "${tot_deplat_direct}" \
    "${tot_demult_anchored}" "${tot_demult_direct}">> "${OUTPUT_SUMMARY}"

    # Length distributions

    awk 'NR % 4 == 2 {print length ($0)}' "${READ1}".new.fastq > seq_lengths_direction_"${threshold}".txt
    awk 'NR % 4 == 2 {print length ($0)}' "${READ1}".anchored.fastq > seq_lengths_achored_"${threshold}".txt
    awk 'NR % 4 == 2 {print length ($0)}' "${READ1}".anchored.directly.fastq > seq_lengths_achored_directly_"${threshold}".txt

done
