# Using Decona on each well (with all sequences in a run)
# Run this pointing to two folders:
# the output folder from running the demultiplexer
# decona's output folder
# usage of this script 
# Rscript --vanilla from_decona_to_asv.r <demult_folder> <decona_folder>

args=commandArgs(trailingOnly = T)

library(tidyverse)



output.folder <- args[1]
decona.folder <- args[2]

input.files <- list.files(path = decona.folder, pattern = "[0-9]+.fa$",  recursive = T )


cluster.members <- function(list.of.fas, decona.folder){
  
  system2('awk', args = c("'NR % 2 ==1 {print $0, FILENAME}'", file.path(decona.folder, list.of.fas), " >> ", file.path(output.folder, "cluster.members.txt")))
  
}

map(input.files, ~cluster.members(.x, decona.folder))

members.list <- read_table(file.path(output.folder, "cluster.members.txt"), col_names = c("seqid", "cluster"))

members.list %>% 
  mutate(cluster = basename(cluster),
         cluster = str_remove(cluster, ".fa$"))  -> members.list

headers.list  <- file.path(output.folder,"hash.list.csv")

read_csv(headers.list, col_names = c("file", "seqid")) %>% 
  mutate(seqid = str_replace(seqid, "^@", ">"),
         file = basename(file)) %>% 
  inner_join(members.list) -> joindf

joindf %>% 
  group_by(file, cluster) %>% 
  summarise(nReads = n()) -> joindf


## Now bring the consensus seqs

consensus.path <- list.files(decona.folder, pattern = "all_medaka_fastas.fasta", recursive = T, full.names = T)

tibble (seqs = insect::readFASTA(consensus.path, bin = F),
        cluster= names(insect::readFASTA(consensus.path)),
        file = dirname(consensus.path)) %>% 
  mutate(cluster = str_remove_all(cluster, "consensus_medaka-|.fasta")) %>% 
  group_by(seqs) %>% 
  mutate(Hash = map_chr(seqs,digest::sha1,serialize = F)) %>% 
  inner_join(joindf) -> all.together


# Generate the two output files


all.together %>%  
  group_by(file, Hash) %>% 
  summarise(nReads = sum( nReads ))-> ASV_table

write_csv(ASV_table, file.path(output.folder, "ASV_table.csv"))

all.together %>% 
  distinct(seqs, Hash) -> hashes

seqinr::write.fasta(sequences = as.list(hashes$seqs),
                    names     = as.list(hashes$Hash),
                    file.out = file.path(output.folder, "hash_key.fasta"))

hashes %>% 
  write_csv(file.path(output.folder,"hash_key.csv"))


