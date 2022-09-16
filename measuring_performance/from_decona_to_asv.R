# decona .clstr nfo

# Prepare this to work with all the decona outputs in a folder
# Because all reads map to the hash.list, and that includes the loci, I thin we can pull 
# this one 

library(tidyverse)



output.folder <- "C:/Users/RG.5015511/Documents/Projects/demultiplexer_nanopore/test_demult/demultiplexed_20220913_1438"

input.file <- list.files(path = output.folder, pattern = "cluster_representatives.clstr", full.names = T, recursive = T )
input.file <- input.file[1]
input.seqs <- list.files(path = output.folder, pattern = "cluster_representatives.fa", full.names = T, recursive = T )
input.seqs <- input.seqs[1]
headers.list  <- file.path(output.folder,"hash.list.csv")
read_table(input.file, col_names = c("Cluster","length" ,"seq", NA,"per") )

read.cdhit.clstr <- function(fname) {
  data.fields <- c("E.Value", "Aln", "Identity")
  read.table(fname, sep = "\t", comment.char = "", quote = "", fill = T, stringsAsFactors = F, col.names = c("Col1", "Col2")) %>%
    separate(Col1, into = c("Seq.Num", "Cluster"), sep = " ", fill = "right") %>%
     fill(Cluster) %>%
     filter(!grepl(">", Seq.Num)) %>%
    separate(Col2, into = c("Seq.Len", "Col2"), sep = "nt, >") %>% # Now something to detect the clusters representative
     separate(Col2, into = c("seqid", "similarity"), sep = "... ") %>%
     mutate(Representative = case_when(similarity == "*" ~ seqid,
                                          TRUE ~ NA_character_),
              similarity = case_when(similarity == "*" ~ "100",
                                   TRUE ~ str_extract(similarity, pattern = "[:digit:]...."))) %>%
    mutate(similarity = as.numeric(similarity)) %>% 
    arrange(Cluster, Representative) %>% 
    fill(Representative) %>% 
    ungroup
}

read.cdhit.clstr(input.file) -> clustered

read_csv(headers.list, col_names = c("file", "seqid")) %>% 
  mutate(seqid = str_remove(seqid, "^@")) %>% 
  inner_join(clustered) -> joindf

joindf %>% 
  group_by(file, Cluster, Representative) %>% 
  summarise(nReads = n()) %>% 
  ungroup() %>% 
  select(-Cluster) -> ASV_table

Clusters <- tibble (seqs = insect::readFASTA(input.seqs, bin = F),
                    Representative= names(insect::readFASTA(input.seqs))) %>% 
  group_by(seqs) %>% 
  mutate(Hash = map_chr(seqs,digest::sha1,serialize = F))

ASV_table %>% 
  left_join(Clusters ) %>% 
  select(file, Hash, nReads) -> ASV_table

ASV_table %>% 
  mutate(file =basename(file))-> ASV_table

seqinr::write.fasta(sequences = as.list(Clusters$seqs),
                    names     = as.list(Clusters$Hash),
                    file.out = file.path(output.folder, "hash_key.fasta"))

write_csv(ASV_table, file.path(output.folder, "ASV_table.csv"))
Clusters %>% 
  select(Hash, seqs) %>% 
  write_csv(file.path(output.folder,"hash_key.csv"))


