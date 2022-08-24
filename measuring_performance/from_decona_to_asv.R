# decona .clstr nfo

input.file <- "~/Documents/test_demult/decona_all_together/all/cluster_representatives.clstr"
input.seqs <- "~/Documents/test_demult/decona_all_together/all/cluster_representatives.fa"
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
    # group_by(Cluster) %>%
    # mutate(Representative = Seq.Name[which(Is.Representative)]) %>%
    # separate_rows(Col2, sep = ",") %>%
    # separate(Col2, into = data.fields, sep = "/", fill = "left", convert = T) %>%
    # mutate(Identity = sub("%", "", Identity) %>% as.numeric) %>%
    # group_by(Seq.Name) %>%
    # mutate(level.rank = paste0(".", 1:n() - 1), level.rank = ifelse(level.rank == ".0", "", level.rank)) %>%
    # pivot_wider(names_from = level.rank, values_from = data.fields, names_sep = "") %>%
    ungroup
}
read.cdhit.clstr(input.file) -> clustered

read_csv("~/Documents/test_demult/decona_all_together/hash.list.csv") %>% 
  mutate(seqid = str_remove(seqid, "^@")) %>%
  left_join(clustered) -> joindf

joindf %>% 
  group_by(file, Cluster, Representative) %>% 
  summarise(nReads = n()) %>% 
  ungroup() %>% 
  select(-Cluster) -> ASV_table

Clusters <- tibble (seqs = insect::readFASTA(input.seqs, bin = F),
                    Representative= names(insect::readFASTA(input.seqs))) %>% 
  group_by(seqs) %>% 
  mutate(Hash = map_chr(seqs, sha1,serialize = F))

ASV_table %>% 
  left_join(Clusters ) %>% 
  select(file, Hash, nReads) -> ASV_table

seqinr::write.fasta(sequences = as.list(Clusters$seqs),
                    names     = as.list(Clusters$Hash),
                    file.out = "~/Documents/test_demult/decona_all_together/hash_key.fasta")

write_csv(ASV_table, "~/Documents/test_demult/decona_all_together/ASV_table.csv")
Clusters %>% 
  select(Hash, seqs) %>% 
  write_csv("~/Documents/test_demult/decona_all_together/hash_key.csv")
