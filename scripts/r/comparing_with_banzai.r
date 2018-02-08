#comparing outputs between banzai and my pipeline

#Load packages

library(devtools)
library(tidyverse)
library(stringr)
library (dada2)
library (Biostrings)
library(digest)
#Load datasets: DUP_table from ASV and derep.map from banzai (Hashes)

ASV_sequences<-read_csv("/Users/rgallego/fastqs_demultiplexed_for_DADA2/demultiplexed_20180207_1550/ASV_table.csv")

banzai_seqs<-read_delim("/Users/rgallego/fastqs_demultiplexed_for_DADA2/banzai_out_20180207_1530/all_lib/derep.csv",delim =  "\t", col_names = c("Sample","Sequence")) %>%
  separate(Sample,into = c("Sample","nReads"),sep = ";size=") 
banzai_seqs$nReads<-as.integer(banzai_seqs$nReads)

#First search for each ASV seq in the banzai output

#Step 1: do all sequences from the ASV reflect a real DUP found in banzai

ASV_sequences %>% arrange(nReads) %>% select(Sequence)  %>% distinct() -> unique.asvs
banzai_seqs   %>% arrange(nReads) %>% select(Sequence)  %>% distinct() -> unique.banzai.seqs

#How to do that? Let's try first with a inner_join
inner_join(unique.asvs, unique.banzai.seqs, by ="Sequence")
#returns 0 - try self union
test<-unique.asvs
inner_join(unique.asvs, test, by ="Sequence")
#That works - now add unique.asvs to unique banzais
test<-rbind(unique.asvs, unique.banzai.seqs)
inner_join(unique.asvs, test, by ="Sequence")

#So maybe the sequences are reverse complement? because their length is similar

unique.asvs %>%
  mutate (revcom=as.character(reverseComplement(DNAStringSet(Sequence)))) %>% 
  select (-Sequence) -> unique.asvs.revcom
inner_join(unique.asvs.revcom, unique.banzai.seqs, by =c("revcom"="Sequence") )
#YESSSSS 162 / 172 sequences are in both datasets.

ASV_sequences %>%
  mutate (revcom=as.character(reverseComplement(DNAStringSet(Sequence)))) %>%
  anti_join(unique.banzai.seqs, by=c("revcom"="Sequence")) ->missing

missing %>% summarize(sum(nReads))
# 42 out of ASV_sequences %>% summarize(sum(nReads)) 20156 seqs were not found - aka  did not pass the quality filters in banzai

## Now try with hashing
for (i in 1:nrow(unique.asvs.revcom)) {   #for each column of the first dataframe
  
  current_seq<-unique.asvs.revcom[[i,"revcom"]]
  
  hashed_seq<-digest(current_seq, algo= "sha1", serialize = F, skip = "auto")
  
  unique.asvs.revcom[i,"revcom_hash"]<-hashed_seq
}

echo=NULL
pecho=NULL
for (i in 1:nrow(unique.asvs.revcom)){
echo[i]<-str_detect(unique.banzai.seqs[1:100,1], unique.asvs.revcom[[i,1]])
pecho[i]<-str_detect(unique.banzai.seqs[1:100,1], unique.asvs[[i,1]])
}
summary(echo)
summary(pecho)
unique.banzai.seqs[[1,1]]%in% unique.banzai.seqs[1:100,1] %>% summary()

unique.asvs.revcom$revcom %in% unique.banzai.seqs[1:100,1] %>% summary()

for (i in 1:ncol(seqtab.nochim)) {   #for each column of the first dataframe
  
  current_seq<-colnames(seqtab.nochim)[i]
  
  hashed_seq<-digest(current_seq, algo= "sha1", serialize = F, skip = "auto")
  
  colnames(seqtab.nochim)[i]<-hashed_seq
}


ASV_hashes_from_banzai<-read.delim("/Users/rgallego/fastqs_demultiplexed_for_DADA2/banzai_out_20180205_1440/all_lib/derep.map", header = F,sep="\t")
str(ASV_hashes_from_banzai)

ASV_hashes_from_banzai$V1 %>% str_replace ( "SHA1=", "") %>% unique()

colnames(seqtab.nochim) %in% (ASV_hashes_from_banzai$V1 %>% str_replace ( "SHA1=", "") %>% unique())

gather(seqtab.nochim.df, key=Sequence, value = nReads, -sample) %>%
  filter(nReads > 0) 
  
