# load original lengths


output_folder <- "~/Documents/test_demult/demultiplexed_20220907_1513/"

# original.lengths <- read_table(file.path(fastq.folder, "seq.lengths.init"),
                               # col_names = c("seq","length.init"))

# Load legnths after redirecting fwd

found.right.primer <- read_table(file.path (fastq.folder, "seq.length.reorientation.primer.found"), col_names = c("seq", "length")) %>% 
  mutate(Step = "2.Found.Fwd.primer.and.removed")

not.found.right.primer <- read_table(file.path (fastq.folder, "seq.length.reorientation.primer.not.found"), col_names = c("seq", "length")) %>% 
  mutate(Step = "2.Fwd.primer.not.found")

# Load lengths after finding the reverse primer

found.right.andleft.primer <- read_table(file.path(fastq.folder, "seq.length.reorientation.both.primers"), col_names = c("seq", "length")) %>% 
  mutate(Step = "3.Found.both.primers.and.removed")

found.right.not.left.primer <- read_table(file.path(fastq.folder, "seq.length.reorientation.only.fwd"), col_names = c("seq", "length")) %>% 
  mutate(Step = "3.Found.Fwd.primer.no.rev")

## Do some checks

First_direction.check <- bind_rows(found.right.primer, not.found.right.primer) %>% left_join(original.lengths)

First_direction.check %>% 
  filter(length.init < 2000) %>% 
  ggplot(aes(y = length, x = Step, fill = Step)) + 
  geom_violindot(fill_dots = "black") +
  theme_modern()+
  scale_fill_material_d()

Second_direction.check <- bind_rows(found.right.andleft.primer, found.right.not.left.primer, not.found.right.primer) %>% left_join(original.lengths)
Second_direction.check %>% 
  filter(length.init < 2000) %>% 
  ggplot(aes(y = length, x = Step, fill = Step)) + 
  geom_violindot(fill_dots = "black") +
  theme_modern()+
  scale_fill_material_d() 


# All in one plot

keepers <- original.lengths %>% filter(length.init < 2000) %>% pull(seq)

good.length <- original.lengths %>% filter (length.init < 2000 & length.init > 1000) %>% pull(seq)

original.lengths %>% 
  mutate(Step = "1.Starting") %>%
  rename(length = length.init) %>% 
  bind_rows(found.right.andleft.primer, found.right.not.left.primer, not.found.right.primer, found.right.primer) %>% 
  filter(seq %in% keepers) -> ready.for.fig1 
  

ready.for.fig1 %>%  
ggplot(aes(y = length, x = Step, fill = Step)) + 
  geom_violindot(fill_dots = "black") +
  theme_modern()+
  scale_fill_material_d() 

ready.for.fig1 %>% 
  mutate(good.length = seq %in% good.length) %>% 
  ggplot(aes(y = length, x = Step, fill = Step)) + 
  geom_violindot(fill_dots = "black") +
  theme_modern()+
  scale_fill_material_d() +
  facet_wrap(~good.length, nrow = 2) +
  guides(fill = "none") 

ggsave(filename = file.path(fastq.folder, "sequence.length.distribution.png"), dpi = "retina", width = 16)  

####### My take from this is that because cutadapt looks 

found.right.andleft.primer.after.rc <- read_table(file.path(fastq.folder, "seq.lengths.rev.looking.from.left"), col_names = c("seq", "length")) %>% 
  mutate(Step = "3.Found.both.primers.and.removed.after.rc") %>% filter(seq %in% keepers)
ready.for.fig1 %>% 
  bind_rows(found.right.andleft.primer.after.rc) %>% 
  
  filter( seq %in% good.length) %>% 
  ggplot(aes(y = length, x = Step, fill = Step)) + 
  geom_violindot(fill_dots = "black") +
  theme_modern()+
  scale_fill_material_d() +
  #facet_wrap(~good.length, nrow = 2) +
  guides(fill = "none")


# Bring in the lengths after each process
output_folder <- "~/Documents/test_demult/demultiplexed_20220909_0947/"

demult.lengths.by.plate <- read_table(file.path(output_folder, "seq.lengths.demult.by.plate.txt"), col_names = c("seq","length.plate", "file")) %>% 
  mutate(file = str_remove(file, "_round1.fastq"),
         file = basename(file))

demult.lengths.adapters <- read_table(file.path(output_folder, "seq.lengths.adapters.txt"), col_names = c("seq","length.plate", "file")) %>% 
  mutate(file = str_remove(file, "_adap.fastq"),
         file = basename(file))

demult.lengths.by.well <- read_table(file.path(output_folder, "seq.lengths.demult.by.well"), col_names = c("seq","length.well", "file"))%>% mutate(file = str_remove(basename(file), ".fastq"))

demult.lengths.primer.removed <- read_table(file.path(output_folder, "seq.lengths.demult.primer.removed"), col_names = c("seq","length.amplicon", "file")) %>% 
  mutate(file = basename(file),
         file = str_remove(file, "_16S_long.fastq")) 

demult.lengths.by.plate %>% 
  mutate(version = "Universal") -> demult.lengths.by.plate.universal

demult.lengths.by.plate <- read_table(file.path(output_folder, "seq.lengths.demult.by.plate.txt"), col_names = c("seq","length.plate", "file")) %>% 
  mutate(file = str_remove(file, "_round1.fastq"),
         file = basename(file)) %>% 
  mutate(version = "Anchored")

bind_rows(demult.lengths.by.plate.universal, demult.lengths.by.plate) 

demult.lengths.by.plate.universal %>%   
  ggplot ( aes(y = length.plate, x = file, fill = file)) +
  geom_violindot(fill_dots = "black") +
  theme_modern()+
  scale_fill_material_d() +
  facet_wrap(~version, nrow = 2) +
  guides(fill = "none")

bind_rows(demult.lengths.by.plate.universal, demult.lengths.by.plate) %>% 
ggplot( aes(x = file, fill= length.plate>1000)) +
  geom_bar(position = "stack")+
  facet_wrap(~version, nrow = 2)

demult.adapters <- read_table(file.path(output_folder, "seq.adapters.txt"), col_names = c("seqID","seq", "length.plate", "file")) %>% 
  mutate(file = str_remove(file, "_round1.fastq"),
         file = basename(file))

demult.adapters %>% 
  group_by(seq) %>% 
  tally(sort = T)

demult.adapters %>% 
  filter (seq == "TCTACACCTA")

output_folder <- "~/Documents/test_demult/demultiplexed_20220909_1006/"
after.pcr.primer <- read_table(file.path(output_folder, "seq.lengths.after.rev.primer.txt"), col_names = c("seq", "length.plate", "file")) %>% 
  mutate(file = str_remove(file, ".anchored.rc.fastq"),
         file = basename(file))
after.pcr.primer %>% 
  ggplot ( aes(y = length.plate, x = file, fill = file)) +
  geom_violindot(fill_dots = "black") +
  theme_modern()+
  scale_fill_material_d() +
  guides(fill = "none")

final.length <- read_table(file.path(output_folder, "seq.lengths.demult.by.well.txt"), col_names = c("seq", "length.plate", "file")) %>% 
  mutate(file = basename(file)) %>% 
  separate(file, into = c("Plate", "Well"), sep = "_Well_")
final.length %>% 
  ggplot ( aes(y = length.plate, x = Plate, fill = Plate)) +
  geom_violindot(fill_dots = "black") +
  theme_modern()+
  scale_fill_material_d() +
  guides(fill = "none")
