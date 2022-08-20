library(tidyverse)



read_csv("~/Documents/test_demult/ouput.csv") %>% 
  pivot_longer(-threshold, names_to = "Step", values_to = "nReads") %>% 
  mutate(Step = fct_relevel(Step, "start", "right.length","direction", "anchored", "direct_anchor")) %>% 
  ggplot(aes(x = Step, y = nReads, group = threshold, color = as.character(threshold))) +
  geom_line()



list.files("~/Projects/demultiplexer_nanopore/", pattern = "seq_") %>% map(~read_table(.x, col_names = c("length")) %>% 
                                                                             mutate(file = .x)) %>% bind_rows() %>% 
  mutate(Anchoring = case_when(str_detect(file, "directly") ~ "One-Step",
                               str_detect(file, "achored") ~ "Two-Step",
                               str_detect(file, "direction") ~ "Two-Step"),
         Step = case_when(str_detect(file, "start") ~ "Start",
                          str_detect(file, "filter") ~ "Filtered",
                          str_detect(file, "deplated") ~ "First_demult",
                          str_detect(file, "direction") ~ "Direction",
                          str_detect(file, "achored") ~ "Anchored"),
         Step = fct_relevel(Step, "Start","Filtered", "Direction","Anchored"),
         threshold = str_extract(file, "[:digit:]..."),
         threshold = str_remove(threshold, "[:punct:]$")) -> df
df %>% 
  filter(Step == "Filtered")


df %>% 
  filter (!is.na(Anchoring)) %>% 
  ggplot(aes(x = length, after_stat(count))) +
  geom_density(aes(color = threshold), alpha = 0.5) +
  facet_grid(Step~Anchoring, scales = "free_x")
