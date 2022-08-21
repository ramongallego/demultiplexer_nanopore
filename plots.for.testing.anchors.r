library(tidyverse)



read_csv("~/Documents/test_demult/ouput.csv") %>% View()
  pivot_longer(-threshold, names_to = "Step", values_to = "nReads") %>% 
  mutate (Anchoring = case_when(Step %in% c("direction", "anchored", "deplated_direction", "demult_twostep") ~ "Two_step",
                                Step %in% c("direct_anchor", "deplated_direct", "demult_onestep")            ~  "One_step"   ),
          Step = case_when(Step %in% c("anchored", "direct_anchor") ~ "Anchoring",
                           str_detect(Step, "deplated") ~             "Deplating",
                           str_detect(Step, "demult")   ~             "Demultiplexing",
                           TRUE                         ~              Step)) %>% 
  mutate(Step = fct_relevel(Step, "start", "right.length","direction", "Anchoring", "Deplating")) ->temp
  temp %>% add_row(temp %>% filter (is.na(Anchoring)) %>% mutate(Anchoring = "One_step")) %>% 
    replace_na(list(Anchoring= "Two_step")) %>% 
    filter (Step != "reversed") %>% 
  ggplot(aes(x = Step, y = nReads, group = threshold, color = as.character(threshold))) +
  geom_line() +
  facet_wrap(~Anchoring, nrow =1) +
    scale_y_continuous(limits = c(0,NA))


### Each Step success rate
  
  temp %>% 
    filter(!str_detect(Step ,"start|reversed")) -> for.success

for.success %>% 
  add_row(for.success %>% filter (is.na(Anchoring)) %>% replace_na(list(Anchoring = "One_step"))) %>% 
  replace_na(list(Anchoring = "Two_step")) %>% 
  rename(Number= Anchoring ) %>% 
  pivot_wider(names_from = Step, values_from = nReads) %>% 
  group_by(threshold, Number) %>% 
summarise(Direction.finder = 100 * direction/right.length, 
          Anchorage = 100 * Anchoring/right.length,
          Deplating = 100 * Deplating/Anchoring,
          Demultiplexing = 100 * Demultiplexing/Anchoring) %>% 
  ggplot(aes())

  
  list.files(path = "~/Documents/test_demult/testing_anchors/", pattern = "seq_", full.names = T) %>% map(~read_table(.x, col_names = c("length")) %>% 
                                                                             mutate(file = basename(.x))) %>% bind_rows() %>%
  mutate(Anchoring = case_when(str_detect(file, "directly") ~ "One-Step",
                               str_detect(file, "achored") ~ "Two-Step",
                               str_detect(file, "one.step") ~ "One-Step",
                               str_detect(file, "two.step") ~ "Two-Step"),
         Step = case_when(str_detect(file, "start") ~ "Start",
                          str_detect(file, "filter") ~ "Filtered",
                          str_detect(file, "deplated") ~ "First_demult",
                          str_detect(file, "direction") ~ "Direction",
                          str_detect(file, "achored") ~ "Anchored",
                          str_detect(file, "demulted") ~ "Second_demult"),
         Step = fct_relevel(Step, "Start","Filtered", "Direction","Anchored"),
         threshold = str_extract(file, "[:digit:]..."),
         threshold = str_remove(threshold, "[:punct:]$")) -> df

df %>% 
  add_row(
df %>% 
  filter(str_detect(file,("_start|_filter"))) %>% 
  mutate(Anchoring = "One-Step")) %>% 
  replace_na(list(Anchoring = "Two-Step")) %>% 
  filter (!is.na(Anchoring)) %>% 
  ggplot(aes(x = length, after_stat(count))) +
  geom_density(aes(color = threshold), alpha = 0.5) +
  facet_grid(Step~Anchoring, scales = "free_x") +
  scale_x_continuous(limits = c(0,2000))
