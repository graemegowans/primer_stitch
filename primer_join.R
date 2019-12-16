#get 60bp overlap primers for making constructs

#load packages
library(tidyverse)

#use dummy data####

#list of parts
lkup <- c("gg1", "gg2", "gg3", "gg4", "gg5")
df <- NULL

#generate random sequences for parts
for (i in 1:length(lkup)) {
 
  df[[i]] <- tibble(
  seq = paste0(sample(c("A", "T", "G", "C"), 
                      200, replace = TRUE), 
               collapse = ""))
}

#make dataframe
df <- map_df(df, bind_rows)

#join with names for final dummy table
lkup <- tibble(name = lkup, seq = df$seq)

#test part_names
test <- tibble(part_name = c("gg1_gg2_gg3", 
                             "gg2_gg3_gg4", 
                             "gg1_gg4", 
                             "gg2_gg3_gg4_gg5",
                             "gg1_gg2_gg3_gg5"))

#or if you want to load own data####

#read parts lists
test <- read_csv("your_parts_list.csv")

#if format is in separate columns
#change a, b, c to the name of the columns
test <- unite(test, part_name, colnames(test)) %>% 
  mutate(part_name = str_replace_all(part_name, "_NA", "")) %>% 
  mutate(part_name = str_replace_all(part_name, "NA_", ""))

#read lookup table of primers
lkup <- read_csv("your_lookup_table.csv")

#start here for generating primers####

#make empty list to store all primers
test_list <- list()

#for every entry in table
for (j in 1:nrow(test)) {
  
  #extract part name  
  xy <- test$part_name[j]
  
  #get parts and make vector
  parts <- xy %>% str_split("_")
  parts <- as.vector(parts[[1]])
  
  #generate primers - one per part
  #make empty list to store
  primers <- list()
  
  for (i in 1:(length(parts)-1)) {
    
    #get the sequences for matching
    a <- tibble(name = parts[i]) %>% left_join(lkup) %>% select(seq) %>% paste0()
    b <- tibble(name = parts[i+1]) %>% left_join(lkup) %>% select(seq) %>% paste0()
    
    #generate 30 bp sequences and paste together
    a_seq <- str_sub(a, start = nchar(a)-29, end = nchar(a))
    b_seq <- str_sub(b, start = 1, end = 30)
    
    #add to table with part name, primer name and sequence
    primers[[i]] <- tibble(
      part_name = xy,
      primer_name = paste0(parts[i], "_", parts[i+1]),
      seq = paste0(a_seq, b_seq))
  }
  
  #flatten primers
  primers <- bind_rows(primers)
  
  #join to input and add to list
  test_list[[j]] <- inner_join(test, primers, by = "part_name")
  
}

#make final table, marking duplicates
primer_list <- bind_rows(test_list) %>% mutate(dup = duplicated(seq))

#show data
primer_list

#export
write_csv("primer_output.csv")


