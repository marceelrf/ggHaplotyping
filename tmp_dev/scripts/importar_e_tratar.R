
# Importar e tratar a base SABE -------------------------------------------

# Dia 18/09/2024

library(tidyverse)
library(glue)

path = "tmp_dev/data/SABE_haplotype_join.xlsx"

sheets <-
  readxl::excel_sheets(path = path)

data_list <-
  sheets %>%
  map(.f = \(x) readxl::read_xlsx(path, sheet = x))

names(data_list) <- sheets


data_list %>%
  map(\(x) x %>% mutate(across(where(is.character),
                               ~ na_if(., "none"))))
