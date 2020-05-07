# SI from https://pubs.acs.org/doi/10.1021/acs.est.8b06017
library(readxl)
library(tidyverse)
library(cfmR)

# read features
xl_path <- '~/Desktop/McCord_EST_SI/Supplemental Table 1_All Chemical Features.xlsx'
sheets <- excel_sheets(xl_path)
features <- map_dfr(sheets, read_excel, path = xl_path, .id = 'sample_location') %>%
  mutate(sample_location = factor(sample_location, levels = c(1,2), labels = sheets))

# read IDs
ids <- read_excel('~/Desktop/McCord_EST_SI/Supplemental Table 2_Assigned Chemical Species.xlsx')

# # merge data
# dat <- full_join(
#   features,
#   mutate(ids, Cpd = as.numeric(str_extract(`Compound ID`, '[1-9]{1,2}')))
#   ) %>%
#   distinct(sample_location, Cpd, Mass, Vol) %>%
#   pivot_wider(id_cols = c(Cpd, Mass), names_from = sample_location, values_from = Vol)

# pivot_wider(
#   features,
#   id_cols = c(Cpd, Mass),
#   names_from = sample_location,
#   values_from =  Vol
# )


# search database by mass
db_file <- '~/Desktop/PFAScreeneR_data/PFAScreeneR.db'
## read mols
mols <-
  cfm_read_db(db_file = db_file,
              table_name = "molecules",
              ID = NULL)
mol_props <-
  cfm_read_db(db_file = db_file,
              table_name = "mol_props",
              ID = NULL)

ppm <- 10
search_rst <- plyr::adply(
  distinct(features, Cpd, Mass),.margins = 1,
  .f = function(d) {
    ## calculate ppm error and filter
    ppm_error <- purrr::map_dbl(
      mol_props$ExactMass,
      .f = function(x) {
        abs(1e6 * ((d[['Mass']] - x) / x))
      }
    )
    mutate(mol_props, ppm_error = ppm_error) %>%
      filter(ppm_error < (ppm / 2))
  }
)

search_summary <- search_rst %>%
  as_tibble() %>%
  group_by(Cpd) %>%
  summarise_at(.vars = vars(ID,MolForm), .funs = list(~length(unique(.)))) %>%
  left_join(x = features) %>%
  select(sample_location, Cpd, ID, MolForm, Vol) %>%
  group_by(sample_location)

search_summary %>%
  ungroup %>%
  arrange(Cpd, ID) %>%
  distinct(Cpd,ID) %>%
  summarise_at(vars(Cpd, ID), list(~sum(!is.na(.))))


summarise_at(search_summary, vars(Vol, ID), .funs = list(~sum(!is.na(.)))) %>%
  rename(`n detected` = Vol, `n matched` = ID) %>%
  mutate(`percent matched` = 100*(`n matched`/`n detected`))


d <- distinct(dat, Cpd, Vol, .keep_all = T) %>%
  arrange(desc(Vol)) %>%
  mutate(fill = ifelse(is.na(`Compound ID`), 'white', 'black')) %>%
  select(Cpd, Vol, fill)

barplot(d$Vol, col = d$fill)
