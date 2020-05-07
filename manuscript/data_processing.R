manuscript_data <- function(){
  require(PFAScreeneR)
  require(cfmR)
  require(tidyverse)
  require(MSnbase)

  pfas_exp_msms <-
    readMgfData(list.files(
      system.file(package = 'PFAScreeneR'),
      pattern = '.mgf',
      full.names = T
    ))

  db_file <- list.files(system.file(package = 'PFAScreeneR'), pattern = '[.]db', full.names = T)
  pol <- ifelse(fData(pfas_exp_msms)$CHARGE == '1-', 'negative', 'positive')
  spec_table_name <- ifelse(pol == 'negative', 'cfm_neg', 'cfm_pos')
  ID <- as.character(fData(pfas_exp_msms)$ID)
  instrument_type <- as.character(fData(pfas_exp_msms)$INSTRUMENT)
  ppm <- ifelse(instrument_type == 'Orbitrap', 10, 50)
  mol_form <- as.character(fData(pfas_exp_msms)$FORMULA)
  exact_mass <- as.numeric(as.character(fData(pfas_exp_msms)$EXACT_MASS))
  nom_mass_ppm <- sapply(exact_mass, function(x) 1e6*(diff(x + c(-0.5,0.5))/x))

  # search by ID
  sim_id_search <- pbapply::pblapply(seq(along = pfas_exp_msms),
                                     function(x) {
                                       cfm_search(
                                         query_spectrum = pfas_exp_msms[[x]],
                                         spec_table_name = spec_table_name[x],
                                         pol = pol[x],
                                         db_file = db_file,
                                         mol_table_name = 'molecules',
                                         molprop_table_name = 'mol_props',
                                         fun = 'dotproduct',
                                         mol_form = NULL,
                                         exact_mass = NULL,
                                         bin = 0.005,
                                         ID = ID[x],
                                         ppm = ppm[x]
                                       )
                                     })

  # search by mol form
  sim_molform_search <- pbapply::pblapply(seq(along = pfas_exp_msms),
                                          function(x) {
                                            cfm_search(
                                              query_spectrum = pfas_exp_msms[[x]],
                                              spec_table_name = spec_table_name[x],
                                              pol = pol[x],
                                              db_file = db_file,
                                              mol_table_name = 'molecules',
                                              molprop_table_name = 'mol_props',
                                              fun = 'dotproduct',
                                              mol_form = mol_form[x],
                                              exact_mass = NULL,
                                              bin = 0.005,
                                              ID = NULL,
                                              ppm = ppm[x]
                                            )
                                          })


  # search by neutral exact mass
  sim_neutralmass_search <- pbapply::pblapply(seq(along = pfas_exp_msms),
                                              function(x) {
                                                cfm_search(
                                                  query_spectrum = pfas_exp_msms[[x]],
                                                  spec_table_name = spec_table_name[x],
                                                  pol = pol[x],
                                                  db_file = db_file,
                                                  mol_table_name = 'molecules',
                                                  molprop_table_name = 'mol_props',
                                                  fun = 'dotproduct',
                                                  mol_form = NULL,
                                                  exact_mass = exact_mass[x],
                                                  bin = 0.005,
                                                  ID = NULL,
                                                  ppm = ppm[x]
                                                )
                                              })

  # search by neutral nominal mass
  sim_nominalneutralmass_search <- pbapply::pblapply(seq(along = pfas_exp_msms),
                                                     function(x) {
                                                       cfm_search(
                                                         query_spectrum = pfas_exp_msms[[x]],
                                                         spec_table_name = spec_table_name[x],
                                                         pol = pol[x],
                                                         db_file = db_file,
                                                         mol_table_name = 'molecules',
                                                         molprop_table_name = 'mol_props',
                                                         fun = 'dotproduct',
                                                         mol_form = NULL,
                                                         exact_mass = exact_mass[x],
                                                         bin = 0.005,
                                                         ID = NULL,
                                                         ppm = nom_mass_ppm[x]
                                                       )
                                                     })

  # search by exact precursor m/z
  sim_precurmz_search <- pbapply::pblapply(seq(along = pfas_exp_msms),
                                           function(x) {
                                             cfm_search(
                                               query_spectrum = pfas_exp_msms[[x]],
                                               spec_table_name = spec_table_name[x],
                                               pol = pol[x],
                                               db_file = db_file,
                                               mol_table_name = 'molecules',
                                               molprop_table_name = 'mol_props',
                                               fun = 'dotproduct',
                                               mol_form = NULL,
                                               exact_mass = NULL,
                                               bin = 0.005,
                                               ID = NULL,
                                               ppm = ppm[x]
                                             )
                                           })

  # search by nominal precursor m/z
  sim_nominalprecurmz_search <- pbapply::pblapply(seq(along = pfas_exp_msms),
                                                  function(x) {
                                                    cfm_search(
                                                      query_spectrum = pfas_exp_msms[[x]],
                                                      spec_table_name = spec_table_name[x],
                                                      pol = pol[x],
                                                      db_file = db_file,
                                                      mol_table_name = 'molecules',
                                                      molprop_table_name = 'mol_props',
                                                      fun = 'dotproduct',
                                                      mol_form = NULL,
                                                      exact_mass = NULL,
                                                      bin = 0.005,
                                                      ID = NULL,
                                                      ppm = nom_mass_ppm[x]
                                                    )
                                                  })
  search_rst <- list(
    `Identifier` = sim_id_search,
    `Mol. form.` = sim_molform_search,
    `Accurate neutral mass` = sim_neutralmass_search,
    `Accurate precur. m/z` = sim_precurmz_search,
    `Nominal neutral mass` = sim_nominalneutralmass_search,
    `Nominal precur. m/z` = sim_nominalprecurmz_search
  ) %>%
    lapply(setNames, nm = ID)
  return(search_rst)
}






