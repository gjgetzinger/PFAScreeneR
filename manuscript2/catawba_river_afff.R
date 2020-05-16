library(network)
library(tidyverse)
library(DBI)
library(ggraph)
library(ggplot2)
library(reticulate)
library(rjson)
library(tidygraph)
library(plyr)
library(reshape2)
library(readr)
library(cfmR)
library(cdutils)

source_python('manuscript2/struc_similarity.py')
source('manuscript2/plotting_functions.R')

# connect to CD result file and PFAS database
cd_rst_file <- '/Users/gordong/Desktop/catawba_river_afff_eg/Catawba River AFFF-(22).cdResult'
cd <- dbConnect(RSQLite::SQLite(), cd_rst_file)
pfas_conn <- dbConnect(RSQLite::SQLite(), '~/Desktop/PFAScreeneR_data/PFAScreeneR.db')

# read consolidated unknown compound items
cpds <- tbl(cd, 'ConsolidatedUnknownCompoundItems') %>%
  as_tibble() %>%
  filter(Checked == 1) %>%
  as_tibble() %>%
  left_join(
    y = tbl(pfas_conn, 'dsstox_mapping') %>%
      as_tibble() %>%
      mutate(Name = substr(InChIKey, 1, 14))
  )

# write cpds with dtx ids to file
filter_at(cpds, vars(contains('dsstox')), any_vars(!is.na(.))) %>%
  group_by(ID) %>%
  select(contains('dsstox'), 'casrn')  %>%
  write_csv(path = '~/Desktop/PFAScreeneR/manuscript2/catawba_river_dsstox_casrn_hits.csv')

# collect pfas data and assign some variables
cf2 <- 49.99681
dat <- cpds %>%
  dplyr::rename(
    CD_ID = ID,
    ID = Name,
    InChIKey_dsstox = InChIKey
  ) %>%
  left_join(
  y = tbl(pfas_conn, 'molecules') %>%
    select(ID, InChIKey, InChI) %>%
    as_tibble()) %>%
  left_join(
    y = tbl(pfas_conn, 'mol_props') %>% select(ID, MolForm) %>% as_tibble()
  ) %>%
  rowwise() %>%
  mutate(
    F_count = as.integer(cfmR::clean_form(MolForm, as_cts = T)[['F']]),
    H_count = as.integer(cfmR::clean_form(MolForm, as_cts = T)[['H']]),
    O_count = ifelse(str_detect(MolForm, 'O'), as.integer(cfmR::clean_form(MolForm, as_cts = T)[['O']]), 0),
    C_count = as.integer(cfmR::clean_form(MolForm, as_cts = T)[['C']])
  ) %>%
  ungroup() %>%
  mutate(
    HC = H_count/C_count,
    OC = O_count/C_count,
    `Kendrick mass (CF2)` = MolecularWeight * (floor(cf2)/cf2),
  ) %>%
  mutate(
     `Nominal Kendrick mass (CF2)` = floor(`Kendrick mass (CF2)`),
  ) %>%
  mutate(
    `Kendrick mass defect (CF2)` = `Kendrick mass (CF2)` - `Nominal Kendrick mass (CF2)`
   )

# For each detected feature, make a network and check for connected nodes, reduce
# list of networks into one network and assign node attributes
rxn_graph <- alply(dat, 1, function(x) {
  # look for biotransformation products/parents
  bio <- tbl(pfas_conn, 'bio_rxns') %>%
    filter_at(vars(contains('Precursor_InChIKey')), .vars_predicate = any_vars(. %in% !!x$InChIKey)) %>%
    select(contains("InChIKey")) %>%
    as_tibble() %>%
    dplyr::rename(Product_InChIKey = InChIKey)
  if (nrow(bio) > 0) {
    bio <- mutate(bio, rxn_type = 'bio') %>% distinct_all()
  } else {
    bio <- mutate(bio, rxn_type = NULL)
  }

  # look for hydrolysis products/parents
  hyd <- tbl(pfas_conn, 'hyd_rxns') %>%
    filter_at(vars(contains('Precursor_InChIKey')), .vars_predicate = any_vars(. %in% !!x$InChIKey)) %>%
    select(contains("InChIKey")) %>%
    as_tibble() %>%
    distinct_all()
  if (nrow(hyd) > 0) {
    hyd <- mutate(hyd, rxn_type = 'hyd') %>% distinct_all()
  } else {
    hyd <- mutate(hyd, rxn_type = NULL)
  }

  # combine transformations and summarise
  rxn <- bind_rows(bio, hyd)
  if (nrow(rxn) > 0) {
    rxn <- group_by_at(rxn, vars(ends_with("InChIKey"))) %>%
      arrange(rxn_type) %>%
      summarize_at(vars(rxn_type), list(~ paste(., collapse = '|'))) %>%
      ungroup() %>%
      distinct_all() %>%
      filter_at(vars(ends_with('InChIKey')), all_vars(. %in% !!dat$InChIKey))
    if (nrow(rxn) > 0) {
      mutate(rxn, CD_ID = x$CD_ID) %>%
        as_tbl_graph()
    } else {
      transmute(
        x,
        Precursor_InChIKey = InChIKey,
        Product_InChIKey = InChIKey,
        rxn_type = NA
      ) %>%
        as_tbl_graph()
    }
  } else {
    transmute(
      x,
      Precursor_InChIKey = InChIKey,
      Product_InChIKey = InChIKey,
      rxn_type = NA
    ) %>%
      as_tbl_graph()
  }
}) %>%
  Reduce('graph_join', .) %>%
  set_vertex_attr(.,
                  name = 'Max area',
                  value = dat[match(get.vertex.attribute(rxn_graph)$name, dat$InChIKey),] %>%
                    pull(MaxArea)
  ) %>%
  set_vertex_attr(.,
                  name = 'retention_time',
                  value = dat[match(get.vertex.attribute(rxn_graph)$name, dat$InChIKey),] %>%
                    pull(RetentionTime)
  ) %>%
  set_vertex_attr(.,
                  name = 'mol_wt',
                  value = dat[match(get.vertex.attribute(rxn_graph)$name, dat$InChIKey),] %>%
                    pull(MolecularWeight)
  ) %>%
  set_vertex_attr(.,
                  name = 'label',
                  value = dat[match(get.vertex.attribute(rxn_graph)$name, dat$InChIKey),] %>%
                    pull(CD_ID)
  ) %>%
  set_vertex_attr(
    .,
    name = 'Fluorine count',
    value = dat[match(get.vertex.attribute(rxn_graph)$name, dat$InChIKey), ] %>%
      pull(F_count)
    ) %>%
  set_vertex_attr(
    .,
    name = 'H/C',
    value = dat[match(get.vertex.attribute(rxn_graph)$name, dat$InChIKey), ] %>%
      pull(HC)
  ) %>%
  set_vertex_attr(
    .,
    name = 'O/C',
    value = dat[match(get.vertex.attribute(rxn_graph)$name, dat$InChIKey), ] %>%
      pull(OC)
  ) %>% set_vertex_attr(
    .,
    name = 'Kendrick mass (CF2)',
    value = dat[match(get.vertex.attribute(rxn_graph)$name, dat$InChIKey), ] %>%
      pull(`Kendrick mass (CF2)`)
  ) %>%
  set_vertex_attr(
    .,
    name = 'Kendrick mass defect (CF2)',
    value = dat[match(get.vertex.attribute(rxn_graph)$name, dat$InChIKey), ] %>%
      pull(`Kendrick mass defect (CF2)`)
  )

# plot network diagrams/graphs for several variables
mz_rt <- plot_graph(rxn_graph, 'retention_time', 'mol_wt', 'Retention Time (min)', 'Molecular Weight (Da)', "identity")
van_krev <- plot_graph(rxn_graph,  "O/C", "H/C",  "O/C", "H/C", "identity")
mass_def <- plot_graph(rxn_graph,  "Kendrick mass (CF2)", "Kendrick mass defect (CF2)",  "Kendrick mass (CF2)", "Kendrick mass defect (CF2)", "identity")

# annotate results with SIRIUS
ms_dir <- 'manuscript2/ms_files/'
dir.create(ms_dir)
map(spec, function(x,y) {
 cat(
    paste(">compound", x$ID),
    paste(">formula", x$featureData$ElementalCompositionFormula),
    paste(">ionization", x$adduct),
    paste(">ms1"),
    matrix(apply(x$ms1, 1,paste, collapse = " "),ncol = 1),
    paste(">ms2"),
    matrix(apply(x$ms2,1,paste,collapse = " "), ncol = 1)
    , sep = '\n', file = paste0(ms_dir, x$ID, '.ms'))
})

#scp -r manuscript2/ms_files/ gjg3@dcc-slogin.oit.duke.edu:~
dir.create('~/sirius_out')
sapply(
  list.files(ms_dir, full.names = T),
  function(x){
    cmd <- sprintf("sirius --ppm-max 2.5 -p orbitrap --output sirius_out/ %s", x)
    system(cmd)
  }
)

#scp -r gjg3@dcc-slogin.oit.duke.edu:~/sirus_out /Users/gordong/Desktop/PFAScreeneR/manuscript2
sirius_rst <- map_dfr(
  list.dirs('manuscript2/sirus_out/', full.names = T, recursive = F),
  function(x){
    CD_ID <- str_extract(x, '(?<=_)[0-9]{1,5}$')
    rst <- read_tsv(list.files(x, '.csv', full.names = T))
    tree <- fromJSON(file = list.files(paste0(x, '/trees'), '.json', full.names = T))
    frags <-  map_dfr(tree$fragments, function(x){as_tibble(x[c('id', 'molecularFormula')])}) %>%
      mutate(CD_ID = CD_ID) %>%
      group_by(CD_ID) %>%
      nest(tree = c(id, molecularFormula))
    mutate(rst,
      CD_ID
    ) %>%
      left_join(frags)
  })


# make a network of SIRIUS annotated fragment ions
mf_graphs <- alply(sirius_rst, 1, function(x) {
  x$tree[[1]] %>%
    mutate(pre_mf = x$formula) %>%
    mutate(pro_mf = molecularFormula) %>%
    mutate(CD_ID = x$CD_ID) %>%
    select(pre_mf, pro_mf, everything()) %>%
    as_tbl_graph() %>%
    set_vertex_attr(.,
                    'Ion type',
                    value = ifelse(pull(., name) %in% x$formula,
                                   'precursor', 'product')) %>%
    set_vertex_attr(.,
                    'label',
                    value = ifelse(pull(., name) %in% x$formula,
                                   x$CD_ID, ""))
}) %>% Reduce("graph_join", .)

f_ct <- sapply(as_tibble(mf_graphs)$name, function(x){
  ifelse(grepl('F', x), cfmR::clean_form(x, as_cts = T)[['F']], 0)
}, USE.NAMES = F)

mf_graphs <- set_vertex_attr(mf_graphs, 'f_ct', value = f_ct)
mf_graphs <- set_edge_attr(mf_graphs, 'rxn_type', value = rep('MS/MS'))

ly <- layout_nicely(mf_graphs)
sirius_graph <- ggraph(mf_graphs, layout = ly) +
  geom_edge_link0(edge_color = 'darkgray', alpha = 0.5) +
  geom_node_point(aes(color = factor(f_ct)), alpha = 0.75) +
  F_colScale +
  scale_shape(guide_legend(title = "Ion type")) +
  geom_node_text(aes(label = label), check_overlap = T, repel = T, color = 'black', size = 2) +
  theme_graph(background = "white") +
  theme(legend.position = 'bottom') +
  labs(subtitle = 'MS/MS fragment trees')

# combine plots to make figure 2
figure2 <- ggarrange(
  mz_rt,
  van_krev,
  mass_def,
  sirius_graph,
  nrow = 2, ncol = 2,
  common.legend = T,
  labels = "AUTO",
  legend = "bottom"
)

ggsave(
  'manuscript2/figures/figure2_catawba_riv_unkns.png',
  figure2,
  width = 7,
  height = 5,
  units = 'in',
  dpi = 300,
  scale = 1.5
)


# determine connections in rxn graph and ms/ms graph
rst <- lapply(dat$CD_ID, function(CD_ID){
  v1 <- rxn_graph %>%
    activate("nodes") %>%
    mutate(id = row_number()) %>%
    filter(label == CD_ID) %>%
    pull(id)

  if(length(v1) > 0){
    to1 <- rxn_graph %>%
      activate("nodes") %>%
      mutate(id = row_number()) %>%
      as_tibble() %>%
      pull(id)

    node_dist1 <- distances(rxn_graph, V(rxn_graph)[v1], V(rxn_graph)[to1], mode = "all")
    inchikeys <- names(node_dist1[1,node_dist1 <= 4])
  } else {
    inchikeys <- NULL
  }

  v2 <- mf_graphs %>%
    activate("nodes") %>%
    mutate(id = row_number()) %>%
    filter(label == CD_ID) %>%
    as_tibble() %>%
    pull(id)

  if(length(v2) > 0){
    to2 <- mf_graphs %>%
      activate("nodes") %>%
      mutate(id = row_number()) %>%
      filter(`Ion type` == 'precursor') %>%
      as_tibble() %>%
      pull(id)

    node_dist2 <- distances(mf_graphs, V(mf_graphs)[v2], V(mf_graphs)[to2], mode = "all")
    node_dist2 <- names(node_dist2[1,node_dist2 <= 4])
    cd_ids <- mf_graphs %>%
      activate("nodes") %>%
      filter(name %in% node_dist2, `Ion type` == 'precursor', label != "")  %>%
      pull(label) %>%
      as.integer()
  } else {
    cd_ids <- NULL
  }

  rst <- dat %>%
    filter(InChIKey %in% inchikeys | CD_ID %in% cd_ids) %>%
    mutate(
      `Mass error (ppm)` = 1e6*((exact_mass - MolecularWeight)/exact_mass)
    ) %>%
    select(CD_ID,
           ID,
           Name = preferred_name,
           `Molecular weight (Da)` = MolecularWeight,
           `Retention time (min)` = RetentionTime,
           `Mol. form.` = MolForm,
           InChI)
  return(rst)
})

rst <- setNames(rst, dat$CD_ID)

# export an annotated table/scheme for CD_ID=3
frame <- within(rst$`3`, {
  smiles = sapply(sapply(InChI, rdkit$Chem$MolFromInchi), rdkit$Chem$MolToSmiles)
  exact_mass = sapply(sapply(InChI, rdkit$Chem$MolFromInchi), rdkit$Chem$rdMolDescriptors$CalcExactMolWt)
  }) %>%
  mutate(`Mass error (ppm)` = 1e6*((exact_mass - `Molecular weight (Da)`)/exact_mass)) %>%
  select(-InChI, -exact_mass)

# keep cpds with 6:2 fluorotelomer backbone
patt <- rdkit$Chem$MolFromSmarts('[S]CCC(F)(C(F)(C(F)(C(F)(C(F)(C(F)(F)F)F)F)F)F)F')
mols <- sapply(frame$smiles, rdkit$Chem$MolFromSmiles)
has_patt <- sapply(mols, function(x){x$HasSubstructMatch(patt)})
frame <- filter(frame, has_patt) %>%
  distinct(ID, .keep_all = T)

# determine R-groups and write a csv table for the manuscript table
reticulate::source_python('manuscript2/pfas_brics.py')
brics <-
  map2_dfr(frame$smiles, frame$ID, function(x,y)
    xx <- pfas_brics(x, '[$([S]);$([S]CCC(F)(F))]~[$([C]);$(CCC(F)(F))]') %>%
      matrix(ncol = 2, dimnames = list(c(),c( 'R1', 'R2'))) %>%
      as_tibble() %>%
      mutate(ID = y)
    )
frame <- left_join(frame, brics, by = "ID")
readr::write_excel_csv(select(frame, -R2), path = 'manuscript2/catawba_river_afff_rst.csv')

# get annotated spectra for plotting as examples
annotated_spectra <- map_dfr(frame$CD_ID, .f = function(id){
  ID <- dat %>% filter(CD_ID == !!id) %>% pull(ID)
  sirius_ann <- sirius_rst %>%
    filter(CD_ID == !!id) %>%
    select(tree) %>%
    unnest(tree) %>%
    rowwise() %>%
    mutate(
      mass = cdutils::get_mass(molecularFormula),
      mz = mass + adducts %>% filter(adduct == '[M-H]-') %>% pull(mass)
    ) %>%
    ungroup()

  cfm_ann <- cfmR::cfm_read_db(
    db_file = db_file,
    table_name = 'cfm_neg',
    ID = ID,
    return_annotation = T
  )[[1]] %>%
    unnest(data) %>%
    distinct(mz, intensity, .keep_all = T)

  ms <- spec[[paste(id)]]$ms2
  mz_range <- sapply(ms$mz, cdutils::get_massRange, 10) %>% t
  mz_range <- cbind(mz_range, ms$mz)

  cfm_spec <- apply(mz_range, 1, function(x){
    filter(cfm_ann, mz > x[2], mz < x[1]) %>%
      select(intensity, smiles) %>%
      mutate(mz = x[3], cfm_intensity = intensity)
  }) %>% bind_rows() %>%
    select(mz, smiles)

  sirius_spec <- apply(mz_range, 1, function(x){
    filter(sirius_ann, mz > x[2], mz < x[1]) %>%
      select(molecularFormula) %>%
      mutate(mz = x[3])
  }) %>%
    bind_rows() %>%
    transmute(mz, sirius_mf = molecularFormula)

  ms %>%
    left_join(y = sirius_spec, by = 'mz') %>%
    left_join(y = cfm_spec, by = 'mz') %>%
    as_tibble() %>%
    mutate_at(vars(smiles, sirius_mf), list(~as.numeric(!is.na(.)))) %>%
    unite(col = "Annotation", sirius_mf, smiles) %>%
    mutate(`Annotation` = factor(
      `Annotation`,
      levels = c('0_0', '1_0', '0_1', '1_1'),
      labels = c('none', 'SIRIUS', 'CFM', 'SIRIUS & CFM')
    )
    ) %>%
    mutate(
      CD_ID = !!id,
      ID = ID,
      intensity = intensity/max(intensity)
    )
})

d <- left_join(annotated_spectra, frame)
d <- unite(d, "label", CD_ID, R1, sep = ',', remove = F) %>%
  mutate(
    CD_ID = factor(CD_ID),
    label = factor(label)
    ) %>%
  mutate(
    label = fct_reorder(label, as.numeric(CD_ID))
  ) %>%
  group_by(CD_ID, label)

# plot the annotated spectra for figure 3
p <- ggplot(d, aes(mz, intensity, col = `Annotation`)) +
  geom_segment(aes(xend = mz, yend = 0), lineend = 'butt') +
  scale_color_brewer(palette = 'Set1') +
  xlab('m/z') +
  ylab('Normalized intensity (%)') +
  facet_wrap(~label, ncol = 3) +
  theme(
   axis.ticks.y = element_blank(),
   legend.position = 'right',
   legend.background = element_blank(),
   legend.key = element_blank()
   )

ggsave(grid.draw(shift_legend(p)),
       filename = 'manuscript2/figures/figure3_catawba_riv_unkns_msms.png',
       width = 7,
       height = 5,
       scale = 1.5
       )



#EXTRA CODE ####
# get areas by site
# samp_names <- cdutils::get_fileNames(CDresultFile = cd_rst_file) %>%
#   pull(Name) %>%
#   str_extract(pattern = '(?<=_)[:alnum:].+(?=.raw)')
#
# area <-
#   matrix(
#     sapply(dat$Area, cdutils::unpackAreas),
#     ncol = length(samp_names),
#     nrow = nrow(dat),
#     dimnames = list(dat$CD_ID, samp_names),
#     byrow = T
#   ) %>%
#   as_tibble(rownames = "CD_ID") %>%
#   filter_at(.vars = vars(-CD_ID), any_vars(. > 1e7)) %>%
#   group_by(CD_ID) %>%
#   pivot_longer(-CD_ID, values_to = "area", names_to = "sample") %>%
#   filter(!grepl('MBLK', sample)) %>%
#   ungroup() %>%
#   mutate(CD_ID = fct_reorder(CD_ID, area, .desc = F)) %>%
#   group_by(CD_ID, sample)
#
# ggplot(area, aes(area, CD_ID, label = CD_ID)) +
#   geom_point(aes(color = sample)) +
#   geom_line(aes(group = CD_ID)) +
#   scale_x_log10() +
#   geom_text_repel(
#     nudge_x = 0.15,
#     direction = 'y',
#     hjust = 0,
#     segment.size = 0.2
#   )
#
# library(ggrepel)
# p <- ggplot(dat, aes(y = CD_ID, x = 1, label = ID)) +
#   geom_point(color = "red") +
#   xlim(1,1.375) +
#   geom_text_repel(
#     nudge_x = 0.15,
#     direction = 'y',
#     hjust = 0,
#     segment.size = 0.2
#   )
#
#
# rxn_degree <- rxn_graph %>%
#   activate("nodes") %>%
#   transmute(InChIKey = name, rxn_degree = centrality_degree(loops = F)) %>%
#   as_tibble()
#
# i <- mf_graphs %>%
#   activate("nodes") %>%
#   mutate(id = row_number()) %>%
#   filter(`Ion type` == 'precursor') %>%
#   pull(id)
#
# msms_degree <- mf_graphs %>%
#   activate("nodes") %>%
#   mutate(id = row_number()) %>%
#   activate("edges") %>%
#   filter(!edge_is_loop()) %>%
#   activate("nodes") %>%
#   mutate(
#     iso = node_is_isolated(),
#     degree_out = centrality_degree(loop = F, mode = "out"),
#     degree_in = centrality_degree(mode = "in", loop = F),
#     adj = node_is_adjacent(i)
#   ) %>%
#   filter((`Ion type` == 'precursor' & degree_out >= 1 & !iso) | (adj & degree_in >= 2)) %>%
#   mutate(
#     iso = node_is_isolated()
#   ) %>%
#   filter(!iso) %>%
#   activate("nodes") %>%
#   #filter(`Ion type` == 'precursor') %>%
#   mutate(CD_ID = as.integer(label))
#
# ggraph(msms_degree, layout = 'nicely') +
#   geom_edge_link0(color = 'lightgray',  alpha = 0.5) +
#   geom_node_point(aes(color = `Ion type`)) +
#   theme_graph()
#
#
#
#
# # make a heatmap of structure similarity
# fp_sim <- matrix(
#   unlist(FpsSimMat(dat$InChI)),
#   ncol = nrow(dat),
#   nrow = nrow(dat),
#   dimnames = list(dat$CD_ID, dat$CD_ID)
# )
#
# m <- 1 - fp_sim
# dist <- as.dist(m)
# hc <- hclust(dist)
# fp_dend <- as.dendrogram(hc)
# fp_m <- fp_sim[hc$order, hc$order]
# fp_m[lower.tri(fp_m)] <- NA
#
# fp_df <- as.data.frame(fp_m)
# fp_df$var1 <- rownames(fp_df)
# fp_df <- na.omit(melt(fp_df, 'var1', variable.name='var2', value.name = 'Tanimoto similarity'))
# fp_df$var1 <- factor(fp_df$var1, levels = hc$labels[hc$order], ordered = T)
# fp_df$var2 <- factor(fp_df$var2, levels = hc$labels[hc$order], ordered = T)
#
# fp_heatmap <- ggplot(fp_df, aes(var1, var2)) +
#   theme_minimal() +
#   geom_tile(aes(fill = `Tanimoto similarity`)) +
#   scale_fill_gradient() +
#   theme(
#     axis.text.y = element_text(size = 6),
#     axis.text.x = element_blank(),
#     axis.title = element_blank(),
#     legend.position = "bottom",
#     panel.grid = element_blank()
#   )
#
#
# # make a heamtap of ms/ms similarity
# spec <- cdutils::get_BestHit(
#   CDresultFile = cd_rst_file,
#   Checked = T,
#   BackgroundStatus = F,
#   MSnStatus = T
# )
#
# mgf_file <- "manuscript2/catawba_river_afff.mgf"
# map(spec, function(x) {
#   sp2 <- new(
#     'Spectrum2',
#     mz = x$ms2$mz,
#     intensity = x$ms2$intensity,
#     centroided = T,
#     precursorMz = x$mz,
#     scanIndex = x$ID
#   )
#   meta_dat <- mutate(x$featureData, CD_ID = x$ID) %>% DataFrame()
#   spectra <- Spectra(sp2, elementMetadata = meta_dat)
#   tmpf <- tempfile()
#   writeMgfData(object = spectra, tmpf)
#   mgf <- readLines(tmpf)
#   write(mgf, file = mgf_file, append = T)
#   unlink(tmpf)
# })
#
# exp <- readMgfData(mgf_file)
# spec_sim <- compareSpectra(exp, fun = 'dotproduct', bin = 0.01)
# dimnames(spec_sim) <- list(as.character(fData(exp)$CD_ID), as.character(fData(exp)$CD_ID))
#
# m <- 1 - spec_sim
# dist <- as.dist(m)
# hc <- hclust(dist)
# spec_dend <- as.dendrogram(hc)
# spec_m <- spec_sim[hc$order, hc$order]
# spec_m[lower.tri(spec_m)] <- NA
#
# spec_df <- as.data.frame(spec_m)
# spec_df$var1 <- rownames(spec_df)
# spec_df <- na.omit(melt(spec_df, 'var1', variable.name='var2', value.name = 'Dotproduct similarity'))
# spec_df$var1 <- factor(spec_df$var1, levels = hc$labels[hc$order], ordered = T)
# spec_df$var2 <- factor(spec_df$var2, levels = hc$labels[hc$order], ordered = T)
#
# spec_heatmap <- ggplot(spec_df, aes(var1, var2)) +
#   theme_minimal() +
#   geom_tile(aes(fill = `Dotproduct similarity`)) +
#   scale_fill_gradient(high = "yellowgreen", low = "forestgreen") +
#   theme(
#     axis.text.y = element_text(size = 6),
#     axis.text.x = element_blank(),
#     axis.title = element_blank(),
#     legend.position = "bottom",
#     panel.grid = element_blank()
#   )

# library(dendextend)
# dlist <- dendlist("Tanimoto" = fp_dend, "Dotproduct" = spec_dend) %>% untangle(method = "step1side")
# tan_dot_tangle <- tanglegram(dlist,
#     highlight_distinct_edges = F,
#     common_subtrees_color_branches = T,
#     highlight_branches_lwd = F,
#     lwd = 1
#     )



# mf <- sapply(as_tibble(mf_graphs)$name, mf_expr, USE.NAMES = F)
# #mf_graphs <- set.vertex.attribute(graph = mf_graphs, name = 'mf', value =  mf)
# form_graph <- mf_graphs %>%
#   activate("nodes") %>%
#   mutate(mf = ifelse(label == "", "", name)) %>%
#    ggraph(layout = ly) +
#   geom_node_text(aes(label = mf), check_overlap = F, repel = F, color = 'black', size = 3, parse = T) +
#   theme_graph(background = "white") +
#   theme(legend.position = 'bottom') +
#   labs(subtitle = 'Molecular formulas')
# rxn_graph2 <- rxn_graph %>%
#   activate("edges") %>%
#   filter(from != to, !is.na(rxn_type)) %>%
#   activate("nodes") %>%
#   mutate(iso = node_is_isolated()) %>%
#   filter(!iso)
#
# trans_graph <- ggraph(rxn_graph2, layout = 'nicely') +
#   geom_node_point(aes(size = `Max area`, color = factor(`Fluorine count`)), alpha = 0.75) +
#   F_colScale +
#   scale_size(trans = "log10") +
#   geom_edge_link(aes(col = factor(rxn_type, levels = c('bio', 'bio|hyd', 'MS/MS'))), arrow = arrow(type = "closed", length = unit(0.05, "inches")), n = 4) +
#   scale_edge_color_manual(values = c(brewer.pal(n = 3, 'Set1'), 'black'), guide_edge_colorbar(title = 'Trans. type')) +
#   geom_node_text(aes(label = label), check_overlap = T, repel = T, color = 'gray66', size = 2) +
#   theme_graph()
