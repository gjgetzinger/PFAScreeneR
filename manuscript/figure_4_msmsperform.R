# plot the performance of spectral library matching by various techniques
library(MSnbase)
library(purrr)
library(tibble)
library(pbapply)

ms <-
  readMgfData(list.files(
    system.file(package = 'PFAScreeneR'),
    pattern = 'mgf',
    full.names = T
  ))

key <- fData(object = ms) %>% as_tibble()
ids <- pull(key, ID) %>% as.character()

search_rst <- list(
  #`Identifier` = sim_id_search,
  `Mol. form.` = sim_molform_search,
  `Accurate neutral mass` = sim_neutralmass_search,
  `Accurate precur. m/z` = sim_precurmz_search,
  `Nominal neutral mass` = sim_nominalneutralmass_search,
  `Nominal precur. m/z` = sim_nominalprecurmz_search
) %>%
  lapply(setNames, nm = ids)

search_summary <- pblapply(names(search_rst), function(s){
  map_dfr(names(search_rst[[s]]), function(id) {
    x <- search_rst[[s]][[id]]
    if (length(x) == 0) {
      enframe(id, name = NULL, value = "ID")
    } else {
      sim <- replace_na(x$dotproduct_similarity, 0)
      if(length(unique(sim))>1) {
        sim_scale <- as.vector(scale(sim, min(sim), diff(range(sim))))
      } else {
        sim_scale <- sim
      }
      x$sim_scale <- sim_scale
      arrange(x, desc(dotproduct_similarity)) %>%
        mutate(
          rank = dense_rank(desc(dotproduct_similarity)),
          correct = ID %in% !!id
          ) %>%
        transmute(
          ID, correct, dotproduct_similarity, sim_scale, rank, query_id = id
        )
    }
  })
}) %>%
  setNames(nm = names(search_rst)) %>%
  bind_rows(.id = 'Search') %>%
  mutate(Search = factor(
    Search,
    levels = c(
      'Mol. form.',
      'Accurate neutral mass',
      'Nominal neutral mass',
      'Accurate precur. m/z',
      'Nominal precur. m/z',
      'Identifier'
    )
  ))

# number of structures by different query
mols_per_query <- search_summary %>%
  filter(Search != 'Identifier') %>%
  group_by(Search, query_id) %>%
  tally(name = 'N mols. per query', sort = T)

# % of test set vs rank of true structure
percent_by_rank <- search_summary %>%
  filter(correct, Search != 'Identifier') %>%
  mutate(rank = factor(
    ifelse(rank >= 20 , 21, rank),
    levels = c(1:21),ordered = T,
    labels = c(as.character(1:20), '>20')
  )) %>%
  group_by(Search, rank) %>%
  tally(name = 'ct') %>%
  mutate(pct = 100 * (cumsum(ct) / length(ids)))


a <- ggplot(group_by(search_summary, Search, query_id) %>% filter(!is.na(query_id)),
            aes(x = sim_scale, y = correct, color = Search)) +
  scale_color_brewer(palette = "Set1", guide = F) +
  geom_boxplot(outlier.shape = 16, outlier.alpha = 0.25, outlier.size = 0.5) +
  theme_minimal() +
  labs(x = 'Scaled dot product similarity', y = "Correct annotation")

b <- ggplot(mols_per_query,
            aes(x = `N mols. per query`, y = Search, color = Search)
            ) +
  scale_color_brewer(palette = "Set1") +
  scale_x_log10() +
  geom_boxplot(show.legend = F) +
  theme_minimal() +
  theme(axis.text.y = element_blank())

c <- ggplot(percent_by_rank, aes(rank, pct, color = Search, group = Search)) +
  scale_color_brewer(palette = "Set1") +
  geom_point() +
  geom_line(show.legend = F) +
  labs(x = 'Rank of true structure', y = 'Fraction of test set (%)') +
  theme_minimal() +
  theme(legend.position = 'none')

figure4 <- ggarrange(a,b,c,
                     nrow = 1, ncol = 3,
                     common.legend = T,
                     legend = "top", labels = "auto")

ggsave(
  'manuscript/figures/figure4_msms_sim_dist.png',
  figure4,
  device = 'png',
  dpi = 300,
  width = 7,
  height = 2, scale = 1.5,
  units = 'in'
)

# eg_spec <- bind_cols(
#   fData(pfas_exp_msms) %>% as_tibble(),
#   enframe(peaksCount(pfas_exp_msms), value = 'peak_count', name = NULL)
# ) %>%
#   left_join(
#     y = filter(search_summary,correct, Search == 'Identifier', rank == 1),
#     by = "ID"
#   ) %>%
#   filter(peak_count > 5, !grepl('in-source', NAME)) %>%
#   group_by(PRECURSOR_TYPE) %>%
#   arrange(desc(dotproduct_similarity), desc(peak_count)) %>%
#   distinct(PRECURSOR_TYPE, .keep_all = T) %>%
#   select(NAME, SMILES, ID, PRECURSOR_TYPE, peak_count, dotproduct_similarity, Search)
#
# lapply(1:nrow(eg_spec), function(x){
#   exp_spec <- pfas_exp_msms[[which(fData(pfas_exp_msms)$ID == pull(eg_spec[x,'ID']))]]
#   pred_spec <- search_rst[[pull(eg_spec[x,'Search'])]][[pull(eg_spec[x,'ID'])]] %>% filter(ID == pull(eg_spec[x,'ID'])) %>% pull(predicted_spectrum)
#   pred_spec <- pred_spec[[1]]
#
#
#   s <- bind_rows(
#     tibble(mz = mz(exp_spec), intensity = intensity(exp_spec)) %>%
#       mutate(intensity = 100*(intensity/max(intensity))),
#     tibble(mz = mz(pred_spec), intensity = intensity(pred_spec)) %>%
#       mutate(intensity = -100*(intensity/max(intensity)))
#   ) %>%
#     mutate(color = c("predicted", "experimental")[as.numeric(factor(sign(intensity)))])
#
#   ggplot(s, aes(x = mz, y = intensity, color = color)) +
#     scale_color_brewer(palette = "Set1") +
#     geom_segment(aes(xend = mz, yend = 0)) +
#     theme_minimal() +
#     theme(legend.position = c(0.1,0.85), legend.title = element_blank())
#
#
# })









