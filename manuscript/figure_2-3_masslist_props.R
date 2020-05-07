# plot list membership and molecular properties of mass list entries
library(DBI)
library(dplyr)
library(eulerr)
library(stringr)
library(ggplot2)
library(tidyr)

dat <- dbConnect(RSQLite::SQLite(), 'inst/PFAScreeneR.db') # change to system.file

## make a venn diagram of list overlap
list_membership <- tbl(dat, 'list_membership') %>%
  as_tibble() %>%
  mutate(
    `List type` = factor(
      List,
      levels = c('HIGGINS', 'PFASTRIER', 'PFASKEMI',
                 'PFASEPA', 'HYD', 'PFASNORMAN',
                 'BIOTRANS', 'PFASLIU', 'OECDPFAS'),
      labels = c('In-house AFFF', 'EPA',
                 'EPA', 'EPA', 'Hyd. products', 'EPA',
                 'Biotrans. products', 'EPA', 'EPA')
      ))

unique_inchikey <- distinct(list_membership, InChIKey) %>% pull()
lists <- split(list_membership$InChIKey, list_membership$`List type`)
m <- sapply(lists, function(x){unique_inchikey %in% x})
m <- m[,c('EPA', 'In-house AFFF', 'Hyd. products', 'Biotrans. products')]
v <- venn(m)

figure2 <- plot(
    v,
    adjust_labels = T,
    legend = list(side = 'bottom', ncol = 2, nrow = 2),
    cex = 0.1
  )

ggsave(
  plot = figure2,
  filename = 'figures/figure2_listvenn.png',
  width = 2,
  height = 2,
  scale = 2,
  units = 'in',
  device = 'png',
  dpi = 300
)


## distributions of some molecular properties (by list)
props <- left_join(list_membership,
          tbl(dat, 'molecules') %>% as_tibble(),
          by = 'InChIKey') %>%
  left_join(
    y = tbl(dat, 'mol_props') %>% as_tibble(),
    by = 'ID'
    ) %>%
  rowwise %>%
  mutate(
    `F count` = as.integer(str_extract(MolForm, '(?<=F)[0-9]{1,100}')), # ext from string gives NA when no alapha after letter
    `C count` = as.integer(str_extract(MolForm, '(?<=C)[0-9]{1,100}'))
    ) %>%
  ungroup() %>%
  mutate(
    `F count` = replace_na(`F count`, replace = 1),
    `C count` = replace_na(`C count`, replace = 1)
    #,`Kendrick mass (CF2)` = ExactMass * (50/49.99681)
  ) %>%
  mutate(
    `Nominal mass` = round(ExactMass)
    #, # trunc? floor? `Nominal Kendrick mass (CF2)` = round(`Kendrick mass (CF2)`)
  ) %>%
  mutate(
    `Mass defect` = ExactMass - `Nominal mass`
    #,`Kendrick mass defect (CF2)` = `Kendrick mass (CF2)` - `Nominal Kendrick mass (CF2)`
  ) %>%
  transmute(
    `List` = `List type`,
    `Exact mass` = ExactMass,
    `C count`,
    `F count`,
    `Mass defect`
    ) %>%
  gather(key = "measure", value = "value",
                     -`List`) %>%
  mutate(measure = factor(
    measure,
    levels = c('Exact mass', 'C count', 'F count', 'Mass defect')
  ))

figure3 <- ggplot(props, aes(value)) +
  facet_grid(`List`~measure, scales = c('free')) +
  geom_histogram(bins = 100) +
  theme_minimal() +
  theme(
    #panel.grid.major.x = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_rect(
      fill = 'lightgray', color = 'lightgray'
      )
  )

ggsave(
  figure3,
  filename = 'manuscript/figures/figure3_prophist.png',
  device = 'png',
  width = 7,
  height = 3.5,
  units = 'in',
  dpi = 300, scale = 1.5
)

