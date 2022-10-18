# need to also install the following packages
if(!require("xfun")) install.packages("xfun")
options(repos = c(
  rtrees = 'https://daijiang.r-universe.dev',
  CRAN = 'https://cloud.r-project.org')
)
xfun::pkg_attach2(c("rtrees", "phylocomr", "microbenchmark", "tidyverse"))
# I cannot use V.PhyloMaker::phylo.maker() and was forced to library it, 
#     which is a little bit annoying.
library(V.PhyloMaker) 

d = read_delim("Data/plants_phylomatic_names.txt", "/", col_names = F) %>%
  set_names(c("family", "genus", "species")) %>%
  mutate(genus = str_extract(genus, "^[^_]+"),
         species = rtrees:::cap_first_letter(species),
         genus = rtrees:::cap_first_letter(genus),
         family = rtrees:::cap_first_letter(family)) %>%
  select(species) %>% distinct()
n_distinct(d$species)

#' Function to test speed
#' 
#' @param n_sp_missing The number of missing species to be binded
#' @param n_sp_in The number of species that were already in the megatree
#' @param n_times How many times to compare? Default of 5 should be enough.
#' @return A data frame with the test methods, and their time used.
#' 
speed_test = function(n_sp_missing, n_sp_in = 500, n_times = 5){
  test_sp = c(sample(sp_in, n_sp_in), sample(sp_missing, n_sp_missing))
  test_sp_df = rtrees::sp_list_df(sp_list = test_sp, taxon = "plant")
  test_sp_df = filter(test_sp_df, !is.na(family))
  test_sp_df2 = mutate(test_sp_df, phylom_sp = paste(family, genus, species, sep = "/"),
                       phylom_sp = casefold(phylom_sp))
  # compare speed
  xx2 = microbenchmark::microbenchmark(
    Phylomatic = phylocomr::ph_phylomatic(taxa = test_sp_df2$phylom_sp, phylo = megatrees::tree_plant_otl),
    rtrees = get_tree(sp_list = test_sp_df, taxon = "plant"),
    V.PhyloMaker = phylo.maker(sp.list = test_sp_df, scenarios = "S1"),
    times = n_times
  )
  
  xx2 = as.data.frame(xx2) %>% 
    mutate(n_sp_in = n_sp_in, n_sp_missing = n_sp_missing)
  
  xx2
}

# all tests have 500 species already in the megatree, then
# start with 50 missing species, and 500, 1000, ..., 5000
speed_out = purrr::map_dfr(c(50, seq(550, 5500, by = 500)-50), speed_test)
# for some reason, phylocomr::ph_phylomatic does not work (no mssing species was inserted)
# I thus removed it from the results
saveRDS(filter(speed_out, expr != "Phylomatic"), "Data/rtrees_speed_out.rds")


## fish ----

taxa_fish = read_csv("Data/taxa_fish.csv") %>% 
  mutate(species = rtrees:::cap_first_letter(taxa_fish),
         species = stringr::str_trim(species)) %>% 
  select(-taxa_fish) %>% 
  distinct()
taxa_fish = rtrees::sp_list_df(taxa_fish, "fish")
taxa_fish = drop_na(taxa_fish, family)
sp_missing = setdiff(taxa_fish$species, megatrees::tree_fish_12k$tip.label)
sp_in = intersect(taxa_fish$species, megatrees::tree_fish_12k$tip.label)
test_sp = unique(c(sample(sp_in, 50), sample(sp_missing, 1500)))
test_sp_df = rtrees::sp_list_df(sp_list = test_sp, taxon = "fish")
test_sp_df = filter(test_sp_df, !is.na(family))


library(FishPhyloMaker)

fpmaker = FishTaxaMaker(test_sp_df$species, allow.manual.insert = FALSE)
str(fpmaker)

system.time(suppressMessages(res_phylo <- FishPhyloMaker(data = fpmaker$Taxon_data_FishPhyloMaker,
                                                         insert.base.node = TRUE, 
                                                         return.insertions = TRUE, 
                                                         progress.bar = TRUE)))
res_phylo$Phylogeny
# user  system elapsed 
# 464.996  13.495 481.044 

system.time(res_rtree <- get_tree(sp_list = test_sp, taxon = "fish"))
# user  system elapsed 
# 2.568   1.685   4.048 

