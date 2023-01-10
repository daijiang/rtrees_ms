# required packages
if(!require(rtrees))
  install.packages('rtrees', repos = c(
    rtrees = 'https://daijiang.r-universe.dev', 
    CRAN = 'https://cloud.r-project.org'
  ))

if(!require(xfun)) install.packages("xfun")

xfun::pkg_attach2("tidyverse", "ape", "rtrees", "picante", "parallel", "furrr", "patchwork")

########### Notes ##############
# The functions below have three scenarios, which was the case before
# v1.0.1 of `rtrees`. Based on the simulations results, there is no 
# need of one of the grafting scenario ("at_or_above_basal") and it was
# removed from `rtrees`.
################################


#' Simulate and estimate the effect of different grafting scenarios on phylogenetic diversity analyses.
#' 
#' This function will simulate 
#' 
#' @param seed The seed for simulations.
#' @param n_sp The number of species in the community.
#' @param n_bind The number of species to drop and then graft back.
#' @param tr The megatree.
#' @param taxon The taxonomic group of the species list.
#' @param n_random_tree The number of trees to generated when the scenario is "random_below_basal".
#' @return A data frame of phylogenetic diversity based on different grafting scenarios.
#' 
eval_bind_loc = function(seed = 1, n_sp = 150, n_bind = 50, 
                         tr = megatrees::tree_fish_12k, 
                         taxon = "fish", n_random_tree = 200){
  set.seed(seed)
  # to work with a smaller tree to save time
  tr_sub = ape::keep.tip(tr, sample(tr$tip.label, 3000))
  # tr_sub = rtrees::add_root_info(tr_sub,
  #                                rtrees::classifications[rtrees::classifications$taxon == taxon,],
  #                                show_warning = FALSE)
  sp_list = sample(tr_sub$tip.label, size = n_sp)
  sp_to_bind  = sample(sp_list, size = n_bind)
  # remove these species from the megatree
  tr1 = ape::drop.tip(tr_sub, tip = sp_to_bind)
  tr1$genus_family_root = NULL
  tr1 = ape::makeNodeLabel(tr1)
  
  # tr1$node.label
  tr2 = rtrees::add_root_info(tr1,
                              rtrees::classifications[rtrees::classifications$taxon == taxon,],
                              show_warning = FALSE)
  # tr_at_basal = rtrees::get_tree(sp_list, tree = tr2, taxon = taxon, scenario = "at_basal_node")
  
  # generate phylogenies
  tr_at_basal = rtrees::get_tree(sp_list, tree = tr2, taxon = taxon, 
                                 scenario = "at_basal_node", tree_by_user = FALSE)
  
  tr_at_or_above_basal = suppressMessages(rtrees::get_tree(sp_list, tree = tr2, taxon = taxon, 
                                                           scenario = "at_or_above_basal", tree_by_user = FALSE))
  # random locations
  tr_random_below_basal = suppressMessages(purrr::map(1:n_random_tree, rtrees::get_tree, 
                                                      sp_list = sp_list, tree = tr2, taxon = taxon, 
                                                      scenario = "random_below_basal", 
                                                      tree_by_user = FALSE, show_grafted = FALSE))
  tr_random_below_basal = tr_random_below_basal[purrr::map_chr(tr_random_below_basal, class) == "phylo"]
  
  if(any(a_na <- map_lgl(tr_random_below_basal, function(x) any(is.na(x$edge.length))))){
    tr_random_below_basal = tr_random_below_basal[-which(a_na)]
  }
  
  tr_sub2 = ape::keep.tip(tr_sub, tr_at_basal$tip.label)
  
  if(!dir.exists("sim_output")) dir.create("sim_output")
  save(tr_sub2, tr_at_basal, tr_at_or_above_basal, tr_random_below_basal, 
       file = paste0("sim_output/trs_", seed, ".RData"))
  
  # purrr::map_chr(tr_random_below_basal, class)
  
  # message("done with rtrees")
  
  # communitdata
  comm = as.data.frame(matrix(1, nrow = 1, ncol = n_sp))
  rownames(comm) = "site1"
  colnames(comm) = sp_list
  comm = comm[, tr_at_basal$tip.label]
  
  # PD
  x_r = lapply(tr_random_below_basal, function(i){
    picante::pd(comm, tree = i)
  }) |> 
    bind_rows() |> 
    colMeans(na.rm = TRUE)
  x_r = dplyr::mutate(as.data.frame(t(x_r)), type = "random_below_basal_ave")
  out_pd = bind_rows(# "true" PD
    dplyr::mutate(picante::pd(comm, tr_sub2), type = "true"),
    dplyr::mutate(picante::pd(comm, tr_at_basal), type = "at_basal_node"),
    dplyr::mutate(picante::pd(comm, tr_at_or_above_basal), type = "at_or_above_basal"),
    x_r
  )
  
  # MPD
  x_r = lapply(tr_random_below_basal, function(i){
    picante::mpd(comm, cophenetic(i))
  }) |> 
    unlist() |> 
    mean(na.rm = T)
  out_mpd = tibble(# "true" MPD
    mpd = c(
      picante::mpd(comm, cophenetic(tr_sub2)),
      picante::mpd(comm, cophenetic(tr_at_basal)),
      picante::mpd(comm, cophenetic(tr_at_or_above_basal)),
      x_r
    ),
    type = c("true", "at_basal_node", "at_or_above_basal", "random_below_basal_ave")
  )
  
  # MNTD
  x_r = lapply(tr_random_below_basal, function(i){
    picante::mntd(comm, cophenetic(i))
  }) |> 
    unlist()
  x_r = mean(x_r[is.finite(x_r)], na.rm = T)
  out_mntd = tibble(# "true" MNTD
    mntd = c(
      picante::mntd(comm, cophenetic(tr_sub2)),
      picante::mntd(comm, cophenetic(tr_at_basal)),
      picante::mntd(comm, cophenetic(tr_at_or_above_basal)),
      x_r
    ),
    type = c("true", "at_basal_node", "at_or_above_basal", "random_below_basal_ave")
  )
  
  out = left_join(out_pd, out_mpd, by = "type") |> 
    left_join(out_mntd, by = "type")
  
  dplyr::select(out, type, sr = SR, pd = PD, mpd, mntd)
}

if(!file.exists("Data/x_pd.rds")){
  plan(multisession, workers = 50)
  
  x_pd = future_map(1:1000, .f = function(i){
    try(eval_bind_loc(seed = i))
  }, .progress = TRUE)
  
  saveRDS(x_pd, "Data/x_pd.rds")
}


# x_pd = readRDS("Data/x_pd.rds")


#' Function to estimate the effects of different grafting scenarios on phylogenetic signal of traits
#' 
#' The phylogenies were already generated and stored above. 
#' The phylogenetic signal calculations were not included in `eval_bind_loc()`
#' because `phylotools::phylosig()` used multiple cores inside. 
#' If I included this part there, nothing will be finished with 
#' limited computer resources ...
#' 
get_phy_sig_sim = function(i){
  f = paste0("sim_output/trs_", i, ".RData")
  load(f)
  
  cat("n_rand =", length(tr_random_below_basal), "\n")
  # phy signal
  trait_sim_BM = ape::rTraitCont(tr_sub2, model = "BM", sigma = 1)
  
  get_phy_signal = function(v, phy, v.name){
    if(is.null(names(v))) names(v) = v.name
    xy_k = phytools::phylosig(tree = phy, x = v, method = "K", test = FALSE)
    xy_l = phytools::phylosig(tree = phy, x = v, method = "lambda", test = FALSE)
    tibble::tibble(method = c("K", "lambda"), 
                   statistics = c(xy_k, xy_l$lambda))
  }
  
  ps_true = get_phy_signal(v = trait_sim_BM, phy = tr_sub2)
  ps_basal = get_phy_signal(v = trait_sim_BM, phy = tr_at_basal)
  ps_above_basal = get_phy_signal(v = trait_sim_BM, phy = tr_at_or_above_basal)
  
  ps_random = map_dfr(tr_random_below_basal, function(i){
    try(get_phy_signal(v = trait_sim_BM, phy = i))
  })
  # sum(purrr::map(ps_random, class) == "try-error")
  # which(purrr::map(ps_random, class) == "try-error")
  ps_random_ave = group_by(ps_random, method) |> 
    summarise(statistics = mean(statistics, na.rm = T))
  
  out_ps1 = bind_rows(# "true" PS
    dplyr::mutate(ps_true, type = "true"),
    dplyr::mutate(ps_basal, type = "at_basal_node"),
    dplyr::mutate(ps_above_basal, type = "at_or_above_basal"),
    dplyr::mutate(ps_random_ave, type = "random_below_basal_ave")
  )
  
  out_ps2 = pivot_wider(out_ps1, names_from = "type", values_from = "statistics") |> 
    mutate(diff_at_basal_node_true = at_basal_node - true,
           diff_at_or_above_basal_true = at_or_above_basal - true,
           diff_random_below_basal_ave_true = random_below_basal_ave - true
    ) |> 
    dplyr::select(method, true, starts_with("diff"))
  
  out_ps2
}

# get_phy_sig_sim(1)

# run it
if(!file.exists("Data/x_ps_1_1000.rds")){
  x_ps_1_1000 = vector("list", 1000)
  for(i in 1:1000){
    cat("i =", i, "\t")
    x_ps_1_1000[[i]] = try(get_phy_sig_sim(i))
    if(i %% 10 == 0) 
      saveRDS(x_ps_1_1000, file = "Data/x_ps_1_1000.rds")
  }
}

## data prep -------------------------

### Phylo diversity
x = readRDS("Data/x_pd.rds")
x = x[-which(map_chr(x, class) == "try-error")]

x_pd = map_dfr(x, function(i) {
  dplyr::select(i, type, pd) |> 
    pivot_wider(names_from = "type", values_from = "pd") |> 
    mutate(diff_at_basal_node_true = at_basal_node - true,
           diff_at_or_above_basal_true = at_or_above_basal - true,
           diff_random_below_basal_ave_true = random_below_basal_ave - true
    ) |> 
    dplyr::select(true, starts_with("diff"))
}, .id = "sim")

x_mpd = map_dfr(x, function(i) {
  dplyr::select(i, type, mpd) |> 
    pivot_wider(names_from = "type", values_from = "mpd") |> 
    mutate(diff_at_basal_node_true = at_basal_node - true,
           diff_at_or_above_basal_true = at_or_above_basal - true,
           diff_random_below_basal_ave_true = random_below_basal_ave - true
    ) |> 
    dplyr::select(true, starts_with("diff"))
}, .id = "sim")

x_mntd = map_dfr(x, function(i) {
  dplyr::select(i, type, mntd) |> 
    pivot_wider(names_from = "type", values_from = "mntd") |> 
    mutate(diff_at_basal_node_true = at_basal_node - true,
           diff_at_or_above_basal_true = at_or_above_basal - true,
           diff_random_below_basal_ave_true = random_below_basal_ave - true
    ) |> 
    dplyr::select(true, starts_with("diff"))
}, .id = "sim")

all_diff = bind_rows(
  x_pd |> 
    set_names(c("sim", "true", "at_basal_node", "at_or_above_basal", "random_below_basal")) |> 
    dplyr::select(-sim, -true) |> 
    pivot_longer(cols = 1:3, names_to = "Scenarios", values_to = "Differences") |> 
    mutate(diversity = "Faith's PD"),
  x_mpd |> 
    set_names(c("sim", "true", "at_basal_node", "at_or_above_basal", "random_below_basal")) |> 
    dplyr::select(-sim, -true) |> 
    pivot_longer(cols = 1:3, names_to = "Scenarios", values_to = "Differences") |> 
    mutate(diversity = "MPD"),
  x_mntd |> 
    set_names(c("sim", "true", "at_basal_node", "at_or_above_basal", "random_below_basal")) |> 
    dplyr::select(-sim, -true) |> 
    pivot_longer(cols = 1:3, names_to = "Scenarios", values_to = "Differences") |> 
    mutate(diversity = "MNTD") 
)

all_raw = bind_rows(
  x_pd |> 
    set_names(c("sim", "true", "at_basal_node", "at_or_above_basal", "random_below_basal")) |> 
    mutate(at_basal_node = true + at_basal_node, 
           at_or_above_basal = true + at_or_above_basal,
           random_below_basal = true + random_below_basal) |> 
    dplyr::select(-sim) |> 
    pivot_longer(cols = 2:4, names_to = "Scenarios", values_to = "Diversity") |> 
    mutate(diversity = "Faith's PD"),
  x_mpd |> 
    set_names(c("sim", "true", "at_basal_node", "at_or_above_basal", "random_below_basal")) |> 
    mutate(at_basal_node = true + at_basal_node, 
           at_or_above_basal = true + at_or_above_basal,
           random_below_basal = true + random_below_basal) |> 
    dplyr::select(-sim) |> 
    pivot_longer(cols = 2:4, names_to = "Scenarios", values_to = "Diversity") |> 
    mutate(diversity = "MPD"),
  x_mntd |> 
    set_names(c("sim", "true", "at_basal_node", "at_or_above_basal", "random_below_basal")) |> 
    mutate(at_basal_node = true + at_basal_node, 
           at_or_above_basal = true + at_or_above_basal,
           random_below_basal = true + random_below_basal) |> 
    dplyr::select(-sim) |> 
    pivot_longer(cols = 2:4, names_to = "Scenarios", values_to = "Diversity") |> 
    mutate(diversity = "MNTD") 
)

### phylo signal
x_ps_1_1000 = readRDS("Data/x_ps_1_1000.rds")

x_ps_1_1000_df = bind_rows(x_ps_1_1000[map(x_ps_1_1000, class) != "try-error"])

all_raw_ps = 
  x_ps_1_1000_df |> 
  set_names(c("method", "true", "at_basal_node", "at_or_above_basal", "random_below_basal")) |> 
  mutate(at_basal_node = true + at_basal_node, 
         at_or_above_basal = true + at_or_above_basal,
         random_below_basal = true + random_below_basal) |> 
  pivot_longer(cols = 3:5, names_to = "Scenarios", values_to = "Phylo_signal")

# Plot for reviewers (3 scenarios) ----
pd_faith_lm = ggplot(filter(all_raw, diversity == "Faith's PD"), 
                     aes(x = true, y = Diversity, color = Scenarios, group = Scenarios)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_smooth(method = "lm") +
  labs(x = "'True' diversity value",
       y = "Estimated diversity value",
       title = "(A) Faith's PD") +
  cowplot::theme_cowplot() +
  theme(legend.position = "none") 

pd_faith_hist = ggplot(filter(all_diff, diversity == "Faith's PD"), 
                       aes(x = Differences, fill = Scenarios)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(alpha = 0.5, position = "identity") +
  labs(y = "Count",
       x = "Estimated - 'true'") +
  theme(legend.position = "none")

pd_faith = pd_faith_lm + inset_element(pd_faith_hist, left = 0.51, bottom = 0.01, right = 0.9, top =  0.38)


pd_mpd_lm = ggplot(filter(all_raw, diversity == "MPD"), 
                   aes(x = true, y = Diversity, color = Scenarios, group = Scenarios)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_smooth(method = "lm") +
  labs(x = "'True' diversity value",
       y = "Estimated diversity value",
       title = "(B) MPD") +
  cowplot::theme_cowplot() +
  theme(legend.position = c(0.1, 0.8)) 

pd_mpd_hist = ggplot(filter(all_diff, diversity == "MPD"), 
                     aes(x = Differences, fill = Scenarios)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(alpha = 0.5, position = "identity") +
  labs(y = "Count",
       x = "Estimated - 'true'") +
  theme(legend.position = "none")

pd_mpd = pd_mpd_lm + inset_element(pd_mpd_hist, left = 0.51, bottom = 0.01, right = 0.9, top =  0.38)

pd_mntd_lm = ggplot(filter(all_raw, diversity == "MNTD"), 
                    aes(x = true, y = Diversity, color = Scenarios, group = Scenarios)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_smooth(method = "lm") +
  labs(x = "'True' diversity value",
       y = "Estimated diversity value",
       title = "(C) MNTD") +
  cowplot::theme_cowplot() +
  theme(legend.position = "none") 

pd_mntd_hist = ggplot(filter(all_diff, diversity == "MNTD"), 
                      aes(x = Differences, fill = Scenarios)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(alpha = 0.5, position = "identity") +
  labs(y = "Count",
       x = "Estimated - 'true'") +
  theme(legend.position = "none")

pd_mntd = pd_mntd_lm + inset_element(pd_mntd_hist, left = 0.51, bottom = 0.01, right = 0.9, top =  0.38)

pd_all = pd_faith + pd_mpd + pd_mntd

# phylogenetic signals
ps_k_lm = ggplot(filter(all_raw_ps, method == "K"), 
                 aes(x = true, y = Phylo_signal, color = Scenarios, group = Scenarios)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_smooth(method = "lm") +
  labs(x = "'True' phylogenetic signal value",
       y = "Phylogenetic signal value (K)",
       title = "(D) K") +
  cowplot::theme_cowplot() +
  theme(legend.position = c(0.6, 0.2)) 

ps_k_hist = dplyr::select(x_ps_1_1000_df, -true) |> 
  pivot_longer(cols = 2:4, names_to = "Scenarios", values_to = "Differences") |> 
  filter(method == "K") |> 
  ggplot(aes(x = Differences, fill = Scenarios)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(alpha = 0.5, position = "identity") +
  labs(y = "Count",
       x = "Estimated - 'true'") +
  theme(legend.position = "none")


ps_k = ps_k_lm + inset_element(ps_k_hist, left = 0.07, bottom = 0.6, right = 0.5, top =  0.98)

ps_lambda_lm = ggplot(filter(all_raw_ps, method == "lambda"), 
                      aes(x = true, y = Phylo_signal, color = Scenarios, group = Scenarios)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_smooth(method = "lm") +
  labs(x = "'True' phylogenetic signal value",
       y = "Phylogenetic signal value (Lambda)",
       title = "(E) Lambda") +
  cowplot::theme_cowplot() +
  theme(legend.position = "none") 

ps_lambda_hist = dplyr::select(x_ps_1_1000_df, -true) |> 
  pivot_longer(cols = 2:4, names_to = "Scenarios", values_to = "Differences") |> 
  filter(method == "lambda") |> 
  ggplot(aes(x = Differences, fill = Scenarios)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(alpha = 0.5, position = "identity") +
  labs(y = "Count",
       x = "Estimated - 'true'") +
  theme(legend.position = "none")

ps_lambda = ps_lambda_lm + inset_element(ps_lambda_hist, left = 0.4, bottom = 0.01, right = 0.83, top =  0.39)

ps_all = ps_k + ps_lambda

pd_ps_all = pd_all / ps_all

ggsave("Figs/pd_ps_scenarios.pdf", plot = pd_ps_all, width = 14, height = 9)
ggsave("Figs/pd_ps_scenarios.png", plot = pd_ps_all, width = 14, height = 9)


# Plot for publication (2 scenarios) ====
pd_faith_lm = ggplot(filter(all_raw, diversity == "Faith's PD", Scenarios != "at_or_above_basal"), 
                     aes(x = true, y = Diversity, color = Scenarios, group = Scenarios)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_smooth(method = "lm") +
  labs(x = "'True' diversity value",
       y = "Estimated diversity value",
       title = "(A) Faith's PD") +
  cowplot::theme_cowplot() +
  theme(legend.position = "none") 

pd_faith_hist = ggplot(filter(all_diff, diversity == "Faith's PD", Scenarios != "at_or_above_basal"), 
                       aes(x = Differences, fill = Scenarios)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(alpha = 0.5, position = "identity") +
  labs(y = "Count",
       x = "Estimated - 'true'") +
  theme(legend.position = "none")

pd_faith = pd_faith_lm + inset_element(pd_faith_hist, left = 0.51, bottom = 0.01, right = 0.9, top =  0.38)


pd_mpd_lm = ggplot(filter(all_raw, diversity == "MPD", Scenarios != "at_or_above_basal"), 
                   aes(x = true, y = Diversity, color = Scenarios, group = Scenarios)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_smooth(method = "lm") +
  labs(x = "'True' diversity value",
       y = "Estimated diversity value",
       title = "(B) MPD") +
  cowplot::theme_cowplot() +
  theme(legend.position = c(0.1, 0.8)) 

pd_mpd_hist = ggplot(filter(all_diff, diversity == "MPD", Scenarios != "at_or_above_basal"), 
                     aes(x = Differences, fill = Scenarios)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(alpha = 0.5, position = "identity") +
  labs(y = "Count",
       x = "Estimated - 'true'") +
  theme(legend.position = "none")

pd_mpd = pd_mpd_lm + inset_element(pd_mpd_hist, left = 0.51, bottom = 0.01, right = 0.9, top =  0.38)

pd_mntd_lm = ggplot(filter(all_raw, diversity == "MNTD", Scenarios != "at_or_above_basal"), 
                    aes(x = true, y = Diversity, color = Scenarios, group = Scenarios)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_smooth(method = "lm") +
  labs(x = "'True' diversity value",
       y = "Estimated diversity value",
       title = "(C) MNTD") +
  cowplot::theme_cowplot() +
  theme(legend.position = "none") 

pd_mntd_hist = ggplot(filter(all_diff, diversity == "MNTD", Scenarios != "at_or_above_basal"), 
                      aes(x = Differences, fill = Scenarios)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(alpha = 0.5, position = "identity") +
  labs(y = "Count",
       x = "Estimated - 'true'") +
  theme(legend.position = "none")

pd_mntd = pd_mntd_lm + inset_element(pd_mntd_hist, left = 0.51, bottom = 0.01, right = 0.9, top =  0.38)

pd_all = pd_faith + pd_mpd + pd_mntd

# phylo signal
ps_k_lm = ggplot(filter(all_raw_ps, method == "K", Scenarios != "at_or_above_basal"), 
                 aes(x = true, y = Phylo_signal, color = Scenarios, group = Scenarios)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_smooth(method = "lm") +
  labs(x = "'True' phylogenetic signal value",
       y = "Phylogenetic signal value (K)",
       title = "(D) K") +
  cowplot::theme_cowplot() +
  theme(legend.position = c(0.1, 0.7)) 

ps_k_hist = dplyr::select(x_ps_1_1000_df, -true) |> 
  pivot_longer(cols = 2:4, names_to = "Scenarios", values_to = "Differences") |> 
  filter(method == "K", Scenarios != "diff_at_or_above_basal_true") |> 
  ggplot(aes(x = Differences, fill = Scenarios)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(alpha = 0.5, position = "identity") +
  labs(y = "Count",
       x = "Estimated - 'true'") +
  theme(legend.position = "none")


ps_k = ps_k_lm + inset_element(ps_k_hist, left = 0.56, bottom = 0.01, right = 0.99, top =  0.39)
# ps_k = ps_k_lm + inset_element(ps_k_hist, left = 0.56, bottom = 0.01, right = 0.95, top =  0.38)

ps_lambda_lm = ggplot(filter(all_raw_ps, method == "lambda", Scenarios != "at_or_above_basal"), 
                      aes(x = true, y = Phylo_signal, color = Scenarios, group = Scenarios)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_smooth(method = "lm") +
  labs(x = "'True' phylogenetic signal value",
       y = "Phylogenetic signal value (Lambda)",
       title = "(E) Lambda") +
  cowplot::theme_cowplot() +
  theme(legend.position = "none") 

ps_lambda_hist = dplyr::select(x_ps_1_1000_df, -true) |> 
  pivot_longer(cols = 2:4, names_to = "Scenarios", values_to = "Differences") |> 
  filter(method == "lambda", Scenarios != "diff_at_or_above_basal_true") |> 
  ggplot(aes(x = Differences, fill = Scenarios)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(alpha = 0.5, position = "identity") +
  labs(y = "Count",
       x = "Estimated - 'true'") +
  theme(legend.position = "none")

ps_lambda = ps_lambda_lm + inset_element(ps_lambda_hist, left = 0.4, bottom = 0.01, right = 0.83, top =  0.39)

ps_all = ps_k + ps_lambda

pd_ps_all = pd_all / ps_all

ggsave("Figs/fig_2.pdf", plot = pd_ps_all, width = 14, height = 9)
