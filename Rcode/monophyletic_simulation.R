library(tidyverse)
library(ape)
library(patchwork)

remotes::install_github("daijiang/rtrees", ref = "largest_cluster") # branch


test_monophy = function(seed = 1, n_sp = 500, n_bind = 200, 
                        tr = megatrees::tree_fish_12k, 
                        taxon = "fish"){
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
  
  # generate phylogenies
  t_inclusive = rtrees::get_tree(sp_list, tree = tr2, taxon = taxon, non_monophyletic = "inclusive")
  
  t_largest = rtrees::get_tree(sp_list, tree = tr2, taxon = taxon, non_monophyletic = "largest_cluster")
  
  tr_sub2 = ape::keep.tip(tr_sub, t_inclusive$tip.label)
  
  # communities
  comm = as.data.frame(matrix(1, nrow = 1, ncol = ape::Ntip(t_inclusive)))
  names(comm) = t_inclusive$tip.label
  
  pd = tibble::tibble(type = c("true", "inclusive", "largest_cluster"),
                      pd = c(picante::pd(comm, tr_sub2)$PD,
                             picante::pd(comm, t_inclusive)$PD,
                             picante::pd(comm, t_largest)$PD),
                      mpd = c(picante::mpd(comm, cophenetic(tr_sub2)),
                              picante::mpd(comm, cophenetic(t_inclusive)),
                              picante::mpd(comm, cophenetic(t_largest))),
                      mntd = c(  picante::mntd(comm, cophenetic(tr_sub2)),
                                 picante::mntd(comm, cophenetic(t_inclusive)),
                                 picante::mntd(comm, cophenetic(t_largest))))
  
  # phylo signals
  # phy signal
  trait_sim_BM = ape::rTraitCont(tr_sub2, model = "BM", sigma = 1)
  
  get_phy_signal = function(v, phy, v.name){
    if(is.null(names(v))) names(v) = v.name
    xy_k = phytools::phylosig(tree = phy, x = v, method = "K", test = FALSE)
    xy_l = phytools::phylosig(tree = phy, x = v, method = "lambda", test = FALSE)
    tibble::tibble(k = xy_k, lambda = xy_l$lambda)
  }
  
  ps_true = get_phy_signal(v = trait_sim_BM, phy = tr_sub2)
  ps_inclusive = get_phy_signal(v = trait_sim_BM, phy = t_inclusive)
  ps_largest = get_phy_signal(v = trait_sim_BM, phy = t_largest)
  ps = dplyr::bind_rows(
    dplyr::mutate(ps_true, type = "true"),
    dplyr::mutate(ps_inclusive, type = "inclusive"),
    dplyr::mutate(ps_largest, type = "largest_cluster")
  )
  
  dplyr::left_join(pd, ps, by = "type")
}

x_monphy_1_1000 = vector("list", 1000)
for(i in 2:1000){
  cat("i =", i, "\t", file = "prog.txt")
  x_monphy_1_1000[[i]] = try(test_monophy(i))
  if(i %% 10 == 0) 
    saveRDS(x_monphy_1_1000, file = "Data/x_monphy_1_1000.rds")
}
file.remove("prog.txt")

x_monphy_1_1000[[634]] # error in phylo signal calculation

# x_monophy = purrr::map_dfr(1:500, test_monophy, .id = "iterations")
# saveRDS(x_monophy, "Data/x_monophy.rds")
# 
# pivot_longer(x_monophy, cols = 3:5, names_to = "diversity") |> 
#   pivot_wider(names_from = "type", values_from = "value") |> 
#   ggplot(aes(x = inclusive, y = largest_cluster)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   geom_abline(intercept = 0, slope = 1) +
#   facet_wrap(~diversity, scales = "free") +
#   labs(x = "Inclusive MRCA for non-monophyletic groups",
#        y = "MRCA of the largest cluster within non-monophyletic groups")
# ggsave(filename = "Figs/monophyletic.pdf", width = 12, height = 5)

x = readRDS("x_monphy_1_1000.rds")
map(x, class)
x2 = x[-634]
x3 = bind_rows(x2, .id = "iteration")
x3$k = as.numeric(x3$k)
x4 = pivot_longer(x3, 3:7) |> 
  pivot_wider(names_from = "type", values_from = "value") |> 
  pivot_longer(cols = c("inclusive", "largest_cluster"), names_to = "monophy")

ggplot(x4, aes(x = true, y = value, color = monophy)) +
  geom_point(alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_smooth(method = "lm") +
  facet_wrap(~name, scales = "free")

diff_monophy = pivot_longer(x3, 3:7) |> 
  pivot_wider(names_from = "type", values_from = "value") |> 
  mutate(diff_inclusive_largest = inclusive - largest_cluster,
         diff_inclusive_true = inclusive - true,
         diff_largest_true = largest_cluster - true) |> 
  dplyr::select(name, diff_inclusive_true, diff_largest_true) |> 
  set_names(c("name", "inclusive", "largest_cluster")) |> 
  pivot_longer(cols = 2:3, names_to = "monophy")

pd_faith_lm = ggplot(filter(x4, name == "pd"), 
                     aes(x = true, y = value, color = monophy, group = monophy)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_smooth(method = "lm") +
  labs(x = "'True' diversity value",
       y = "Estimated diversity value",
       title = "(A) Faith's PD") +
  cowplot::theme_cowplot() +
  theme(legend.position = "none") 

pd_faith_hist = ggplot(filter(diff_monophy, name == "pd"), 
                       aes(x = value, fill = monophy)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(alpha = 0.5, position = "identity") +
  labs(y = "Count",
       x = "Estimated - 'true'") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45))

pd_faith = pd_faith_lm + inset_element(pd_faith_hist, left = 0.01, bottom = 0.6, right = 0.4, top =  0.98)

pd_mpd_lm = ggplot(filter(x4, name == "mpd"), 
                   aes(x = true, y = value, color = monophy, group = monophy)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_smooth(method = "lm") +
  labs(x = "'True' diversity value",
       y = "Estimated diversity value",
       title = "(B) MPD",
       color = "Non-monophyletic") +
  cowplot::theme_cowplot() +
  theme(legend.position = c(0.6, 0.3)) 

pd_mpd_hist = ggplot(filter(diff_monophy, name == "mpd"), 
                     aes(x = value, fill = monophy)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(alpha = 0.5, position = "identity") +
  labs(y = "Count",
       x = "Estimated - 'true'") +
  theme(legend.position = "none")

pd_mpd = pd_mpd_lm + inset_element(pd_mpd_hist, left = 0.01, bottom = 0.6, right = 0.4, top =  0.98)


pd_mntd_lm = ggplot(filter(x4, name == "mntd"), 
                    aes(x = true, y = value, color = monophy, group = monophy)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_smooth(method = "lm") +
  labs(x = "'True' diversity value",
       y = "Estimated diversity value",
       title = "(C) MNTD") +
  cowplot::theme_cowplot() +
  theme(legend.position = "none") 

pd_mntd_hist = ggplot(filter(diff_monophy, name == "mntd"), 
                      aes(x = value, fill = monophy)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(alpha = 0.5, position = "identity") +
  labs(y = "Count",
       x = "Estimated - 'true'") +
  theme(legend.position = "none")

pd_mntd = pd_mntd_lm + inset_element(pd_mntd_hist, left = 0.01, bottom = 0.6, right = 0.4, top =  0.99)

pd_all = pd_faith + pd_mpd + pd_mntd

ps_k_lm = ggplot(filter(x4, name == "k"), 
                 aes(x = true, y = value, color = monophy, group = monophy)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_smooth(method = "lm") +
  labs(x = "'True' phylogenetic signal",
       y = "Estimated phylogenetic signal",
       title = "(D) K") +
  cowplot::theme_cowplot() +
  theme(legend.position = "none") 

ps_k_hist = ggplot(filter(diff_monophy, name == "k", value > -2), 
                   aes(x = value, fill = monophy)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(alpha = 0.5, position = "identity") +
  labs(y = "Count",
       x = "Estimated - 'true'") +
  theme(legend.position = "none")

ggplot(filter(diff_monophy, name == "k", value > -2), 
       aes(x = value, fill = monophy)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(alpha = 0.5, position = "identity") +
  facet_wrap(~monophy)

ps_k = ps_k_lm + inset_element(ps_k_hist, left = 0.01, bottom = 0.6, right = 0.4, top =  0.98)

ps_lambda_lm = ggplot(filter(x4, name == "lambda"), 
                      aes(x = true, y = value, color = monophy, group = monophy)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_smooth(method = "lm") +
  labs(x = "'True' phylogenetic signal",
       y = "Estimated phylogenetic signal",
       title = "(E) Lambda") +
  cowplot::theme_cowplot() +
  theme(legend.position = "none") 

ps_lambda_hist = ggplot(filter(diff_monophy, name == "lambda"), 
                        aes(x = value, fill = monophy)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(alpha = 0.5, position = "identity") +
  labs(y = "Count",
       x = "Estimated - 'true'") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45))

ggplot(filter(diff_monophy, name == "lambda"), 
       aes(x = value, fill = monophy)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_histogram(alpha = 0.5, position = "identity") +
  facet_wrap(~monophy)

ps_lambda = ps_lambda_lm + inset_element(ps_lambda_hist, left = 0.01, bottom = 0.01, right = 0.4, top =  0.49)

ps_all = ps_k + ps_lambda

pd_ps_all = pd_all / ps_all

ggsave("Figs/monophyletic.pdf", plot = pd_ps_all, width = 14, height = 9)
ggsave("Figs/monophyletic.png", plot = pd_ps_all, width = 14, height = 9)
