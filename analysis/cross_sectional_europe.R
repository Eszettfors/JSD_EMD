library(MASS)
library(tidyverse)
library(rnaturalearth)
library(sf)
library(treemapify)
library(xtable)
library(scales)
library(vegan)
library(factoextra)
library(ggrepel)
library(plotly)
library(ape)
library(GGally)
library(pheatmap)
library(xtable)
library(DescTools)
library(dendextend)
library(patchwork)
library(plotly)
library(randomcoloR)

set.seed(42)

# read data
df_div = read_csv("data/processed/diversity_cs.csv")
df_pw = read_csv("data/processed/page_views_cs.csv")
df_pairs = read_csv("data/processed/pairwise_distances.csv")

# resolve problem with Namibia
df_div = df_div %>%
  mutate(country_code = ifelse(is.na(country_code), "NA", country_code))

df_pw = df_pw %>%
  mutate(country_code = ifelse(is.na(country_code), "NA", country_code))

df_pairs = df_pairs %>%
  mutate(country_1 = ifelse(is.na(country_1), "NA", country_1),
         country_2 = ifelse(is.na(country_2), "NA", country_2))
  
# adjust theme
set_theme(theme_bw(base_size = 10))

# european countries
europe_countries = df_div %>%
  filter(continent == "Europe") %>%
  pull(country_code)

# add Greenland and Turkey, cyrpus, armenia, georgia, azerbaijan
europe_countries = c(europe_countries, "TR", "CY", "AM", "AZ", "GE")

# filter
df_div = df_div %>%
  filter(country_code %in% europe_countries)

df_pw = df_pw %>%
  filter(country_code %in% europe_countries)

df_pairs = df_pairs %>%
  filter(country_1 %in% europe_countries, country_2 %in% europe_countries)

# prepare geo data
df_geo = ne_countries(scale = "medium")
df_geo = df_geo %>%
  select(iso_a2_eh, name, pop_est, subregion) %>%
  rename("country_code" = iso_a2_eh)

df_geo = st_transform(df_geo, 3035)


# function for calculating JSD
get_jsd = function(p, q){
  p = p / sum(p) + 0.00000001
  q = q / sum(q) + 0.00000001
  
  M = 0.5 * (p + q)
  
  jsd = 0.5 * sum(p * log2(p / M) + q * log2(q / M))
  return(sqrt(jsd))
}

# some tests
A = c(0.25, 0.5, 0.25, 0, 0)
B = c(0.5, 0, 0, 0.1, 0.4)
get_jsd(A, B)

A2 = c(0.5, 0.5, 0,0,0)
get_jsd(A2, B)


get_jsd(c(0,1), c(1,0))

# Descriptive statistics
df_pw = df_pw %>%
  group_by(country_name) %>%
  mutate(prop = page_views / sum(page_views))

df_pw %>%
  group_by(country_name) %>%
  summarize(total_page_views = sum(page_views)) %>%
  ggplot(aes(y = reorder(country_name, total_page_views), x = total_page_views)) + 
  geom_bar(stat = "identity", fill = "lightblue", colour = "black") + 
  labs(y = "", x = "Total number of Wikipedia page views (2015-2024)")

compare_prob_dist = function(langs, prob_lim = 0.01, cols = c("forestgreen", "lightblue")){
  df_pw %>%
    filter(country_name %in% langs) %>%
    filter(prop > prob_lim) %>%
    ggplot(aes(y = prop, x = reorder(language_name, prop), fill = country_name)) + 
    geom_bar(stat = "identity", position = "dodge", colour = "black") + 
    geom_text(aes(label = paste0(as.character(round(prop,2)))), position = position_dodge(width = 0.9), vjust = -0.2) +
    labs(y = "Probability", x = "",
         fill = "Country") + 
    theme_bw(base_size = 11) + 
    theme(axis.text.x = element_text(angle = 15),
          legend.position = "inside",
          legend.position.inside = c(0.2, 0.6),
          legend.text = element_text(size = 13)) + 
    ylim(c(0,1)) + 
    scale_fill_manual(values = cols)
}

vatican_vs_italy = compare_prob_dist(c("Vatican City", "Italy"), cols = c("forestgreen", "yellow"))
swiss_vs_lux = compare_prob_dist(c("Switzerland", "Luxembourg"), cols = c("lightblue", "red"))
serb_vs_croat = compare_prob_dist(c("Serbia", "Croatia"), cols = c("white", muted("blue")))
compare_prob_dist(c("Serbia", "Croatia", "Bosnia and Herzegovina"), cols = c("red", "white", "blue"))
prob_plots = vatican_vs_italy / serb_vs_croat
ggsave("analysis/plots/prob_plots.jpeg", prob_plots, dpi = 300, width = 8, height = 6)

# pairwise distances
Desc(df_pairs$jsd)

Desc(df_pairs$jsd_phyl_dist)

Desc(df_pairs$naive_div_dist)

Desc(df_pairs$phyl_div_dist)

# map
naive_div = df_div %>%
  right_join(df_geo, join_by(country_code)) %>%
  st_as_sf() %>%
  ggplot(aes(geometry = geometry, fill = naive_div_q_1)) + 
  geom_sf(colour = "black", size = 0.2) + 
  scale_fill_gradient(low = ("lightblue"), high = muted("blue")) + 
  labs(fill = "Eff. no. \nLanguages",
       title = "A) Naive diversity, q = 1") +
  theme(legend.position = c(0.1, 0.5),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.background = element_blank()) +
  coord_sf(xlim = c(2500000, 7300000), ylim = c(1500000, 5200000))


phyl_div = df_div %>%
  right_join(df_geo, join_by(country_code)) %>%
  st_as_sf() %>%
  ggplot(aes(geometry = geometry, fill = phyl_div_q_1)) + 
  geom_sf(colour = "black", size = 0.2) + 
  scale_fill_gradient(low = ("lightblue"), high = muted("blue")) + 
  labs(fill = "Eff. no. \nLanguages",
       title = "B) Phylogenetic diversity, q = 1") +
  theme(legend.position = c(0.1, 0.5),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.background = element_blank()) +
  coord_sf(xlim = c(2500000, 7300000), ylim = c(1500000, 5200000))

# dominant language
df_dom_lang = df_pw %>%
  group_by(country_code) %>%
  slice_max(n = 1, prop)

dom_lang = df_dom_lang %>%
  right_join(df_geo, join_by(country_code)) %>%
  st_as_sf() %>%
  ggplot(aes(geometry = geometry, fill = language_name)) + 
  geom_sf(colour = "black", size = 0.2) + 
  labs(fill = element_blank(),
       title = "C) Most Widely Used Wikipedia Language") +
  theme(legend.position = c(0.85, 0.65),
        axis.ticks = element_blank(),
        legend.key.size = unit(0.4, "line"),
        axis.text = element_blank(),
        legend.background = element_blank()) + 
  coord_sf(xlim = c(2500000, 7300000), ylim = c(1500000, 5200000))

dom_2_lang = df_pw %>%
  group_by(country_code) %>%
  slice_max(n = 2, prop) %>%
  slice_min(n = 1, prop) %>%
  right_join(df_geo, join_by(country_code)) %>%
  st_as_sf() %>%
  ggplot(aes(geometry = geometry, fill = language_name)) + 
  geom_sf(colour = "black", size = 0.2) + 
  labs(fill = element_blank(),
       title = "D) Second Most Widely Used Wikipedia Language") +
  theme(legend.position = c(0.80, 0.65),
        axis.ticks = element_blank(),
        legend.key.size = unit(0.4, "line"),
        axis.text = element_blank(),
        legend.background = element_blank()) + 
  coord_sf(xlim = c(2500000, 7300000), ylim = c(1500000, 5200000))


# make sure colorus are the same across plots
all_languages = union(dom_lang$data$language_name, dom_2_lang$data$language_name)
colors = c(Spanish = "#F1BF00", English = "#CE1124", Armenian = "purple", German = "grey", Swedish = "#006AA7", French = "#000091", Bulgarian = "#00966E",
           Azerbaijani = "#0092BC", Italian = "#CDDFA0", Russian = muted("purple"), Georgian = "#BA274A", Czech = "#BA9D27", Slovak = "#B2ECE1", Finnish = "#002F6C", 
           Turkish = "pink", Polish = "#EE85B5", Albanian = "#BA0C0C", Croatian = "#6BAA75", Icelandic = "#02529C", Lithuanian = "#046A38", Romanian = "#FCD116",
           Portuguese = "#55868C", Macedonian = "#253133", Ukrainian = "lightblue", Slovenian = "#A8DCD1", Serbian = "#73BFB8", Dutch = "orange", Hungarian = "#77AF9C",
           Norwegian = "#EFF1F3", Estonian = "#0072CE", Danish = "#ED6A5A", Greek = "#0D5EAF")

color_palette = colors[all_languages]

dom_lang = dom_lang +
  scale_fill_manual(values = color_palette, name = element_blank())
print(dom_lang)
dom_2_lang = dom_2_lang +
  scale_fill_manual(values = color_palette, name = element_blank())
print(dom_2_lang)

# Display the plots
print(dom_lang)
print(dom_2_lang)

div_maps = (naive_div + phyl_div) / (dom_lang + dom_2_lang)

ggsave("analysis/plots/diversity_maps_europe.png", div_maps, width = 12 , height = 8, dpi = 300)


# plot diversity against one another
diversity_cor = cor.test(df_div$naive_div_q_1, df_div$phyl_div_q_1, method = "spearman")

phyl_vs_naive_div = df_div %>%
  left_join(df_geo) %>%
  mutate(subregion = ifelse(country_name == "Gibraltar", "Southern Europe", subregion)) %>%
  ggplot(aes(y = rank(phyl_div_q_1), x = rank(naive_div_q_1), fill = subregion)) + 
  geom_point() + 
  geom_label_repel(aes(label = country_name)) + 
  geom_abline(intercept = 0, slope = 1) +
  labs(y = "Phylogenetic diversity, q = 1 (rank)",
       x = "Naive diversity, q = 1 (rank)",
       title = paste0("Relationship between phylogenetic and naive diveristy, q = 1, Sperman's Rho = ", round(diversity_cor$estimate, 2))) + 
  theme(legend.position = c(0.9,0.2))

ggsave("analysis/plots/phyl_vs_naive_div.png", phyl_vs_naive_div, width = 8, height = 8, dpi = 300)

# map difference in diversity
df_div %>%
  right_join(df_geo, join_by(country_code)) %>%
  mutate(div_diff = phyl_div_q_1 - naive_div_q_1) %>%
  st_as_sf() %>%
  ggplot(aes(geometry = geometry, fill = div_diff)) + 
  geom_sf(colour = "black", size = 0.2) + 
  scale_fill_gradient2(low = ("lightblue"), mid = "white") + 
  labs(fill = "Eff. no. \nLanguages",
       title = "A) diveristy difference") +
  theme(legend.position = c(0.1, 0.5),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.background = element_blank()) +
  coord_sf(xlim = c(2500000, 7300000), ylim = c(1500000, 5200000))


# table of countries according to diversity
df_div %>%
  select(country_name, naive_div_q_1, phyl_div_q_1) %>%
  mutate(rank_naive = rank(-naive_div_q_1),
         rank_phyl = rank(-phyl_div_q_1)) %>%
  arrange(-naive_div_q_1)

# treemap
get_treemap = function(country = "Sweden"){
  df_pw %>%
    filter(country_name %in% country) %>%
    ggplot(aes(area = prop, fill = language_name)) + 
    geom_treemap(col = "black") + 
    geom_treemap_text(aes(label = paste0(language_name, "\n", round(prop * 100, 2), "%")),
                      grow = TRUE, place = "center",
                      colour = "white") + 
    theme(legend.position = "none",
          strip.text = element_text(size = 18),
          axis.text.x = element_blank()) + 
    facet_wrap(~country_name)
}

get_treemap(c("Vatican City", "Italy"))

# similarities ####
jsd_mat = read_rds("data/processed/jsd_mat.rds")
phyl_jsd_mat = read_rds("data/processed/phyl_jsd_dist_mat.rds")
naive_div_dist_mat = read_rds("data/processed/shannon_dist_mat.rds")
phyl_div_dist_mat = read_rds("data/processed/phyl_div_dist_mat.rds")

length(europe_countries)

# filter to europe
jsd_mat = jsd_mat[europe_countries, europe_countries]
phyl_jsd_mat = phyl_jsd_mat[europe_countries, europe_countries]
naive_div_dist_mat = naive_div_dist_mat[europe_countries, europe_countries]
phyl_div_dist_mat = phyl_div_dist_mat[europe_countries, europe_countries]

#### correlate the distance matrices #####
ggpairs(df_pairs %>%
          select(!c(country_1, country_2)))

hex_fn = function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_hex(bins = 30, ...) +
    scale_fill_viridis_c() +
    theme_minimal(base_size = 14) + 
    geom_smooth(method = "lm")
}

corrs = ggpairs(
  data = df_pairs %>% select(!c(country_1, country_2)) %>%
    rename("JSD" = jsd, 
           "Phylo. EMD" = jsd_phyl_dist,
           "Naive Diveristy Distance" = naive_div_dist,
           "Phylo. Diveristy Distance" = phyl_div_dist),
  lower = list(continuous = wrap(hex_fn)),
  upper = list(continuous = wrap("cor", size = 8)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.5, fill = "lightblue")),
) + 
  theme_bw(base_size = 14) + 
  theme(
    strip.text = element_text(size = 10),
    panel.grid = element_blank()
  )
print(corrs)

ggsave("analysis/plots/corrs.png", corrs, width = 10, height = 8, dpi = 300)

# most and least similar countries
df_pairs = df_pairs %>%
  left_join(df_div %>%
              select(country_code, country_name), join_by("country_1" == "country_code")) %>%
  left_join(df_div %>%
              select(country_code, country_name), join_by("country_2" == "country_code")) %>%
  rename(country_name_1 = country_name.x,
         country_name_2 = country_name.y)

df_pairs = df_pairs %>%
  mutate(diff_jsd = jsd - jsd_phyl_dist)

Desc(df_pairs$jsd) # mean dist 0.70
Desc(df_pairs$jsd_phyl_dist) # mean dist 0.54
Desc(df_pairs$diff_jsd) # mean 0.16
  

#### clustering #####

#### Heatmaps ####
# jsd
country_names = df_div$country_name
country_codes = df_div$country_code

lookup = setNames(df_div$country_name, df_div$country_code)

rownames(jsd_mat) = lookup[rownames(jsd_mat)]
colnames(jsd_mat) = lookup[colnames(jsd_mat)]

heatmap_jsd = pheatmap::pheatmap(jsd_mat, angle_col = 45,
                   fontsize_row = 6,
                   fontsize_col = 6,
                   main = "Heatmap of pairwise JSD, hierarchical clustering") 
ggsave("analysis/plots/heatmap_jsd.png", heatmap_jsd, width = 12, height = 8, dpi = 300)

# EMD
rownames(phyl_jsd_mat) = lookup[rownames(phyl_jsd_mat)]
colnames(phyl_jsd_mat) = lookup[colnames(phyl_jsd_mat)]

heatmap_phyl = pheatmap::pheatmap(phyl_jsd_mat, angle_col = 45,
                                  fontsize_row = 6,
                                  fontsize_col = 6, main = "Heatmap of pairwise Phylogenetic EMD, hierarchical clustering")
ggsave("analysis/plots/heatmap_phyl.png", heatmap_phyl, width = 12, height = 8, dpi = 300)

##### hierarchical clustering #####

## jsd
d_jsd = as.dist(jsd_mat)
hc_jsd_single = hclust(d_jsd, method = "single")
hc_jsd_average = hclust(d_jsd, method = "average")
hc_jsd_complete = hclust(d_jsd, method = "complete")
hc_jsd_ward = hclust(d_jsd, method = "ward.D2")

c_jsd_single = cophenetic(hc_jsd_single)
c_jsd_average = cophenetic(hc_jsd_average)
c_jsd_complete = cophenetic(hc_jsd_complete)
c_jsd_ward = cophenetic(hc_jsd_ward)

cor.test(d_jsd, c_jsd_single)
cor.test(d_jsd, c_jsd_average)
cor.test(d_jsd, c_jsd_complete)
cor.test(d_jsd, c_jsd_ward)

plot(hc_jsd_average)

fviz_nbclust(jsd_mat, FUNcluster = hcut, method = "silhouette", k.max = 30) + 
  geom_vline(xintercept = c(7, 12, 18), colour = muted("blue"), linetype = "dashed")

get_dend_k = function(dist_mat, k, method = "average"){
  
  hc =  hclust(as.dist(dist_mat), method = method)
  
  p = fviz_dend(
    hc,
    k = k,
    rect = TRUE,
    rect_border = "jco",
    rect_fill = TRUE,
    cex = 0.7,
    lwd = 0.8,
    color_labels_by_k = TRUE,
    horiz = TRUE
  ) +
    theme_minimal(base_size = 14)
  return(p)
}

get_dend_k(jsd_mat, k = 7, method = "average") +
  labs(title = "Hierarchical Clustering, JSD, k = 7",
       y = "Jensen-Shannon Distance")
dendrogram_jsd_12 = get_dend_k(jsd_mat, k = 12, method = "average") +
  labs(title = "Hierarchical Clustering, JSD, k = 12",
       y = "Jensen-Shannon Distance")
get_dend_k(jsd_mat, k = 18, method = "average") + 
  labs(title = "Hierarchical Clustering, JSD, k = 18",
       y = "Jensen-Shannon Distance")
get_dend_k(jsd_mat, k = 24, method = "average") + 
  labs(title = "Hierarchical Clustering, JSD, k = 24",
       y = "Jensen-Shannon Distance")


ggsave("analysis/plots/dendrogram.png",dendrogram_jsd_12, height = 12, width = 6, dpi = 300)


### emd
d_phyl = as.dist(phyl_jsd_mat)
hc_phyl = hclust(d_phyl, method = "average")

hc_phyl_single = hclust(d_phyl, method = "single")
hc_phyl_average = hclust(d_phyl, method = "average")
hc_phyl_complete = hclust(d_phyl, method = "complete")
hc_phyl_ward = hclust(d_phyl, method = "ward.D2")

c_phyl_single = cophenetic(hc_jsd_single)
c_phyl_average = cophenetic(hc_jsd_average)
c_phyl_complete = cophenetic(hc_jsd_complete)
c_phyl_ward = cophenetic(hc_jsd_ward)

cor.test(d_phyl, c_phyl_single)
cor.test(d_phyl, c_phyl_average)
cor.test(d_phyl, c_phyl_complete)
cor.test(d_phyl, c_phyl_ward)


fviz_nbclust(phyl_jsd_mat, FUNcluster = hcut, method = "silhouette", k.max = 30) + 
  geom_vline(xintercept = c(7, 21, 26), colour = muted("blue"), linetype = "dashed")

get_dend_k(phyl_jsd_mat, k = 7, method = "average") + 
  labs(title = "Hierarchical Clustering, Phylogenetic EMD, k = 7",
       y = "Phylogenetic EMD")
get_dend_k(phyl_jsd_mat, k = 13, method = "average")  + 
  labs(title = "Hierarchical Clustering, Phylogenetic EMD, k = 13",
       y = "Phylogenetic EMD")
get_dend_k(phyl_jsd_mat, k = 21, method = "average") + 
  labs(title = "Hierarchical Clustering, Phylogenetic EMD, k = 21",
       y = "Phylogenetic EMD")
get_dend_k(phyl_jsd_mat, k = 26, method = "average")  + 
  labs(title = "Hierarchical Clustering, Phylogenetic EMD, k = 26",
       y = "Phylogenetic EMD")


dendrogram_phyl_13 = fviz_dend(
  hc_phyl,
  k = 13,
  rect = TRUE,
  rect_border = "jco",
  rect_fill = TRUE,
  cex = 0.7,
  lwd = 0.8,
  color_labels_by_k = TRUE,
  horiz = TRUE
) +
  labs(
    title = "B) Hierarchical Clustering, Phylogenetic EMD",
    y = "Phylogenetic EMD"
  ) +
  theme_minimal(base_size = 14) + 
  theme(
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),  
    panel.grid.minor.y = element_blank())
print(dendrogram_phyl_13)

ggsave("analysis/plots/dendrogram_phyl.png", dendrogram_phyl_13, height = 12, width = 6, dpi = 300)


cor_cophenetic(hc_jsd_average, hc_phyl_average) # r = 0.844


## population size vs diversity
pop_div = df_div %>%
  left_join(df_dom_lang) %>%
  left_join(df_geo) %>%
  rename("dominant_language" = language_name)
  


## add estiamte for gibraltar
pop_div = pop_div %>%
  mutate(pop_est = case_when(country_code == "GI" ~ 35520,
                             TRUE ~ pop_est))


pop_div %>%
  ggplot(aes(y = rank(naive_div_q_1), x = rank(pop_est))) + 
  geom_point(aes(fill = dominant_language)) + 
  geom_label(aes(label = country_name, fill = dominant_language)) + 
  geom_smooth(method = "lm") + 
  labs(y = "Naive Diversity, q = 1, Rank",
       fill = "Dominant Language",
       x = "Estimated Population (2019), Rank")

cor.test(pop_div$naive_div_q_1, pop_div$pop_est, method = "spearman")


pop_div %>%
  ggplot(aes(y = rank(phyl_div_q_1), x = rank(pop_est))) + 
  geom_point(aes(fill = dominant_language)) + 
  geom_label(aes(label = country_name, fill = dominant_language)) + 
  geom_smooth(method = "lm") + 
  labs(y = "Phylogenetic Diversity, q = 1, Rank",
       fill = "Dominant Language",
       x = "Estimated Population (2019), Rank")

cor.test(pop_div$phyl_div_q_1, pop_div$pop_est, method = "spearman")



# without tiny countries sharing official language of 
micro_countries = c("Jersey", "Guernsey", "Malta", "Isle of Man", "Gibraltar",
                      "Vatican City", "Åland", "Andorra", "Monaco", "San Marino", "Liechtenstein")

pop_div_no_micro = pop_div %>%
  filter(!country_name %in% micro_countries)

pop_div_no_micro %>%
  ggplot(aes(y = rank(naive_div_q_1), x = rank(pop_est))) + 
  geom_point(aes(fill = dominant_language)) + 
  geom_label(aes(label = country_name, fill = dominant_language)) + 
  geom_smooth(method = "lm") + 
  labs(y = "Naive Diversity, q = 1, Rank",
       fill = "Dominant Language",
       x = "Estimated Population (2019), Rank")
  
cor.test(pop_div_no_micro$naive_div_q_1, pop_div_no_micro$pop_est, method = "spearman")

##
pop_div_no_micro %>%
  ggplot(aes(y = rank(phyl_div_q_1), x = rank(pop_est))) + 
  geom_point(aes(fill = dominant_language)) + 
  geom_label(aes(label = country_name, fill = dominant_language)) + 
  geom_smooth(method = "lm") + 
  labs(y = "Phylogenetic Diversity, q = 1, Rank",
       fill = "Dominant Language",
       x = "Estimated Population (2019), Rank")

cor.test(pop_div_no_micro$phyl_div_q_1, pop_div_no_micro$pop_est, method = "spearman")
