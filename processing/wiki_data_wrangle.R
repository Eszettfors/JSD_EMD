library(tidyverse)
library(vegan)
library(transport)

country_data = read_csv("data/raw/country_data.csv")
language_data = read_csv("data/raw/language_data.csv")
df_wiki = read_csv("data/raw/country_year_lang_views.csv")

# fix namibia
country_data = country_data %>%
  mutate(country_code = ifelse(is.na(country_code), "NA", country_code))

df_wiki = df_wiki %>%
  mutate(country_code = ifelse(is.na(country_code), "NA", country_code))

# make data complete, not present = 0 page views
df_wiki = df_wiki %>%
  complete(country_code, year, ISO6393) %>%
  mutate(page_views = ifelse(is.na(page_views), 0, page_views))

df_pw = country_data %>%
  full_join(df_wiki) %>%
  full_join(language_data)

df_pw = df_pw %>%
  mutate(country_code = ifelse(is.na(country_code), "NA", country_code))

write_csv(df_pw, "data/processed/page_views_ts.csv")

# cross sectional data

# summarize page views across all years
df_pw_cs = df_pw %>%
  group_by(country_code, country_name, ISO6393, language_name) %>%
  summarize(region_wb = unique(region_wb),
            continent = unique(continent),
            family = unique(family),
            page_views = sum(page_views)) %>%
  group_by(country_code) %>%
  mutate(prop = page_views / sum(page_views)) %>%
  ungroup() %>%
  arrange(country_code, ISO6393)


### phylogenetic distances ####
df_classification = language_data %>%
  select(ISO6393, language_name, classification) %>%
  arrange(ISO6393)

get_phyl_dist_matrix = function(vec_1, vec_2){
  # this function takes two vectors of ISO codes, joins them and calculates a phylogenetic similarity matrix for their union
  
  # define the union of the two vectors
  union_vec = union(vec_1, vec_2)
  
  # get the relevant ISO codes
  ISO = df_classification %>%
    filter(ISO6393 %in% union_vec) %>%
    select(ISO6393, classification)
  
  dist_mat = ISO %>%
    separate_rows(classification, sep = "/") %>%
    pivot_wider(names_from = classification,
                values_from = classification) %>%
    select(!ISO6393) %>%
    mutate(across(everything(), ~ case_when(!is.na(.) ~ 1,
                                            TRUE ~ 0))) %>%
    as.matrix() %>%
    vegdist(method = "jaccard") %>%
    as.matrix()
  
  
  colnames(dist_mat) = ISO$ISO6393
  rownames(dist_mat) = ISO$ISO6393
  
  return(dist_mat)
}
get_phyl_dist_matrix("rus", "lat")
phyl_dist_mat = get_phyl_dist_matrix(language_data$ISO6393, language_data$ISO6393)
phyl_sim_mat = 1 - phyl_dist_mat


# calculate diversity
get_shannon_entropy = function(pop_vec){
  # calculates the shannon entropy given a vector with populations
  
  # remove any zero
  pop_vec = pop_vec[pop_vec != 0]
  
  prop_vec = pop_vec / sum(pop_vec)
  
  # calculate entropy
  shannon = -1 * sum(prop_vec*log(prop_vec))
  return(shannon)
}


get_exp_shannon = function(pop_vec){
  # wrapper function to turn shannon entropy into a hill number
  exp_shannon = exp(get_shannon_entropy(pop_vec))
  return(exp_shannon)
}

get_shannon_diversity = function(prop, sim_m){
  # calculates diversity for q = 1 ergo shannon given a vector with proportions
  # remobve instances with 0
  
  idx_zero = prop != 0
  
  prop = prop[idx_zero]
  
  sim_m = sim_m[idx_zero, idx_zero]
  
  # for each proportion, get the expected similarity to all other proportions
  expected = log(sim_m %*% prop)
  
  # for each proportion, multiply by expected similarity to all other proportions
  # and derive entropy
  E = -1 * sum(prop * expected)
  
  # exponentiate entropy to get diversity
  D = exp(E)
  
  return(D)
}

df_diversity_cs = df_pw_cs %>%
  group_by(country_code, country_name) %>%
  summarize(continent = unique(continent), 
            region_wb = unique(region_wb),
            naive_div_q_1 = get_exp_shannon(pop_vec = prop),
            phyl_div_q_1 = get_shannon_diversity(prop = prop, sim_m = phyl_sim_mat))

# export
write_csv(df_pw_cs, "data/processed/page_views_cs.csv")
write_csv(df_diversity_cs, "data/processed/diversity_cs.csv")


# calculate distance matrices #####
df_pw_cs %>%
  count(country_code, ISO6393) %>%
  filter(n > 1)


# JSD ####
get_jsd = function(p, q) {
  p = p / sum(p) + 0.000001
  q = q / sum(q) + 0.000001
  
  
  M = 0.5 * (p + q)
  
  jsd = 0.5 * sum(p * log2(p / M) + q * log2(q / M))
  
  return(sqrt(jsd))
  
}

pw_wide = df_pw_cs %>%
  ungroup() %>%
  select(country_code, language_name, prop) %>%
  pivot_wider(names_from = language_name, values_from = prop, values_fill = 0) 

countries = pw_wide$country_code
jsd_mat = matrix(nrow = length(countries),
                 ncol = length(countries))
colnames(jsd_mat) = countries
rownames(jsd_mat) = countries

for (c1 in countries) {
  c1p = as.numeric(pw_wide[pw_wide$country_code == c1, 2:ncol(pw_wide)])
  for (c2 in countries){
    c2p = as.numeric(pw_wide[pw_wide$country_code == c2, 2:ncol(pw_wide)])
    JSD = get_jsd(c1p, c2p)
    
    jsd_mat[c1, c2] = JSD
  }
}
write_rds(jsd_mat, "data/processed/jsd_mat.rds")



# turn into pairwise similarities
df_jsd = data.frame(list(jsd_mat)[[1]])
df_jsd$country_1 = countries
df_jsd = df_jsd %>%
  pivot_longer(cols = !c(country_1), names_to = "country_2", values_to = "jsd") %>%
  filter(country_1 != country_2) %>%
  arrange(country_1) %>%
  filter(country_1 < country_2)

# export
write_rds(jsd_mat, "jsd_mat.rds")
write_csv(df_jsd, "pairwise_jsd")

# diversity distances
get_div_dist = function(x, y){
  div_dist = abs(x - y)
  return(div_dist)
}

# naive shannon
countries = df_diversity_cs$country_code
naive_div_dist_mat = matrix(nrow = length(countries),
                      ncol = length(countries))
colnames(naive_div_dist_mat) = countries
rownames(naive_div_dist_mat) = countries

for (c1 in countries) {
  div_c1 = df_diversity_cs %>%
    filter(country_code == c1) %>%
    pull(naive_div_q_1)
  for (c2 in countries){
    div_c2 = df_diversity_cs %>%
      filter(country_code == c2) %>%
      pull(naive_div_q_1)
    
    distance = get_div_dist(div_c1, div_c2)
    naive_div_dist_mat[c1, c2] = distance
  }
}

df_naive_div_dist = data.frame(list(naive_div_dist_mat)[[1]])
df_naive_div_dist$country_1 = countries
df_naive_div_dist = df_naive_div_dist %>%
  pivot_longer(cols = !c(country_1), names_to = "country_2", values_to = "naive_div_dist") %>%
  filter(country_1 != country_2) %>%
  arrange(country_1) %>%
  filter(country_1 < country_2)

# write
write_rds(naive_div_dist_mat, "data/processed/richness_dist_mat.rds")

# phyl diversity
phyl_div_dist_mat = matrix(nrow = length(countries),
                           ncol = length(countries))
colnames(phyl_div_dist_mat) = countries
rownames(phyl_div_dist_mat) = countries


for (c1 in countries) {
  div_c1 = df_diversity_cs %>%
    filter(country_code == c1) %>%
    pull(phyl_div_q_1)
  for (c2 in countries){
    div_c2 = df_diversity_cs %>%
      filter(country_code == c2) %>%
      pull(phyl_div_q_1)
    
    distance = get_div_dist(div_c1, div_c2)
    phyl_div_dist_mat[c1, c2] = distance
  }
}

df_phyl_div_dist = data.frame(list(phyl_div_dist_mat)[[1]])
df_phyl_div_dist$country_1 = countries
df_phyl_div_dist = df_phyl_div_dist %>%
  pivot_longer(cols = !c(country_1), names_to = "country_2", values_to = "phyl_div_dist") %>%
  filter(country_1 != country_2) %>%
  arrange(country_1) %>%
  filter(country_1 < country_2)

# export
write_rds(phyl_div_dist_mat, "data/processed/phyl_div_dist_mat.rds")


# phylogenetic wasserstein
phyl_dist_order = colnames(phyl_dist_mat)
pw_wide = df_pw_cs %>%
  ungroup() %>%
  arrange(ISO6393) %>%
  select(country_code, ISO6393, prop) %>%
  pivot_wider(names_from = ISO6393, values_from = prop, values_fill = 0)

order = colnames(pw_wide %>%
                   select(!country_code))
phyl_dist_mat = phyl_dist_mat[order, order]


get_phyl_jsd_dist = function(country_code_1, country_code_2){
  c1 = as.numeric(pw_wide[pw_wide$country_code == country_code_1, 2:ncol(pw_wide)])
  
  c2 = as.numeric(pw_wide[pw_wide$country_code == country_code_2, 2:ncol(pw_wide)])
  
  plan = transport::transport(
    c1,
    c2,
    costm = phyl_dist_mat,
    method = "shortsimplex"
  )
  dist = sum(plan$mass * phyl_dist_mat[cbind(plan$from, plan$to)])
  return(dist)
}

# loop and populate matrix
distance = c()
phyl_jsd_mat = matrix(nrow = length(countries),
                      ncol = length(countries))
colnames(phyl_jsd_mat) = countries
rownames(phyl_jsd_mat) = countries
for (c1 in countries) {
  for (c2 in countries){

    
    distance = get_phyl_jsd_dist(c1, c2)
    phyl_jsd_mat[c1, c2] = distance
  }
}

View(phyl_jsd_mat)

df_phyl_jsd_dist = data.frame(list(phyl_jsd_mat)[[1]])
df_phyl_jsd_dist$country_1 = colnames(df_phyl_jsd_dist)
df_phyl_jsd_dist = df_phyl_jsd_dist %>%
  mutate(country_1 = ifelse(country_1 == "NA.", "NA", country_1)) %>%
  pivot_longer(cols = !c(country_1), names_to = "country_2", values_to = "jsd_phyl_dist") %>%
  filter(country_1 != country_2) %>%
  arrange(country_1) %>%
  filter(country_1 < country_2) %>%
  mutate(country_1 = ifelse(country_1 == "NA.", "NA", country_1))
head(df_phyl_jsd_dist)

write_rds(phyl_jsd_mat, "data/processed/phyl_jsd_dist_mat.rds")

# join pairwise distances and export
df_distances = df_jsd %>%
  full_join(df_phyl_jsd_dist) %>%
  full_join(df_naive_div_dist) %>%
  full_join(df_phyl_div_dist)

write_csv(df_distances, "data/processed/pairwise_distances.csv")


