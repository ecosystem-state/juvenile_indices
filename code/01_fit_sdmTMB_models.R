library(sdmTMB)
library(ggplot2)

species = read.csv("data/species_list_revised.csv")
names(species) = tolower(names(species))
species = dplyr::rename(species,
                        common_name = common.name,
                        scientific_name = scientific.name,
                        juv_threshold = max.length.cm)

for(sp_num in 1:nrow(species)){
  comm_name = species$common_name[sp_num]
  d = readRDS(file=paste0("output/", sub(" ", "_", comm_name), "_expanded_2023.rds"))

  # filter out years with no data
  drop <- dplyr::group_by(d, year) %>%
    dplyr::summarize(s = sum(juv_cpue_kg_km2)) %>%
    dplyr::filter(s == 0)
  if(nrow(drop) > 0) {
    d = dplyr::filter(d, year %in% drop$year == FALSE)
  }

  d = add_utm_columns(d, ll_names = c("longitude_dd","latitude_dd"))
  mesh = make_mesh(d, xy_cols = c("X","Y"),
                   cutoff = 50)
  d$yearf <- as.factor(d$year)
  d$log_depth <- scale(log(d$depth_hi_prec_m))

  fit = sdmTMB(juv_cpue_kg_km2 ~ -1 + yearf + log_depth + I(log_depth^2),
               mesh = mesh,
               family = tweedie(),
               spatial = "on",
               spatiotemporal = "iid",
               time="year",
               data = d)

  grid = readRDS("data/wc_grid.rds")
  grid$log_depth <- log(grid$depth)
  grid$XY <- paste(grid$X, grid$Y)
  new_grid = expand.grid(XY = unique(grid$XY),
                         year = unique(fit$data$year))
  new_grid = dplyr::left_join(new_grid, grid)
  new_grid$yearf <- as.factor(new_grid$year)
  pred = predict(fit, new_grid, return_tmb_object = TRUE)
  index = get_index(pred, bias_correct = TRUE)
  index$common_name <- comm_name
  saveRDS(index, file=paste0("output/", sub(" ", "_", comm_name), "_index_2023.rds"))

  if(sp_num == 1) {
    d <- index
  } else {
    d <- rbind(d, index)
  }
}

saveRDS(d, "data/all_juvenile_indices.rds")

# Plot de-meaned trends on top of eachother
dplyr::filter(d, !is.na(se)) %>%
  dplyr::group_by(common_name) %>%
  dplyr::mutate(x = log_est - mean(log_est)) %>%
  ggplot(aes(year, x, col = common_name)) +
  #geom_ribbon(aes(ymin=log_est - 2*se, ymax = log_est + 2*se)) +
  geom_line()
#facet_wrap(~common_name, scale="free_y")
