
###############################################################################
# 
# Script title: Variable selection
#
# Article title: Dynamic landscape genetics predicts disturbance-driven
#                declines in genetic connectivity in an endangered forest 
#                carnivore
#
# Article authors: Eric Palm, Erin Landguth, Katherine Zeller, Rebecca Green,
#                  Kristine Pilgrim, Sarah Stock, Craig Thompson, 
#                  Michael Schwartz, Jody Tucker
#
# Script author: Eric Palm
#
###############################################################################

#' ============================================================================
#' LCP tools 
#' ============================================================================

# Load required packages
sapply(c("raster", "sp", "gdistance", "terra", "sf", "foreach", "purrr",
         "tibble", "dplyr"),
       require, character.only = T)

source("code/01_extract_covariates.R")
source("code/fit_tools.R")

#' ---- User inputs  ----
period_years <- tibble::tibble(
  period    = c("early", "late"),
  year_min  = c(2006, 2016),   # 
  year_max  = c(2015, 2023)    # 
)
covs_dir      <- "covs"
waterbody_shp <- "shp/waterbody_over_100_ha.shp"
#' ---------------------------------------------------------------

stacks <- build_dynamic_stacks()

#' Static single-file covariate rasters, loaded once and reused across
#' every year/subgroup/lag/iteration call that needs them.
static_rasters <- list(
  snowpack_early = terra::rast(file.path(covs_dir, "snowpack_early_150.tif")),
  temp_min_early = terra::rast(file.path(covs_dir, "temp_min_early_150.tif")),
  precip_early   = terra::rast(file.path(covs_dir, "precip_early_150.tif")),
  rivers   = terra::rast(file.path(covs_dir, "rivers_150.tif")),
  water   = terra::rast(file.path(covs_dir, "water_150.tif")),
  roads   = terra::rast(file.path(covs_dir, "roads_150.tif"))
)

waterbody <- terra::vect(waterbody_shp)

#' ----------------------------------------------------------------------
#' Tool: consistent sex/region tagging from a subgroup's pairs table.
#' ----------------------------------------------------------------------
subgroup_tag <- function(sub_pairs) {
  list(sex = tolower(sub_pairs$sex[1]), region = sub_pairs$region[1])
}

#' ----------------------------------------------------------------------
#' Build one year's prediction raster stack, restricted to whatever
#' predictor_cols are included. euc_geog doesn't
#' contribute to the resistance surface's spatial pattern.
#' ----------------------------------------------------------------------
build_prediction_stack <- function(year, predictor_cols, euc_geog_value, domain) {
  layers <- list()
  for (v in intersect(predictor_cols, names(stacks))) {
    yrs <- stacks[[v]]$years
    layer_idx <- match(if (dynamic_var_settings[[v]]$before_start == "clamp") max(year, min(yrs)) else year, yrs)
    if (is.na(layer_idx)) {
      if (!(year < min(yrs))) stop("requested year falls outside the ", v, " stack's years")
      layers[[v]] <- static_rasters[[1]] %>% terra::setValues(dynamic_var_settings[[v]]$fill)
    } else {
      layers[[v]] <- terra::rast(stacks[[v]]$paths[layer_idx])
    }
  }
  for (v in intersect(predictor_cols, names(static_rasters))) layers[[v]] <- static_rasters[[v]]
  if ("euc_geog" %in% predictor_cols) {
    layers$euc_geog <- static_rasters[[1]] %>% terra::setValues(euc_geog_value)
  }
  terra::rast(layers) %>% `names<-`(names(layers)) %>% terra::crop(terra::vect(domain), mask = TRUE)
}

predict_resistance_surface <- function(fit, predictor_cols, year, euc_geog_value, domain) {
  env <- build_prediction_stack(year, predictor_cols, euc_geog_value, domain)
  terra::predict(env[[predictor_cols]], fit) %>%
    terra::crop(terra::vect(domain), mask = TRUE) %>%
    `names<-`(as.character(year))
}

#' ----------------------------------------------------------------------
#' LCP line computation -- gdistance 
#' Per-year transition matrices; per-transect shortest paths via
#' foreach %dopar%, buffered 500m, water-erased. Uses whatever parallel
#' backend the CALLER already registered -- see the CLUSTER MANAGEMENT
#' note in 01_extract_covariates.R.
#' ----------------------------------------------------------------------
compute_lcp_lines <- function(resistance_stack_by_year, sub_pairs, lag, domain) {

  r_norm <- resistance_stack_by_year %>%
    purrr::map(function(r) {
      max_val <- terra::global(r, "max", na.rm = TRUE)$max
      min_val <- terra::global(r, "min", na.rm = TRUE)$min
      terra::app(r, function(x) (x - min_val) / (max_val - min_val) + .01) %>%
        terra::mask(r) %>%
        terra::mask(waterbody, inverse = TRUE, updatevalue = 50)
    })

  years <- names(resistance_stack_by_year)

  tr_list <- r_norm %>%
    purrr::map(function(r) {
      tempie <- gdistance::transition(raster::raster(r), function(x) 1 / mean(x), directions = 8) %>%
        gdistance::geoCorrection(type = "c", multpl = FALSE)
      min_pos <- min(tempie@transitionMatrix@x[tempie@transitionMatrix@x > 0])
      tempie@transitionMatrix@x[tempie@transitionMatrix@x < 0] <- min_pos
      tempie
    })

  # Mask by large waterbodies to remove them from extractions
  masker_sf <- waterbody %>% terra::intersect(terra::vect(domain)) %>% sf::st_as_sf()

  masker_geom <- sf::st_make_valid(sf::st_geometry(masker_sf))
  masker_geom <- if (is.na(sf::st_crs(masker_geom))) sf::st_set_crs(masker_geom, 5070) else sf::st_transform(masker_geom, 5070)

  foreach::foreach(
    i = seq_len(nrow(sub_pairs)), .combine = rbind, .verbose = FALSE,
    .packages = c("sp", "sf", "raster", "dplyr", "gdistance")
  ) %dopar% {
    idx   <- sub_pairs$index[i]
    long1 <- sub_pairs$long1[i]; lat1 <- sub_pairs$lat1[i]
    long2 <- sub_pairs$long2[i]; lat2 <- sub_pairs$lat2[i]

    tr_year <- sub_pairs$year[i] - lag - (as.integer(min(years)) - 1)
    tr_year <- max(tr_year, 1)

    gdistance::shortestPath(tr_list[[tr_year]],
                             sp::SpatialPoints(cbind(long1, lat1)),
                             sp::SpatialPoints(cbind(long2, lat2)),
                             output = "SpatialLines") %>%
      sf::st_as_sf() %>%
      sf::st_set_crs(5070) %>%
      sf::st_buffer(500) %>%
      (function(b) {
        near <- masker_geom[sf::st_intersects(b, masker_geom)[[1]]]
        if (length(near) == 0) b else sf::st_difference(b, sf::st_union(near))
      })() %>%
      dplyr::mutate(index = idx)
  }
}

#' ----------------------------------------------------------------------
#' Static-surface variant of compute_lcp_lines(), for the static-vs-
#' dynamic comparison (05_static_canopy_comparison.R). Takes
#' ONE resistance raster (not a per-year stack) and builds ONE transition
#' matrix shared by every pair
#' ----------------------------------------------------------------------
compute_lcp_lines_static <- function(resistance, sub_pairs, domain) {

  max_val <- terra::global(resistance, "max", na.rm = TRUE)$max
  min_val <- terra::global(resistance, "min", na.rm = TRUE)$min
  r_norm <- terra::app(resistance, function(x) (x - min_val) / (max_val - min_val) + .01) %>%
    terra::mask(resistance) %>%
    terra::mask(waterbody, inverse = TRUE, updatevalue = 50)

  tr <- gdistance::transition(raster::raster(r_norm), function(x) 1 / mean(x), directions = 8) %>%
    gdistance::geoCorrection(type = "c", multpl = FALSE)
  min_pos <- min(tr@transitionMatrix@x[tr@transitionMatrix@x > 0])
  tr@transitionMatrix@x[tr@transitionMatrix@x < 0] <- min_pos

  masker_sf <- waterbody %>% terra::intersect(terra::vect(domain)) %>% sf::st_as_sf()

  masker_geom <- sf::st_make_valid(sf::st_geometry(masker_sf))
  masker_geom <- if (is.na(sf::st_crs(masker_geom))) sf::st_set_crs(masker_geom, 5070) else sf::st_transform(masker_geom, 5070)

  foreach::foreach(
    i = seq_len(nrow(sub_pairs)), .combine = rbind, .verbose = FALSE,
    .packages = c("sp", "sf", "raster", "dplyr", "gdistance")
  ) %dopar% {
    idx   <- sub_pairs$index[i]
    long1 <- sub_pairs$long1[i]; lat1 <- sub_pairs$lat1[i]
    long2 <- sub_pairs$long2[i]; lat2 <- sub_pairs$lat2[i]

    gdistance::shortestPath(tr,
                             sp::SpatialPoints(cbind(long1, lat1)),
                             sp::SpatialPoints(cbind(long2, lat2)),
                             output = "SpatialLines") %>%
      sf::st_as_sf() %>%
      sf::st_set_crs(5070) %>%
      sf::st_buffer(500) %>%
      (function(b) {
        near <- masker_geom[sf::st_intersects(b, masker_geom)[[1]]]
        if (length(near) == 0) b else sf::st_difference(b, sf::st_union(near))
      })() %>%
      dplyr::mutate(index = idx)
  }
}
