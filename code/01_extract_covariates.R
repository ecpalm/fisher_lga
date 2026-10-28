
###############################################################################
# 
# Script title: Covariate extraction
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
#' Shared covariate extraction -- used for iteration 0
#' (straight-line transects) and every LCP iteration.
#'
#' Lines are buffered by BUFFER_M (default 500m) before extraction.
#'
#' Dynamic covariates: canopy cover (cc_lf), burn severity (dnbr), binary
#' any-burn (burn_01), binary high-severity burn (burn_high_01), and years
#' since most recent burn (years_since_burn). Each has its own folder of
#' yearly rasters; all use the same lag from the pair's year.
#' ============================================================================

# Load required packages
sapply(c("foreach", "doParallel", "sf", "exactextractr", "terra", "purrr",
         "tibble", "stringr", "dplyr"),
       require, character.only = T)

#' ----------------------------------------------------------------------
#' Build straight-line LINESTRINGs (iteration 0 only)
#' ----------------------------------------------------------------------
lines_from_pairs <- function(gen_table, crs = 5070) {
  stopifnot(all(c("index", "lat1", "long1", "lat2", "long2") %in% names(gen_table)))
  
  purrr::map(seq_len(nrow(gen_table)), function(i) {
    sf::st_linestring(rbind(
      c(gen_table$long1[i], gen_table$lat1[i]),
      c(gen_table$long2[i], gen_table$lat2[i])
    ))
  }) %>%
    sf::st_sfc(crs = crs) %>%
    sf::st_sf(index = gen_table$index, geometry = .)
}

#' ----------------------------------------------------------------------
#' Static covariates
#'
#' `lines`: sf object with an `index` column matching input data. Either
#' straight lines or least cost paths.
#' ----------------------------------------------------------------------
extract_static_covariates <- function(lines, buffer_m = 500,
                                      covs_dir = "covs", already_buffered = FALSE) {
  stopifnot("index" %in% names(lines))
  
  var_names <- c("snowpack_early", "temp_min_early", "precip_early", "tri_45_riley", 
                 "roads", "water", "rivers")
  static_paths <- file.path(covs_dir, c("snowpack_early_150.tif", "temp_min_early_150.tif",
                                        "precip_early_150.tif", "tri_45_riley_150.tif", 
                                        "roads_150.tif", "water_150.tif", "rivers_150.tif"))
  n_layers <- length(static_paths)
  
  buffered <- if (already_buffered) lines else sf::st_buffer(lines, buffer_m)
  
  #' Parallelization
  out <- foreach::foreach(
    i = seq_len(nrow(buffered)), .combine = rbind,
    .packages = c("sf", "exactextractr", "terra")
  ) %dopar% {
    env_local <- terra::rast(static_paths)
    vals <- exactextractr::exact_extract(env_local, buffered[i, ])[[1]]
    vals <- vals[, seq_len(n_layers), drop = FALSE]
    c(index = buffered$index[i], colMeans(vals, na.rm = TRUE))
  }
  
  as.data.frame(out) %>%
    `names<-`(c("index", var_names)) %>%
    tibble::as_tibble()
}

#' ----------------------------------------------------------------------
#' Yearly raster stack for ONE dynamic variable. 
#' ----------------------------------------------------------------------
build_year_stack <- function(covs_dir, subdir, pattern = "\\.tif$") {
  paths <- list.files(file.path(covs_dir, subdir), full.names = TRUE, pattern = pattern)
  years <- as.integer(stringr::str_extract(basename(paths), "(19|20)\\d{2}"))
  stopifnot("couldn't parse a year from every file name in the folder" = !anyNA(years),
            "duplicate years in folder" = !anyDuplicated(years))
  o <- order(years)
  list(paths = paths[o], years = years[o])
}

#' One entry per dynamic variable: name = column name in the output table,
#' value = folder under covs_dir. 
build_dynamic_stacks <- function(covs_dir = "covs") {
  folders <- c(cc_lf = "cc_150", dnbr = "dnbr_150", burn_01 = "burn_01_150",
               burn_high_01 = "burn_high_01_150")
  lapply(folders, function(f) build_year_stack(covs_dir, f))
}

#' Per-variable handling of (a) target years BEFORE the raster series starts,
#' and (b) transects with no valid pixels (NaN mean).
#'  - cc_lf: clamp to the earliest year.
#'  - burn variables: before the series starts, "no burn"
dynamic_var_settings <- list(
  cc_lf            = list(before_start = "clamp", fill = NA_real_),
  dnbr             = list(before_start = "fill",  fill = 0),
  burn_01          = list(before_start = "fill",  fill = 0),
  burn_high_01     = list(before_start = "fill",  fill = 0)
)

#' ----------------------------------------------------------------------
#' Dynamic covariates (any set of variables), lagged from each pair's year.
#'
#' `lines`      : sf LINESTRING with an `index` column matching gen_table.
#' `gen_table`  : must have `index` and `year` (pairwise year used for lag).
#' `stacks`     : output of build_dynamic_stacks().
#' `vars`       : which dynamic variables to extract.
#' `fixed_year` : optional -- use this single year for every row.
#' ----------------------------------------------------------------------
extract_dynamic <- function(lines, gen_table, lag, stacks,
                            vars = names(stacks), settings = dynamic_var_settings,
                            fixed_year = NULL, buffer_m = 500, already_buffered = FALSE) {
  stopifnot("index" %in% names(lines),
            all(c("index", "year") %in% names(gen_table)),
            all(vars %in% names(stacks)), all(vars %in% names(settings)))
  
  buffered  <- if (already_buffered) lines else sf::st_buffer(lines, buffer_m)
  target_yr <- if (is.null(fixed_year)) {
    gen_table$year[match(buffered$index, gen_table$index)] - lag
  } else rep(fixed_year, nrow(buffered))
  
  ## layer index per row x variable (NA = before series start, handled by fill)
  layer_idx <- sapply(vars, function(v) {
    yrs <- stacks[[v]]$years
    yr  <- if (settings[[v]]$before_start == "clamp") pmax(target_yr, min(yrs)) else target_yr
    idx <- match(yr, yrs)
    bad <- is.na(idx) & !(yr < min(yrs))
    if (any(bad)) stop("years outside the ", v, " stack (after its last year or missing): ",
                       paste(unique(yr[bad]), collapse = ", "))
    idx
  })
  if (is.null(dim(layer_idx))) {
    layer_idx <- matrix(layer_idx, nrow = nrow(buffered), dimnames = list(NULL, vars))
  }
  
  paths_by_var <- lapply(stacks[vars], `[[`, "paths")
  fills        <- vapply(settings[vars], function(s) as.numeric(s$fill), numeric(1))
  
  out <- foreach::foreach(
    i = seq_len(nrow(buffered)), .combine = rbind,
    .packages = c("sf", "exactextractr", "terra")
  ) %dopar% {
    vals_i <- vapply(vars, function(v) {
      k <- layer_idx[i, v]
      if (is.na(k)) return(fills[[v]])
      r <- terra::rast(paths_by_var[[v]][k])
      x <- exactextractr::exact_extract(r, buffered[i, ])[[1]]
      m <- mean(x[["value"]], na.rm = TRUE)
      if (is.nan(m)) fills[[v]] else m
    }, numeric(1))
    c(index = buffered$index[i], vals_i)
  }
  
  as.data.frame(out) %>% tibble::as_tibble()
}

#' ----------------------------------------------------------------------
#' Static + dynamic in one call, joined by index column.
#' ----------------------------------------------------------------------
extract_all_covariates <- function(lines, gen_table, lag, stacks, vars = names(stacks),
                                   buffer_m = 500, covs_dir = "covs",
                                   already_buffered = FALSE) {
  static_tbl  <- extract_static_covariates(lines, buffer_m, covs_dir, already_buffered)
  dynamic_tbl <- extract_dynamic(lines, gen_table, lag, stacks, vars,
                                 buffer_m = buffer_m, already_buffered = already_buffered)
  
  static_tbl %>%
    dplyr::inner_join(dynamic_tbl, by = "index") %>%
    dplyr::inner_join(gen_table, by = "index")
}

#' ----------------------------------------------------------------------
#' Static-year counterpart, for the static-vs-dynamic comparison
#' (05_static_canopy_comparison.R). Every row uses the SAME fixed year.
#' ----------------------------------------------------------------------
extract_all_covariates_static <- function(lines, gen_table, year, stacks,
                                          vars = c("cc_lf"),
                                          buffer_m = 500, covs_dir = "covs",
                                          already_buffered = FALSE) {
  static_tbl  <- extract_static_covariates(lines, buffer_m, covs_dir, already_buffered)
  dynamic_tbl <- extract_dynamic(lines, gen_table, lag = 0, stacks, vars,
                                 fixed_year = year, buffer_m = buffer_m,
                                 already_buffered = already_buffered)
  
  static_tbl %>%
    dplyr::inner_join(dynamic_tbl, by = "index") %>%
    dplyr::inner_join(gen_table, by = "index")
}
