
###############################################################################
# 
# Script title: Randomized shortest paths to predict genetic connectivity
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
#' Randomized shortest paths (RSP) for the early-period subgroups for 2 
#' representative years (2012 = pre-disturbance, 2023 = current) 
#' ============================================================================

# Load required packages
sapply(c("raster", "sp", "gdistance", "terra", "sf", "units", "smoothr",
         "stringr", "foreach", "doParallel", "tibble", "purrr", "dplyr"),
       require, character.only = T)

source("code/lcp_tools.R")
source("code/cv_diagnostics.R") # compute_fold_metrics(), for find_best_model()

#' ---- User inputs ---- 
response_col    <- "GD"
date_type       <- "max"
model_dir       <- "models/"
out_dir         <- "outputs/rsp/"
sdm_dir         <- "covs/hart_habitat/raw"
waterbody_shp   <- "shp/waterbody_over_100_ha.shp"
regions_shp     <- "shp/regions_final.gpkg"

REGION_COL      <- "region"

subgroups       <- c("male_north_early", "male_south_early",
                      "female_north_early", "female_south_early")

subgroups_to_run <- subgroups
#' ----------------------------------------------------------------

years           <- c(2012, 2023)     # 2012 = pre-disturbance, 2023 = current
theta_values    <- c(0.0001, 0.001, 0.01) # 3 exploration parameter values

target_n_pairs  <- 5000  # target number of pairs of start/end locations
                         # sampled from the spatial union of both years' habitat
                         # patches; each year's own final count is whatever of 
                         # those pairs survives filtering into that year's patches 

PRINT_DIAGNOSTICS_ONLY <- TRUE   # set TRUE to avoid running the RSP step if you
                                 # only want final pair numbers -- lets you check those
                                   # numbers without waiting on or redoing the RSP computation.

#' Distance-scaled crop padding for each pair of start/end points to speed up RSP 
crop_buffer_base_m <- 5000  # fixed minimum padding, meters
crop_buffer_frac    <- 0.2  # additional padding as a fraction of pair_distance

mm_cache_path      <- "results/rsp_normalization_range.rds"
FORCE_RECOMPUTE_MM <- FALSE   # set TRUE to rebuild the shared normalization range 
min_area_km2    <- 1
smooth_amt      <- 1.5
FORCE_RECOMPUTE <- FALSE
noCLS           <- 100
seed            <- 1234
#' ---------------------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fold_results      <- readRDS("results/fisher_cv_folds.rds")
selection_results <- readRDS("results/fisher_variable_selection.rds")
pairs <- readRDS("data/pairs_fisher_geomasked.rds") %>% 
  dplyr::mutate(index = dplyr::row_number())

regions_sf <- sf::st_read(regions_shp, quiet = TRUE)

#' One polygon (north or south) per subgroup's tag$region -- used as the
#' domain for resistance prediction/patches/point sampling/output extent
get_region_poly <- function(region_tag) {
  match_rows <- regions_sf %>% dplyr::filter(tolower(.data[[REGION_COL]]) == tolower(region_tag))
  if (nrow(match_rows) != 1) {
    stop("Expected exactly 1 polygon in ", regions_shp, " matching region = '", region_tag,
         "' (column '", REGION_COL, "'), found ", nrow(match_rows), ". Available values: ",
         paste(unique(regions_sf[[REGION_COL]]), collapse = ", "))
  }
  match_rows
}

min_max_dist <- tibble::tibble(
  sex      = c("male", "female"),
  min_dist = c(0, 0),
  max_dist = c(50000, 30000)
)

lower_df <- function(x) {
  mat <- as.matrix(x)
  ind <- which(lower.tri(mat, diag = FALSE), arr.ind = TRUE)
  nn <- dimnames(mat)
  tibble::tibble(id_1 = nn[[1]][ind[, 1]], id_2 = nn[[2]][ind[, 2]], d = mat[ind])
}

out_path <- function(label, year, theta) {
  paste0(out_dir, "rsp_", label, "_", year, "_theta_", theta, ".tif")
}

#' ----------------------------------------------------------------------
#' Best (iteration, lag) model per subgroup, by pooled RMSE 
#' ----------------------------------------------------------------------
find_best_model <- function(tag, period) {
  files <- list.files(
    model_dir,
    pattern = paste0("^lcp_[1-9][0-9]*_", tag$sex, "_", tag$region, "_lag_[0-9]+_",
                      date_type, "_gbm_", period, "\\.rds$"),
    full.names = TRUE
  )
  stopifnot("no dynamic LCP models found for this subgroup -- run 03/04 first" = length(files) > 0)

  results <- files %>%
    purrr::map_dfr(function(f) {
      fit  <- readRDS(f)
      rmse <- compute_fold_metrics(dplyr::semi_join(fit$pred, fit$bestTune))$global$rmse
      m    <- stringr::str_match(basename(f), "^lcp_(\\d+)_.*_lag_(\\d+)_")
      tibble::tibble(iteration = as.integer(m[, 2]), lag = as.integer(m[, 3]), rmse = rmse, file = f)
    })

  best <- results %>%
    dplyr::arrange(rmse, iteration, lag) %>%
    dplyr::slice(1)

  stopifnot(
    "find_best_model() produced 0 rows -- every candidate file failed to score an RMSE" =
      nrow(best) == 1,
    "find_best_model()'s chosen file path is blank/NA" =
      is.character(best$file) && length(best$file) == 1 && !is.na(best$file) && nzchar(best$file),
    "find_best_model()'s chosen file doesn't exist on disk" =
      file.exists(best$file)
  )

  best
}

#' ----------------------------------------------------------------------
#' Habitat patches -- raw SDM predictions loaded once per year, then
#' resampled/cropped/thresholded/smoothed once per subgroup-year.
#' ----------------------------------------------------------------------
sdm_cache <- years %>%
  purrr::set_names(as.character) %>%
  purrr::map(~ terra::rast(file.path(sdm_dir, paste0("Processed_Continuous_SDM_", .x, ".tif"))))

build_patches <- function(sdm_raw, template) {
  crop_ext <- template %>%
    terra::classify(cbind(-Inf, Inf, 1)) %>%
    terra::as.polygons() %>%
    terra::buffer(-255)

  sdm_raw %>%
    terra::resample(template) %>%
    terra::crop(crop_ext, mask = TRUE) %>%
    terra::classify(cbind(-Inf, 4999.5, NA)) %>%
    terra::classify(cbind(4999.5, Inf, 1)) %>%
    terra::as.polygons() %>%
    sf::st_as_sf() %>%
    sf::st_cast("POLYGON") %>%
    dplyr::mutate(area_km2 = as.numeric(units::set_units(sf::st_area(.), "km2"))) %>%
    dplyr::filter(area_km2 >= min_area_km2) %>%
    smoothr::smooth(method = "ksmooth", smooth = smooth_amt) %>%
    dplyr::mutate(poly_id = dplyr::row_number())
}

#' ----------------------------------------------------------------------
#' Point sampling with automated target-connection-count matching
#' ----------------------------------------------------------------------
sample_pairs <- function(patches, min_dist, max_dist, target_n, seed,
                          initial_n = target_n * 3, growth = 1.5, max_tries = 6) {
  n <- initial_n
  valid_pairs <- NULL
  pts_raw <- NULL
  for (try in seq_len(max_tries)) {
    set.seed(seed + try)   # vary per try 
    pts_raw <- patches %>% sf::st_sample(n) %>% sf::st_sf() %>% dplyr::mutate(id = dplyr::row_number())

    poly_lookup <- pts_raw %>%
      sf::st_join(patches) %>%
      sf::st_drop_geometry() %>%
      tibble::as_tibble() %>%
      dplyr::select(id, poly_id)

    valid_pairs <- pts_raw %>%
      sf::st_distance() %>%
      units::drop_units() %>%
      lower_df() %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
      dplyr::filter(d > min_dist, d <= max_dist) %>%
      dplyr::inner_join(poly_lookup %>% dplyr::rename(id_1 = id, poly_id_1 = poly_id), by = "id_1") %>%
      dplyr::inner_join(poly_lookup %>% dplyr::rename(id_2 = id, poly_id_2 = poly_id), by = "id_2") %>%
      dplyr::filter(poly_id_1 != poly_id_2)

    if (nrow(valid_pairs) >= target_n) {
      set.seed(seed)   # reseed so the final subsample is reproducible 
      return(list(pts = pts_raw, pairs = valid_pairs %>% dplyr::slice_sample(n = target_n)))
    }
    message("    only ", nrow(valid_pairs), " valid pairs from ", n, " points -- retrying with more points")
    n <- ceiling(n * growth)
  }
  stop("Could not reach target_n_pairs = ", target_n, " after ", max_tries, " tries (last had ",
       nrow(valid_pairs), " valid pairs from ", n, " points) -- raise initial_n/growth/max_tries, or lower target_n.")
}

#' ----------------------------------------------------------------------
#' Filter one shared set of start/end location pairs down to those that are 
#' valid in a specific year 
#' ----------------------------------------------------------------------
filter_pairs_to_year <- function(pts_sf, pairs_df, patches_this_year) {
  poly_lookup <- pts_sf %>%
    sf::st_join(patches_this_year) %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble() %>%
    dplyr::select(id, poly_id)

  pairs_df %>%
    dplyr::select(id_1, id_2, d) %>%
    dplyr::inner_join(poly_lookup %>% dplyr::rename(id_1 = id, poly_id_1 = poly_id), by = "id_1") %>%
    dplyr::inner_join(poly_lookup %>% dplyr::rename(id_2 = id, poly_id_2 = poly_id), by = "id_2") %>%
    dplyr::filter(!is.na(poly_id_1), !is.na(poly_id_2), poly_id_1 != poly_id_2)
}

#' ----------------------------------------------------------------------
#' RSP for one subgroup-year's fixed set of point pairs, at one theta 
#' ----------------------------------------------------------------------
#' n_total: the FIXED denominator to divide the summed passage rasters by. A
#' fixed shared denominator means a year that lost more connectivity routes
#' correctly shows lower average passage almost everywhere it lost them.
run_rsp <- function(resistance_raster, pts_sf, pairs_ids, theta,
                     crop_buffer_base_m, crop_buffer_frac, noCLS, n_total) {

  st_pts     <- pts_sf %>% sf::as_Spatial()
  pts_combn  <- pairs_ids %>% dplyr::select(id_1, id_2) %>% as.matrix()
  r_full     <- raster::raster(resistance_raster)   # shipped to the cluster ONCE below

  cl <- parallel::makeCluster(noCLS)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl, varlist = "r_full", envir = environment())

  results <- tryCatch({
    foreach::foreach(
      i = seq_len(nrow(pts_combn)), .verbose = FALSE,
      .packages = c("gdistance", "raster", "sp", "terra")
    ) %dopar% {

      p1 <- st_pts@coords[pts_combn[i, 1], 1:2]
      p2 <- st_pts@coords[pts_combn[i, 2], 1:2]
      locations <- sp::SpatialPoints(rbind(p1, p2), proj4string = sp::CRS(sp::proj4string(st_pts)))

      pair_dist <- sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2)
      buffer_m  <- crop_buffer_base_m + crop_buffer_frac * pair_dist

      bbox_pair <- raster::extent(
        min(p1[1], p2[1]) - buffer_m, max(p1[1], p2[1]) + buffer_m,
        min(p1[2], p2[2]) - buffer_m, max(p1[2], p2[2]) + buffer_m
      )
      r_crop <- raster::crop(r_full, bbox_pair)

      tr <- gdistance::transition(r_crop, function(x) 1 / mean(x), directions = 8) %>%
        gdistance::geoCorrection(type = "c", multpl = FALSE)
      min_pos <- min(tr@transitionMatrix@x[tr@transitionMatrix@x > 0])
      tr@transitionMatrix@x[tr@transitionMatrix@x < 0] <- min_pos

      pass <- gdistance::passage(tr, origin = locations[1], goal = locations[2], theta = theta)

      #' place the cropped result back onto the FULL grid (0 outside the
      #' crop) so every pair's output sums cell-for-cell into one composite
      raster::extend(pass, r_full, value = 0) %>% raster::values()
    }
  }, finally = parallel::stopCluster(cl))

  passages_sum <- Reduce(`+`, results) / n_total

  #' re-mask against resistance_raster's own NA pattern 
  terra::rast(resistance_raster) %>%
    terra::setValues(passages_sum) %>%
    terra::mask(resistance_raster)
}

#' ----------------------------------------------------------------------
#' Pass 1: predict RAW (not normalized) resistance for all subgroups/years
#' ----------------------------------------------------------------------
subgroup_setup <- subgroups_to_run %>%
  purrr::set_names() %>%
  purrr::map(function(label) {
    folds <- fold_results[[label]]
    sel   <- selection_results[[label]]
    if (is.null(folds) || is.null(sel)) { message("Skipping ", label, " (no folds/selection)"); return(NULL) }

    sub_pairs      <- pairs %>% dplyr::filter(index %in% folds$orig_index)
    tag            <- subgroup_tag(sub_pairs)
    period         <- sub("^.*_(early|late)$", "\\1", label)
    sex            <- sub("^(male|female)_.*$", "\\1", label)

    message("== ", label, " ==")
    best     <- find_best_model(tag, period)
    best_fit <- readRDS(best$file)
    message("  best model: iteration ", best$iteration, ", lag ", best$lag, " (RMSE = ", round(best$rmse, 3), ")")

    predictor_cols <- setdiff(colnames(best_fit$trainingData), ".outcome")

    list(tag = tag, sex = sex, fit = best_fit, predictor_cols = predictor_cols,
         domain = get_region_poly(tag$region), euc_geog_value = stats::median(sub_pairs$euc_geog),
         dist_bounds = min_max_dist %>% dplyr::filter(sex == .env$sex))
  }) %>%
  purrr::compact()

raw_resistance <- subgroup_setup %>%
  purrr::imap(function(s, label) {
    years %>%
      purrr::set_names(as.character) %>%
      purrr::map(function(yr) {
        message(label, " | ", yr, ": predicting raw resistance surface")
        predict_resistance_surface(s$fit, s$predictor_cols, yr, s$euc_geog_value, s$domain)
      })
  })

#' SHARED NORMALIZATION 
if (!FORCE_RECOMPUTE_MM && file.exists(mm_cache_path)) {

  mm <- readRDS(mm_cache_path)
  message("Using cached normalization range from ", mm_cache_path, ": [",
          round(mm$min, 4), ", ", round(mm$max, 4),
          "] -- delete this file or set FORCE_RECOMPUTE_MM <- TRUE to rebuild it.")

} else {

  mm <- raw_resistance %>%
    purrr::flatten() %>%
    purrr::map_dfr(~ tibble::tibble(min = terra::global(.x, "min", na.rm = TRUE)$min,
                                     max = terra::global(.x, "max", na.rm = TRUE)$max)) %>%
    dplyr::summarise(min = min(min), max = max(max))

    dir.create(dirname(mm_cache_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(mm, mm_cache_path)
    message("Computed and cached normalization range to ", mm_cache_path, ": [",
            round(mm$min, 4), ", ", round(mm$max, 4), "]")
}

#' ----------------------------------------------------------------------
#' Pass 2: normalize (shared min/max) -> patches -> sample pairs -> RSP per theta. 
#' Patches are first built for both years; point pairs are sampled once per subgroup 
#' from the union of both years' habitat patches and reused across each theta value for
#' both years.
#' ----------------------------------------------------------------------
all_failures          <- character(0)
patches_by_region_year <- list()   

for (label in names(subgroup_setup)) {

  s   <- subgroup_setup[[label]]
  tag <- s$tag

  r_norm_by_year <- years %>%
    purrr::set_names(as.character) %>%
    purrr::map(function(yr) {
      raw_r <- raw_resistance[[label]][[as.character(yr)]]
     
      terra::app(raw_r, function(x) (x - mm$min) / (mm$max - mm$min) + .01) %>%
        terra::mask(waterbody, inverse = TRUE, updatevalue = 50) %>%
        terra::mask(raw_r)
    })

  message(label, ": building habitat patches (", paste(years, collapse = ", "), ")")
  patches_by_year <- years %>%
    purrr::set_names(as.character) %>%
    purrr::map(~ build_patches(sdm_cache[[as.character(.x)]], raw_resistance[[label]][[as.character(.x)]]))

  #' Per-subgroup patch count/area, printed for both years. 
  years %>% purrr::walk(function(yr) {
    p <- patches_by_year[[as.character(yr)]]
    message("    ", label, " | ", yr, " patches: ", nrow(p), " (", round(sum(p$area_km2), 1), " km2 total)")
  })
  patches_by_region_year[[tag$region]] <- patches_by_year

  #' Union-based sampling 
  patches_union <- dplyr::bind_rows(
    patches_by_year[[as.character(years[1])]] %>% sf::st_geometry() %>% sf::st_sf(),
    patches_by_year[[as.character(years[2])]] %>% sf::st_geometry() %>% sf::st_sf()
  ) %>%
    sf::st_union() %>%
    sf::st_cast("POLYGON") %>%
    sf::st_sf() %>%
    dplyr::mutate(poly_id = dplyr::row_number())

  message(label, " | union of ", paste(years, collapse = "/"), ": sampling point pairs (target = ", target_n_pairs, ")")
  union_sampled <- sample_pairs(patches_union, s$dist_bounds$min_dist, s$dist_bounds$max_dist,
                                 target_n_pairs, seed = seed)
  message(label, " | union: ", nrow(union_sampled$pairs), " point pairs sampled")

  sampled_by_year <- years %>%
    purrr::set_names(as.character) %>%
    purrr::map(function(yr) {
      list(pts   = union_sampled$pts,
           pairs = filter_pairs_to_year(union_sampled$pts, union_sampled$pairs,
                                         patches_by_year[[as.character(yr)]]))
    })

  #' Surviving pair count per year 
  years %>% purrr::walk(function(yr) {
    message("    ", label, " | ", yr, ": ", nrow(sampled_by_year[[as.character(yr)]]$pairs),
            " of ", nrow(union_sampled$pairs), " union-sampled pairs valid",
            " -- raise target_n_pairs if either year comes out too small to map meaningfully.")
  })

  if (PRINT_DIAGNOSTICS_ONLY) {
    message(label, ": PRINT_DIAGNOSTICS_ONLY is TRUE -- skipping RSP computation for this subgroup.")
    next
  }

  for (year in years) {

    output_files <- theta_values %>% purrr::map_chr(~ out_path(label, year, .x))
    if (!FORCE_RECOMPUTE && all(file.exists(output_files))) {
      message(label, " | ", year, ": all theta outputs already exist, skipping (set FORCE_RECOMPUTE <- TRUE to redo).")
      next
    }

    tryCatch({

      sampled <- sampled_by_year[[as.character(year)]]
      if (nrow(sampled$pairs) == 0) {
        message(label, " | ", year, ": 0 valid pairs -- skipping (nothing to map). Raise target_n_pairs and rerun.")
        next
      }

      for (theta in theta_values) {
        this_out <- out_path(label, year, theta)
        if (!FORCE_RECOMPUTE && file.exists(this_out)) {
          message("    theta ", theta, ": already computed, skipping.")
          next
        }
        message("    theta ", theta, " (", year, "): running RSP over ", nrow(sampled$pairs), " pairs")
        passages_overlay <- run_rsp(r_norm_by_year[[as.character(year)]], sampled$pts, sampled$pairs, theta,
                                     crop_buffer_base_m, crop_buffer_frac, noCLS,
                                     n_total = nrow(union_sampled$pairs))
        terra::writeRaster(passages_overlay, this_out, overwrite = TRUE)
      }

    }, error = function(e) {
      message("FAILED: ", label, " | ", year, " -- ", conditionMessage(e))
      all_failures <<- c(all_failures, paste0(label, " | ", year))
    })

    gc()
  }
}

#' Region-combined patch summary 
message("\n=== Habitat patch summary, combined across regions ===")
years %>% purrr::walk(function(yr) {
  combined <- patches_by_region_year %>%
    purrr::map_dfr(~ .x[[as.character(yr)]] %>% sf::st_drop_geometry())
  message("  ", yr, ": ", nrow(combined), " patches total, ", round(sum(combined$area_km2), 1), " km2 total")
})

if (length(all_failures) > 0) {
  message("\n=== ", length(all_failures), " (subgroup, year) combination(s) failed this run ===")
  message(paste("  -", all_failures, collapse = "\n"))
} else {
  message("\nAll subgroup/year/theta combinations completed without error.")
}
