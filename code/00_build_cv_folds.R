
###############################################################################
# 
# Script title: Create spatial cross validation folds by sex-region-period
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
#' Creating transect-based kNNDM spatial CV folds for gradient-boosting machine
#' with pairwise genetic distance as response variable
#' Separate CV folds per sex x region x period subgroup.
#'
#'
#' Necessary inputs:
#'   pairs: data frame, one row per candidate pairwise transect
#'          Columns: index, sex, id_1, id_2, long1, lat1, long2, lat2, region,
#'          euc_geog, date_1, date_2, year, 
#'
#' predpoints are sampled from a per-subgroup minimum convex polygon (MCP) around that
#' subgroup's actual sample locations, buffered by mcp_buffer_m (default 10 km)
#'
#' ============================================================================

# Load packages
sapply(c("dplyr", "tidyr", "purrr", "tibble", "stringr", "caret", "sf",
         "units", "twosamples", "lubridate", "CAST"),
       require, 
       character.only = T)

set.seed(1234)

pairs <- readRDS("data/pairs_fisher_geomasked.rds") %>%
  dplyr::mutate(index = dplyr::row_number())

#' ----------------------------------------------------------------------
#' Calculate medoid spatial coordinate for fisher with multiple genetic samples
#' ----------------------------------------------------------------------
get_medoid <- function(coords) {
  if (nrow(coords) == 1) return(coords)
  d <- as.matrix(stats::dist(coords))
  coords[which.min(rowSums(d)), , drop = FALSE]
}

#' ----------------------------------------------------------------------
#' Minimum convex polygon around a set of sf points 
#' ----------------------------------------------------------------------
mcp <- function(x, percentile = 100) {
  centroid <- sf::st_centroid(sf::st_union(x))
  dist <- as.numeric(sf::st_distance(x, centroid))
  within_range <- dist <= stats::quantile(dist, percentile / 100)
  sf::st_convex_hull(sf::st_union(x[within_range, ]))
}

#' ----------------------------------------------------------------------
#' Transect-based kNNDM: uses point-to-nearest-TRANSECT
#' distance rather than point-to-point distance, so sampling unit is a 
#' pair of locations rather than a single location.
#' ----------------------------------------------------------------------
knndm_transect <- function(tpoints, pt2transect, pred2transect, transect_ends,
                            k = 10, maxp = 0.5, linkf = "ward.D2", nk_len = 100) {

  tcoords <- sf::st_coordinates(tpoints)
  n_animals <- nrow(tcoords)
  gij <- apply(pred2transect, 1, min)

  cv_distances_transect <- function(clust_k) {
    vapply(seq_len(n_animals), function(i) {
      f <- clust_k[i]
      train_cols <- !(clust_k[transect_ends[, 1]] == f | clust_k[transect_ends[, 2]] == f)
      if (!any(train_cols)) return(NA_real_)
      min(pt2transect[i, train_cols])
    }, numeric(1))
  }

  hc <- stats::hclust(stats::dist(tcoords), method = linkf)

  clustgrid <- data.frame(
    nk = unique(as.integer(round(exp(seq(log(k), log(n_animals - 2), length.out = nk_len)))))
  )
  clustgrid$w_stat <- NA_real_
  clustgroups <- list()
  pcacoords <- stats::prcomp(tcoords, center = TRUE, scale. = FALSE, rank. = 1)

  for (nk in clustgrid$nk) {

    clust_nk <- stats::cutree(hc, k = nk)
    tabclust <- as.data.frame(table(clust_nk))
    tabclust$clust_k <- NA_integer_

    centr <- vapply(tabclust$clust_nk, function(x) {
      m <- matrix(colMeans(tcoords[clust_nk %in% x, , drop = FALSE]), nrow = 1)
      colnames(m) <- colnames(tcoords)
      stats::predict(pcacoords, m)
    }, numeric(1))
    tabclust$centrpca <- centr
    tabclust <- tabclust[order(tabclust$centrpca), ]

    clust_i <- 1
    for (i in seq_len(nrow(tabclust))) {
      if (tabclust$Freq[i] >= n_animals / k) {
        tabclust$clust_k[i] <- clust_i
        clust_i <- clust_i + 1
      }
    }
    remaining <- setdiff(seq_len(k), unique(tabclust$clust_k))
    if (length(remaining) > 0) {
      na_idx <- is.na(tabclust$clust_k)
      tabclust$clust_k[na_idx] <- rep(remaining, ceiling(nk / length(remaining)))[seq_len(sum(na_idx))]
    }

    id_map <- merge(data.frame(pos = seq_along(clust_nk), clust_nk = clust_nk), tabclust, by = "clust_nk")
    id_map <- id_map[order(id_map$pos), ]
    clust_k <- id_map$clust_k

    if (!any(table(clust_k) / length(clust_k) > maxp)) {
      gjstar_i <- cv_distances_transect(clust_k)
      clustgrid$w_stat[clustgrid$nk == nk] <- twosamples::wass_stat(stats::na.omit(gjstar_i), gij)
      clustgroups[[paste0("nk", nk)]] <- clust_k
    }
  }

  k_final <- clustgrid$nk[which.min(clustgrid$w_stat)]
  clust <- clustgroups[[paste0("nk", k_final)]]

  list(clusters = stats::setNames(clust, tpoints$id),
       w_stat = min(clustgrid$w_stat, na.rm = TRUE), q = k_final)
}

#' ----------------------------------------------------------------------
#' Try several candidate number of folds, k, and keep whichever gives the
#' lowest W statistic as per CAST::knndm(). 
#' ----------------------------------------------------------------------
choose_k <- function(tpoints, pt2transect, pred2transect, transect_ends,
                      k_candidates, maxp = 0.5, linkf = "ward.D2", nk_len = 100) {

  attempts <- k_candidates %>%
    purrr::map(function(k) {
      res <- tryCatch(
        knndm_transect(tpoints, pt2transect, pred2transect, transect_ends,
                        k = k, maxp = maxp, linkf = linkf, nk_len = nk_len),
        error = function(e) NULL
      )
      if (is.null(res)) return(NULL)
      tibble::tibble(k = k, w_stat = res$w_stat, q = res$q)
    }) %>%
    purrr::compact() %>%
    dplyr::bind_rows()

  list(table = attempts, best_k = attempts$k[which.min(attempts$w_stat)])
}

#' ----------------------------------------------------------------------
#' Build folds for one subgroup
#' ----------------------------------------------------------------------
build_subgroup_folds <- function(sex, region, period, min_dist, max_dist, label,
                                  crs = 5070, mcp_buffer_m = 5000, mcp_percentile = 100) {

  message("Building folds: ", label)

  in_period <- if (period == "early") {
    lubridate::year(pairs$date_1) < 2016 & lubridate::year(pairs$date_2) < 2016
  } else {
    lubridate::year(pairs$date_1) >= 2016 & lubridate::year(pairs$date_2) >= 2016
  }

  orig_index <- pairs %>%
    dplyr::mutate(.in_period = in_period) %>%
    dplyr::filter(.data$sex == .env$sex, .data$region == .env$region,
                  .data$euc_geog > .env$min_dist, .data$euc_geog <= .env$max_dist,
                  .data$.in_period) %>%
    dplyr::pull(index)

  if (length(orig_index) == 0) {
    warning("No eligible pairs for ", label)
    return(NULL)
  }

  pairs_sub <- pairs %>% dplyr::filter(index %in% orig_index)

  repeat_locs <- dplyr::bind_rows(
    dplyr::transmute(pairs_sub, id = id_1, x = long1, y = lat1),
    dplyr::transmute(pairs_sub, id = id_2, x = long2, y = lat2)
  ) %>%
    dplyr::distinct()

  anchors <- repeat_locs %>%
    dplyr::group_by(id) %>%
    dplyr::group_modify(~ tibble::as_tibble(get_medoid(as.matrix(.x[, c("x", "y")])))) %>%
    dplyr::ungroup()

  spread_by_animal <- repeat_locs %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(spread_m = if (dplyr::n() < 2) 0 else max(stats::dist(cbind(x, y))), .groups = "drop")

  tpoints <- anchors %>%
    sf::st_as_sf(coords = c("x", "y"), crs = crs, remove = FALSE)

  #' Model domain for this subgroup: MCP around every actual sample location
  #' in this sex x region x period subgroup (including repeated locations within
  #' animals), buffered by mcp_buffer_m. 
  model_domain <- repeat_locs %>%
    sf::st_as_sf(coords = c("x", "y"), crs = crs) %>%
    mcp(percentile = mcp_percentile) %>%
    sf::st_buffer(mcp_buffer_m) %>%
    sf::st_as_sf()

  predpoints <- sf::st_sample(model_domain, size = 3000, type = "regular")

  #' straight-line transect geometry from the endpoint coordinates
  transect_lines <- purrr::map(seq_len(nrow(pairs_sub)), function(i) {
    sf::st_linestring(rbind(c(pairs_sub$long1[i], pairs_sub$lat1[i]),
                             c(pairs_sub$long2[i], pairs_sub$lat2[i])))
  }) %>%
    sf::st_sfc(crs = crs)

  transect_ends <- cbind(match(pairs_sub$id_1, tpoints$id), match(pairs_sub$id_2, tpoints$id))

  pt2transect   <- sf::st_distance(tpoints, transect_lines)    %>% units::drop_units()
  pred2transect <- sf::st_distance(predpoints, transect_lines) %>% units::drop_units()

  n_animals <- nrow(tpoints)
  k_candidates <- c(4, 5, 6, 8, 10) %>% purrr::keep(~ .x <= max(2, floor(n_animals / 3)))
  if (length(k_candidates) == 0) k_candidates <- max(2, floor(n_animals / 3))

  k_search <- choose_k(tpoints, pt2transect, pred2transect, transect_ends, k_candidates)
  message("  k candidates tried: ", paste(k_search$table$k, collapse = ", "),
          " -- chose k = ", k_search$best_k, " (lowest W)")

  res <- knndm_transect(tpoints, pt2transect, pred2transect, transect_ends, k = k_search$best_k)
  clust <- res$clusters
  n_folds <- length(unique(clust))

  c1 <- clust[match(pairs_sub$id_1, names(clust))]
  c2 <- clust[match(pairs_sub$id_2, names(clust))]
  coin <- sample(c(1L, 2L), nrow(pairs_sub), replace = TRUE)
  test_cluster <- ifelse(c1 == c2, c1, ifelse(coin == 1L, c1, c2))

  #' create fold_indices to feed into caret::trainControl(index =, indexOut =) 
  fold_indices <- purrr::map(seq_len(n_folds), function(f) {
    touches_f <- (c1 == f) | (c2 == f)
    list(train = which(!touches_f), test = which(test_cluster == f))
  })

  list(
    index        = purrr::map(fold_indices, "train") %>% purrr::set_names(paste0("Fold", seq_len(n_folds))),
    indexOut     = purrr::map(fold_indices, "test")  %>% purrr::set_names(paste0("Fold", seq_len(n_folds))),
    orig_index   = orig_index,
    clusters     = clust,
    n_folds      = n_folds,
    w_stat       = res$w_stat,
    q            = res$q,
    k_search     = k_search$table,
    spread_by_animal = spread_by_animal,
    model_domain = model_domain,
    n_pairs      = nrow(pairs_sub),
    n_animals    = n_animals
  )
}

#' ----------------------------------------------------------------------
#' Run all 8 subgroups, save one .rds
#' ----------------------------------------------------------------------
subgroup_config <- tidyr::expand_grid(
  sex    = c("F", "M"),
  region = c("north", "south"),
  period = c("early", "late")
) %>%
  dplyr::mutate(
    min_dist = dplyr::if_else(sex == "F", 2500, 3500),
    max_dist = dplyr::if_else(sex == "F", 30000, 50000),
    label    = paste(dplyr::if_else(sex == "F", "female", "male"), region, period, sep = "_")
  )

fold_results <- subgroup_config %>%
  purrr::pmap(build_subgroup_folds) %>%
  purrr::set_names(subgroup_config$label)

dir.create("results", recursive = TRUE, showWarnings = FALSE)
saveRDS(fold_results, "results/fisher_cv_folds.rds")
