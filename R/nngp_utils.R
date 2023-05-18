
exp_cov <- function(dist, phi, r) {
    phi * exp(-dist / r)
}

euc_dist <- function(x1, y1, x2, y2) {
    sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

order_coords <- function(coords, space_id) {
    coords_dist <- euc_dist(coords[, "x"], 
                           coords[, "y"], 
                           max(coords[, "x"]), 
                           max(coords[, "y"]))
    coord_ordering <- order(coords_dist)
    coords_ordered <- coords[coord_ordering, ]
    coords_dist_ordered <- coords_dist[coord_ordering]
    space_id_ordered <- space_id[coord_ordering]

    coord_reverse_ordering <- match(space_id, space_id_ordered)

    coord_list <- list("ordered_coords" = coords_ordered,
                       "coord_ordering" = coord_ordering,
                       "coord_reverse_ordering" = coord_reverse_ordering)
    return(coord_list)
}

get_neighbors <- function(ordered_coords, m) {
    neighbors <- list()
    for (j in 1:(nrow(ordered_coords) - 1)) {
        dist_to_i <- euc_dist(ordered_coords[j, "x"], 
                              ordered_coords[j, "y"], 
                              ordered_coords[(j + 1):nrow(ordered_coords), "x"], 
                              ordered_coords[(j + 1):nrow(ordered_coords), "y"])  
        #get the m neighbor index wrt full coords
        neighbors[[j]] <- order(dist_to_i)[1:(min(m, length(dist_to_i)))] + j
        
    }
    neighbors[[nrow(ordered_coords)]] <- numeric(0)
    return(neighbors)
}

rnngp <- function(ordered_coords, neighbors, phi, r) {
    y <- rep(0, nrow(ordered_coords))
    y[nrow(ordered_coords)] <- stats::rnorm(1, 0, sqrt(exp_cov(0, phi, r)))
    for (i in (nrow(ordered_coords) - 1):1) {
        neighbor_i <- neighbors[[i]]
        neighbor_coords <- ordered_coords[neighbor_i, , drop = FALSE]
        neighbor_y <- y[neighbor_i]
        y_cov <- exp_cov(0, phi, r)
        cross_cov <- exp_cov(euc_dist(neighbor_coords[, "x"], neighbor_coords[, "y"], 
                                      ordered_coords[i, "x"], ordered_coords[i, "y"]), 
                             phi, r)
        neighbor_cov <- exp_cov(as.matrix(stats::dist(neighbor_coords, 
                                                      upper = TRUE, 
                                                      diag = TRUE)),
                                phi, r)
        cross_neighbor_cov <- t(cross_cov) %*% solve(neighbor_cov) 
        conditional_mean <- cross_neighbor_cov %*% neighbor_y
        conditional_cov <- y_cov - cross_neighbor_cov %*% cross_cov
        y[i] <- stats::rnorm(1, 
                             conditional_mean,
                             sqrt(conditional_cov))
    }
    return(y)
}



dnngp <- function(y, ordered_coords, neighbors, phi, r, log = FALSE) {
    d <- rep(0, nrow(ordered_coords))
    d[nrow(ordered_coords)] <- stats::dnorm(y[nrow(ordered_coords)], 
                                            0, 
                                            sqrt(exp_cov(0, phi, r)), 
                                            log = log)
    for (i in 1:(nrow(ordered_coords) - 1)) {
        neighbor_i <- neighbors[[i]]
        neighbor_coords <- ordered_coords[neighbor_i, , drop = FALSE]
        neighbor_y <- y[neighbor_i]
        y_cov <- exp_cov(0, phi, r)
        cross_cov <- exp_cov(euc_dist(neighbor_coords[, "x"], neighbor_coords[, "y"], 
                                      ordered_coords[i, "x"], ordered_coords[i, "y"]), 
                             phi, r)
        neighbor_cov <- exp_cov(as.matrix(stats::dist(neighbor_coords, 
                                                      upper = TRUE, 
                                                      diag = TRUE)),
                                phi, r)
        cross_neighbor_cov <- t(cross_cov) %*% solve(neighbor_cov) 
        conditional_mean <- cross_neighbor_cov %*% neighbor_y
        conditional_cov <- y_cov - cross_neighbor_cov %*% cross_cov
        d[i] <- stats::dnorm(y[i], 
                             conditional_mean,
                             sqrt(conditional_cov),
                             log = log)
    }
    return(d)
}
