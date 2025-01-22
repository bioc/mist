rmBad <- function(k, dat_ori, ptime_all) {
    dat <- dat_ori[k,]
    # remove 0 total count
    ind_remove <- is.nan(dat)
    dat <- dat[!ind_remove]
    ptime <- ptime_all[!ind_remove]
    segments <- seq(0, 1, by = 0.025) # timepoints are defined for every 0.02.
    indices <- lapply(seq_len(length(segments) - 1), function(i) {
        lower_bound <- segments[i]
        upper_bound <- segments[i + 1]
        if (upper_bound == 1) {
            which(ptime >= lower_bound & ptime <= upper_bound)
        } else {
            which(ptime >= lower_bound & ptime < upper_bound)
        }
    }) # get indices for each timepoints
    n_t <- lengths(indices)
    return(all(n_t %in% c(0, 1)) * 1)
}

# Function to check for Inf or NA
contains_inf_or_na <- function(x) {
    any(is.infinite(x) | is.na(x))
}

calculate_integral <- function(params_A, params_B = c(0, 0, 0, 0, 0, 0, 0, 0, 0)) {
    t_points <- seq(0, 1, length.out = 100) # Generate 100 points between 0 and 1
    differences <- numeric(length(t_points)) # Initialize a vector to store differences
    coef_A <- params_A[2:5] # Exclude intercept from curve A
    coef_B <- params_B[2:5] # Exclude intercept from curve B

    for (i in seq_len(length(t_points))) {
        t <- t_points[i]
        z_t <- c(t, t^2, t^3, t^4) # Polynomial terms, exclude '1' for the intercept term

        # Compute inverse logit transformed values without intercepts
        curve_A <- 1 / (1 + exp(-(sum(z_t * coef_A) + params_A[1])))
        curve_B <- 1 / (1 + exp(-(sum(z_t * coef_B) + params_B[1])))

        # Store the difference between the curves
        differences[i] <- curve_A - curve_B
    }

    # Calculate the optimal shift C
    C_optimal <- mean(differences)

    # Adjust the differences by subtracting C
    adjusted_differences <- abs(differences - C_optimal)
    # Compute the integral using the trapezoidal rule
    integral <- sum((adjusted_differences[-1] + adjusted_differences[-length(adjusted_differences)]) / 2 * diff(t_points))

    return(integral) # Return the computed integral
}


########
T_gamma <- function(shape, scale, lower, upper) {
    p_lower <- pgamma(lower, shape = shape, scale = scale)
    p_upper <- pgamma(upper, shape = shape, scale = scale)
    smallvalue <- 1e-08
    if (((p_lower > 1 - smallvalue) & (p_upper > 1 - smallvalue)) |
        ((p_lower < smallvalue) & (p_upper < smallvalue))) {
        return(lower * ((p_lower > 1 - smallvalue) & (p_upper > 1 - smallvalue)) +
            upper * ((p_lower < smallvalue) & (p_upper < smallvalue)))
    } else {
        random_p <- runif(1, p_lower, p_upper)

        return(qgamma(random_p, shape, scale = scale))
    }
}
####
# Function to run Bayesian estimation
run_bayesian_estimation <- function(k, dat_ready, ptime_all) {
    tryCatch(
        {
            T_gamma <- function(shape, scale, lower, upper) {
                p_lower <- pgamma(lower, shape = shape, scale = scale)
                p_upper <- pgamma(upper, shape = shape, scale = scale)
                smallvalue <- 1e-08
                if (((p_lower > 1 - smallvalue) & (p_upper > 1 - smallvalue)) |
                    ((p_lower < smallvalue) & (p_upper < smallvalue))) {
                    return(lower * ((p_lower > 1 - smallvalue) & (p_upper > 1 - smallvalue)) +
                        upper * ((p_lower < smallvalue) & (p_upper < smallvalue)))
                } else {
                    random_p <- runif(1, p_lower, p_upper)

                    return(qgamma(random_p, shape, scale = scale))
                }
            }
            poly_d <- 4
            stage_num <- 4
            Num_t <- 40
            beta_0 <- rnorm(1, mean = 0, sd = 3)
            poly_scale <- 0.025
            lambda2_0 <- 1
            tau2_0 <- 1
            sigma2 <- rep(2, stage_num)
            dat <- dat_ready[k, ]
            # remove 0 total count
            ind_remove <- is.nan(dat)
            dat <- dat[!ind_remove]
            ptime <- ptime_all[!ind_remove]
            ### correction
            dat[dat == 0] <- runif(length(dat[dat == 0]), 0.00001, 0.01)
            dat[dat == 1] <- runif(length(dat[dat == 1]), 0.95, 0.99999)

            segments <- seq(0, 1, by = 0.025) # timepoints are defined for every 0.02.
            indices <- lapply(seq_len(length(segments) - 1), function(i) {
                lower_bound <- segments[i]
                upper_bound <- segments[i + 1]
                if (upper_bound == 1) {
                    which(ptime >= lower_bound & ptime <= upper_bound)
                } else {
                    which(ptime >= lower_bound & ptime < upper_bound)
                }
            }) # get indices for each timepoints
            ### consider same values within a sigment to be one value
            dat_list <- lapply(indices, function(x) as.numeric(dat[x]))
            rm(indices, ptime)

            n_t <- lengths(dat_list)

            t_ex <- which(n_t %in% c(0, 1))
            z_t <- poly(poly_scale * (seq_len(Num_t)), poly_d, raw = TRUE, simple = TRUE) # zt
            u <- lapply(dat_list, function(x) {
                if (length(x) == 0) {
                    return(NULL)
                } else {
                    return(car::logit(x))
                }
            })
            u_eachT <- unlist(lapply(u, sum))
            stages <- (seq(Num_t) - 1) %/% (Num_t / stage_num) + 1
            shape_ig <- tapply(n_t, stages, sum) / 2
            rm(dat_list)

            lambda2 <- rep(lambda2_0, poly_d)
            tau2 <- tau2_0


            beta_mu_m <- matrix(nrow = 5000, ncol = poly_d)
            sigma2_m <- matrix(nrow = 5000, ncol = 4)
            beta0_v <- NULL
            for (i in seq_len(5000)) {
                u_eachT_beta0 <- unlist(lapply(u, function(x) (sum(x - beta_0))))
                s1 <- lapply(u, function(x) (x - beta_0))
                v_beta_mu <- chol2inv(chol(t(z_t) %*% diag(n_t) %*% diag(1 / rep(sigma2, each = Num_t / stage_num)) %*% z_t +
                    diag(1 / (lambda2 * tau2)))) ## V_beta_mu

                m_beta_mu <- v_beta_mu %*% colSums(diag(u_eachT_beta0) %*% (diag(1 / rep(sigma2, each = Num_t / stage_num)) %*% z_t)) ## m_beta_mu

                beta_mu <- as.vector(mvtnorm::rmvnorm(1, m_beta_mu, (t(v_beta_mu) + v_beta_mu) / 2))

                s2 <- z_t %*% beta_mu
                s3 <- mapply(function(s1, s2) s1 - s2, s1, s2, SIMPLIFY = FALSE)
                s4 <- vapply(s3, function(x) sum((x^2) / 2), numeric(1))
                scale_ig <- tapply(s4, stages, sum)

                sigma2 <- MCMCpack::rinvgamma(4, shape = shape_ig, scale = scale_ig)

                ### horseshoe for lambda
                eta <- 1 / lambda2
                upsi <- runif(poly_d, 0, 1 / (1 + eta))
                ub <- (1 - upsi) / upsi
                Fub <- 1 - exp(-(beta_mu^2 / (2 * tau2)) * ub)
                up <- runif(poly_d, 0, Fub)
                eta <- -log(1 - up) / (beta_mu^2 / (2 * tau2))
                lambda2 <- 1 / eta

                ### horseshoe for tau
                et <- 1 / tau2
                utau <- runif(1, 0, 1 / (1 + et))
                ubt <- (1 - utau) / utau
                scale_tau <- sum((beta_mu^2 / (2 * lambda2)))
                shape_tau <- (poly_d + 1) / 2
                upper_tau <- ubt
                et <- T_gamma(shape = shape_tau, scale = 1 / scale_tau, lower = 0, upper = upper_tau)
                tau2 <- 1 / et
                ### beta0
                v_beta0 <- 1 / ((1 / 9) + as.numeric(n_t %*% (1 / rep(sigma2, each = Num_t / stage_num))))
                m_beta0 <- v_beta0 * as.numeric(t(u_eachT - s2) %*% (1 / rep(sigma2, each = Num_t / stage_num))) ## m_beta_mu
                beta_0 <- rnorm(1, mean = m_beta0, sd = sqrt(v_beta0))
                ### Record Betas
                ### Record Betas
                beta_mu_m[i, ] <- m_beta_mu
                sigma2_m[i, ] <- sigma2
                beta0_v <- c(beta0_v, m_beta0)
            }


            beta_mu_mean <- apply(beta_mu_m[4000:5000, ], 2, mean)
            sigma2_mean <- apply(sigma2_m[4000:5000, ], 2, mean)
            beta0_mean <- mean(beta0_v[4000:5000])
            rm(beta_mu_m, sigma2_m, s1, s2, s3, z_t, beta0_v)
            return(c(beta0_mean, beta_mu_mean, sigma2_mean))
        },
        error = function(e) {
            # Return NA or another suitable value for all expected outputs
            return(rep(NA, 9)) # Update `length_of_expected_output` accordingly
        }
    )
}
