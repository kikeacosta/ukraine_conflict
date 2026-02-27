options(scipen = 9999)
try(dev.off(), silent = T)

# installing and loading required packages ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# install pacman to streamline further package installation
if (!require("pacman", character.only = TRUE)) {
  install.packages("pacman", dep = TRUE)
  if (!require("pacman", character.only = TRUE)) {
    stop("Package pacman not found")
  }
}
library(pacman)

packages_CRAN <- c(
  "tidyverse",
  "countrycode",
  "lubridate",
  "readxl",
  "ungroup",
  "mgcv",
  "data.table",
  "MortalityLaws",
  "purrr",
  "demography",
  "patchwork",
  "fields",
  "forecast",
  "ggrepel",
  "ggh4x"
)

# Install required CRAN packages if not available yet
if (sum(!p_isinstalled(packages_CRAN)) > 0) {
  install_these <- packages_CRAN[!p_isinstalled(packages_CRAN)]
  for (i in 1:length(install_these)) {
    install.packages(install_these[i], dependencies = "Depends")
  }
}

# Load the required CRAN/github packages
p_load(packages_CRAN, character.only = TRUE)

copy_this <- function(x, row.names = FALSE, col.names = TRUE, ...) {
  write.table(
    x,
    file = paste0("clipboard-", object.size(x)),
    sep = "\t",
    row.names = row.names,
    col.names = col.names,
    ...
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# functions for processing data ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # function to split CDR into mx for mortality crises
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dts <- ucdp2
# pop <- exps
# yrs <- 2023
crd_to_mx2 <- function(dts_in, pop_in, as) {
  ev_type <- c("Conflict", "Conflict: combatant")

  # adjusting UN patterns
  # no mortality among combatants for children under 10
  as2 <-
    as %>%
    filter(eventype %in% ev_type, sex != "t") %>%
    mutate(
      rr = ifelse(eventype == "Conflict: combatant" & age < 10, 0, rr),
      rr_ll = ifelse(eventype == "Conflict: combatant" & age < 10, 0, rr_ll),
      rr_ul = ifelse(eventype == "Conflict: combatant" & age < 10, 0, rr_ul)
    )

  # assuming civilians as deaths in "conflict" scenario
  # and combatants in the "Conflict: combatant" scenario
  # first assuming a CDR of 1/K

  conf <-
    dts_in %>%
    mutate(
      eventype = case_when(
        role == "civilians" ~ "Conflict",
        role == "combatants" ~ "Conflict: combatant"
      )
    ) %>%
    left_join(as2, by = "eventype", relationship = "many-to-many") %>%
    left_join(pop_in, by = join_by(year, sex, age)) %>%
    mutate(
      dx_sd = rr * exposure / 1000,
      dx_sd_ll = rr_ll * exposure / 1000,
      dx_sd_ul = rr_ul * exposure / 1000
    )

  # adjustment of death rates to match the total
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dts3 <-
    conf %>%
    group_by(year, eventype) %>%
    mutate(dx_sd_sum = sum(dx_sd)) %>%
    ungroup() %>%
    mutate(
      dx = dx_sd * dts / dx_sd_sum,
      dx_ll = dx_sd_ll * dts / dx_sd_sum,
      dx_ul = dx_sd_ul * dts / dx_sd_sum
    ) %>%
    select(
      everything(),
      -starts_with("rr"),
      -starts_with("dx_sd"),
      -starts_with("inc")
    )

  dts4 <-
    dts3 %>%
    bind_rows(
      dts3 %>%
        reframe(
          dts = sum(dts),
          exposure = mean(exposure),
          dx = sum(dx),
          dx_ll = sum(dx_ll),
          dx_ul = sum(dx_ul),
          .by = c(year, sex, age)
        ) %>%
        # ungroup() %>%
        mutate(role = "all", eventype = "all")
    ) %>%
    mutate(
      mx = dx / exposure,
      mx_ll = dx_ll / exposure,
      mx_ul = dx_ul / exposure
    ) %>%
    select(year, role, sex, age, everything()) %>%
    arrange(year, role, sex, age)

  return(dts4)
}

# ungrouping ages
ung_age_x <- function(chunk) {
  dt_in <-
    tibble(age = chunk$age, dx = chunk$dx, dx_mt = dx * 1e5) %>%
    mutate(dx_mt = ifelse(dx_mt == 0, 1, dx_mt))
  nl <- 26
  dxs <- pclm(x = dt_in$age, y = dt_in$dx_mt, nlast = nl)$fitted
  fit <- tibble(age = 0:100, dx = dxs / 1e5)

  out <-
    chunk %>%
    select(year, sex) %>%
    unique() %>%
    left_join(fit, by = character())
  return(out)
}

# Funtions for estimating exposures and conflict mortality rates ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

est_exposures <- function(chunk) {
  chunk |>
    left_join(migs2, by = join_by(sex, age, year)) %>%
    left_join(exp_mort, by = join_by(sex, age, year)) %>%
    left_join(asfr2, by = join_by(sex, age, year)) %>%
    replace_na(list(fx = 0)) %>%
    mutate(
      pop2 = pop - 0.5 * ems,
      dx = pop2 * mx,
      exposure = pop2 - 0.5 * dx,
      pop = ifelse(age == 0, 0.5 * sum(fx * exposure), pop),
      pop2 = ifelse(age == 0, pop - 0.5 * ems, pop2),
      dx = ifelse(age == 0, pop2 * mx, dx),
      exposure = ifelse(age == 0, pop2 - 0.5 * dx, exposure)
    ) %>%
    select(year, sex, age, ems, dx, exposure, pop_ini = pop)
}

est_age_sex_conf <- function(exposures, cnf, rl) {
  pop_in <-
    exposures |>
    mutate(
      age = case_when(
        age == 0 ~ 0,
        age %in% 1:4 ~ 1,
        age >= 75 ~ 75,
        TRUE ~ age - age %% 5
      )
    ) %>%
    reframe(exposure = sum(exposure), .by = c(year, sex, age))

  yr <- unique(exposures$year)

  cnf_in <-
    cnf %>%
    filter(year == yr)

  cnf_age_sex <- crd_to_mx2(cnf_in, pop_in, as)

  cnf_age_sex_ung <-
    cnf_age_sex %>%
    group_by(role, year, sex) %>%
    do(ung_age_x(chunk = .data)) %>%
    ungroup() %>%
    rename(cnf = dx) %>%
    filter(role == rl) %>%
    select(-role)

  return(cnf_age_sex_ung)
}

est_pop_end <- function(exposures, cnf_y) {
  exposures %>%
    left_join(cnf_y, by = join_by(year, sex, age)) %>%
    # removing net migration, expected deaths, and conflict deaths
    mutate(pop_end = pop_ini - ems - dx - cnf)
}

est_summary <- function(pop_end) {
  sum <-
    pop_end %>%
    select(year, sex, age, conflict = cnf, expected = dx, pop = exposure) %>%
    mutate(all = conflict + expected) %>%
    mutate(across(c(all, conflict, expected), as.numeric)) %>%
    pivot_longer(
      c(all, conflict, expected),
      names_to = "cause",
      values_to = "dx"
    ) %>%
    mutate(mx = dx / pop)

  return(sum)
}

end_to_ini <- function(pop_end) {
  pop_end %>%
    select(year, sex, age, pop = pop_end) %>%
    mutate(age = ifelse(age == 100, 100, age + 1), year = year + 1) %>%
    summarise(pop = sum(pop), .by = c(year, sex, age)) %>%
    bind_rows(tibble(age = 0, sex = c("f", "m"), pop = 0)) |>
    fill(year) |>
    arrange(sex, age)
}

# function ofr estimating a life table ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lifetable <- function(dt_in) {
  x <- dt_in$age
  mx <- dt_in$mx
  sex <- unique(dt_in$sex)
  year <- unique(dt_in$year)

  m <- length(x)
  n <- c(diff(x), NA)
  ax <- rep(0, m)
  if (x[1] != 0 | x[2] != 1) {
    ax <- n / 2
    ax[m] <- 1 / mx[m]
  } else {
    if (sex == "f") {
      if (mx[1] < 0.01724) {
        ax[1] <- 0.14903 - 2.05527 * mx[1]
      } else if (mx[1] >= 0.01724 & mx[1] < 0.06891) {
        ax[1] <- 0.04667 + 3.88089 * mx[1]
      } else {
        ax[1] <- 0.31411
      }
    }
    if (sex == "m") {
      if (mx[1] < 0.02300) {
        ax[1] <- 0.14929 - 1.99545 * mx[1]
      } else if (mx[1] >= 0.02300 & mx[1] < 0.08307) {
        ax[1] <- 0.02832 + 3.26021 * mx[1]
      } else {
        ax[1] <- 0.29915
      }
    }
    ax[-1] <- n[-1] / 2
    ax[m] <- 1 / mx[m]
  }
  qx <- n * mx / (1 + (n - ax) * mx)
  qx[m] <- 1
  px <- 1 - qx
  lx <- cumprod(c(1, px)) * 100000
  dx <- -diff(lx)
  Lx <- n * lx[-1] + ax * dx
  lx <- lx[-(m + 1)]
  Lx[m] <- lx[m] / mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  dt_out <- data.frame(age = x, mx, ex)
  return(dt_out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Funtions for forecasting mortality rates ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

frcst <- function(data, t = 2012:2019) {
  x <- unique(data$age)
  m <- length(x)

  n <- length(t)
  t.fore <- (t[n] + 1):2025
  n.fore <- length(t.fore)
  ## LC routine ------

  ## extract deaths and exposures of interest
  y <- data %>%
    # filter(sex==my.sex,region==my.reg,source==my.source) %>%
    filter(year %in% t) %>%
    select(dts) %>%
    pull()
  e <- data %>%
    # filter(sex==my.sex,region==my.reg,source==my.source) %>%
    filter(year %in% t) %>%
    select(pop) %>%
    pull()
  Y <- matrix(y, m, n)
  E <- matrix(e, m, n)
  LMX <- log(Y / E)

  # ## check evolution of e0
  # e0 <- apply(exp(LMX), 2, e0.mx, x=x, sex="F")
  # # plot(t,e0)

  ## plotting
  # matplot(x,LMX,t="l",lty=1,col=rainbow(n))
  # matplot(t,t(LMX),t="l",lty=1,col=rainbow(m))

  ## starting LC parameters
  start.pars <- LC_starting_pars(Dth = Y, Exp = E)
  Alpha <- start.pars$Alpha
  Beta <- start.pars$Beta
  Kappa <- start.pars$Kappa
  One <- start.pars$One
  Eta <- Alpha %*% t(One) + Beta %*% t(Kappa)
  ## compute fitted deaths
  D.fit <- E * exp(Eta)

  ## fit LC
  ## starting the iteration
  for (iter in 1:1000) {
    Alpha.old <- Alpha
    Beta.old <- Beta
    Kappa.old <- Kappa
    ## update Alpha
    temp <- Update.alpha(Alpha, Beta, Kappa, One, Dth = Y, Exp = E, D.fit)
    D.fit <- temp$D.fit
    Alpha <- temp$Alpha
    ## update Beta
    temp <- Update.beta(Alpha, Beta, Kappa, One, Dth = Y, Exp = E, D.fit)
    D.fit <- temp$D.fit
    Beta <- temp$Beta
    ## update Kappa
    temp <- Update.kappa(Alpha, Beta, Kappa, One, Dth = Y, Exp = E, D.fit)

    D.fit <- temp$D.fit
    Kappa <- temp$Kappa
    ## tolerance criterion
    crit <- max(
      max(abs(Alpha - Alpha.old)),
      max(abs(Beta - Beta.old)),
      max(abs(Kappa - Kappa.old))
    )
    # cat(iter,crit,"\n")
    if (crit < 1e-04) break
  }

  ## adding the constraints
  sum.Beta <- sum(Beta)
  Beta <- Beta / sum.Beta
  Kappa <- Kappa * sum.Beta

  ## plotting LC parameters
  # par(mfrow=c(1,3))
  # plot(x,Alpha)
  # plot(x,Beta);abline(h=0)
  # plot(t,Kappa);abline(h=0)
  # par(mfrow=c(1,1))

  ## forecasting kappa index

  ## set a RW model with drift
  modK <- Arima(
    ts(Kappa, start = t[1]),
    order = c(0, 1, 0),
    include.drift = TRUE
  )
  ## forecast Kappa until 2040
  predK <- forecast(modK, h = n.fore)
  ## plotting the time-series and forecast
  # plot(predK)
  ## forecast Kappa
  KappaF <- predK$mean
  ## forecast Eta
  OneF <- matrix(1, nrow = n.fore, ncol = 1)
  LMX.fore <- Alpha %*% t(OneF) + Beta %*% t(KappaF)

  ## plotting
  # cols <- viridis(n+n.fore)
  # matplot(x,LMX,t="l",lty=1,col=cols[1:n])
  # matlines(x,LMX.fore,lty=1,col=cols[1:n.fore+n])

  # ## life expectancy
  # e0.fore <- apply(exp(LMX.fore),2, e0.mx, x=x, sex="F")

  ## plotting e0
  # plot(t,e0,ylim=range(e0,e0.fore),xlim=range(t,t.fore))
  # points(t.fore,e0.fore,col=4,lwd=2,pch=16)

  ## plotting rates
  # whi.age <- 80
  # plot(t,LMX[which(x==whi.age),],xlim=range(t,t.fore),
  #      ylim=range(LMX[which(x==whi.age),],LMX.fore[which(x==whi.age),]))
  # points(t.fore,LMX.fore[which(x==whi.age),],col=4,lwd=2,pch=16)

  my.sex <- unique(data$sex)
  # my.reg <- unique(data$region)
  # my.source <- unique(data$source)

  ## results in tibble
  data.with.fore <- tibble(
    year = rep(t.fore, each = m),
    sex = my.sex,
    age = rep(x, n.fore),
    # region=my.reg,
    # source=paste0("lc_",my.source,"_",t[n]),
    mx = c(1e5 * exp(LMX.fore))
  )

  ## adding results to data
  return(data.with.fore)
}


## functions to fit Poisson LC model

## Starting values of LC model
LC_starting_pars <- function(Dth, Exp) {
  ## dimensions
  m <- nrow(Dth)
  n <- ncol(Dth)
  ## starting values
  One <- matrix(1, nrow = n, ncol = 1)
  Fit <- log((Dth + 1) / (Exp + 2))
  ## for Alpha, take mean of log death rates
  Alpha <- apply(Fit / n, 1, sum, na.rm = T)
  ## for Beta, take alpha and normalize to sum to one
  Beta <- matrix(1 * Alpha, ncol = 1)
  sum.Beta <- sum(Beta)
  Beta <- Beta / sum.Beta
  ## for Kappa, standardize a series from n to 1
  Kappa <- matrix(seq(n, 1, by = -1), nrow = n, ncol = 1)
  Kappa <- Kappa - mean(Kappa)
  Kappa <- Kappa / sqrt(sum(Kappa * Kappa))
  ## return list
  out <- list(Alpha = Alpha, Beta = Beta, Kappa = Kappa, One = One)
}

## Update Alpha
Update.alpha <- function(Alpha, Beta, Kappa, One, Dth, Exp, D.fit) {
  difD <- Dth - D.fit
  Alpha <- Alpha +
    difD %*% One / ifelse((D.fit %*% One) == 0, 1e-8, D.fit %*% One)
  Eta <- Alpha %*% t(One) + Beta %*% t(Kappa)
  D.fit <- Exp * exp(Eta)
  list(Alpha = Alpha, D.fit = D.fit)
}
## Update Beta
Update.beta <- function(Alpha, Beta, Kappa, One, Dth, Exp, D.fit) {
  difD <- Dth - D.fit
  Kappa2 <- Kappa * Kappa
  Beta <- Beta +
    difD %*% Kappa / ifelse((D.fit %*% Kappa2) == 0, 1e-8, D.fit %*% Kappa2)
  Eta <- Alpha %*% t(One) + Beta %*% t(Kappa)
  D.fit <- Exp * exp(Eta)
  list(Beta = Beta, D.fit = D.fit)
}
## Update Kappa
Update.kappa <- function(Alpha, Beta, Kappa, One, Dth, Exp, D.fit) {
  difD <- Dth - D.fit
  Beta2 <- Beta * Beta
  Kappa <- Kappa + t(difD) %*% Beta / (t(D.fit) %*% Beta2)
  Kappa <- Kappa - mean(Kappa)
  Kappa <- Kappa / sqrt(sum(Kappa * Kappa))
  Kappa <- matrix(Kappa, ncol = 1)
  Eta <- Alpha %*% t(One) + Beta %*% t(Kappa)
  D.fit <- Exp * exp(Eta)
  list(Kappa = Kappa, D.fit = D.fit)
}
## function for constructing a classic (& rather general) lifetable
## from mortality rates
lifetable.mx <- function(x, mx, sex = "m", ax = NULL) {
  m <- length(x)
  n <- c(diff(x), NA)
  if (is.null(ax)) {
    ax <- rep(0, m)
    if (x[1] != 0 | x[2] != 1) {
      ax <- n / 2
      ax[m] <- 1 / mx[m]
    } else {
      if (sex == "f") {
        if (mx[1] >= 0.107) {
          ax[1] <- 0.350
        } else {
          ax[1] <- 0.053 + 2.800 * mx[1]
        }
      }
      if (sex == "m") {
        if (mx[1] >= 0.107) {
          ax[1] <- 0.330
        } else {
          ax[1] <- 0.045 + 2.684 * mx[1]
        }
      }
      ax[-1] <- n[-1] / 2
      ax[m] <- 1 / mx[m]
    }
  }
  qx <- n * mx / (1 + (n - ax) * mx)
  qx[m] <- 1
  px <- 1 - qx
  lx <- cumprod(c(1, px)) * 100000
  dx <- -diff(lx)
  Lx <- n * lx[-1] + ax * dx
  lx <- lx[-(m + 1)]
  Lx[m] <- lx[m] / mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  return.df <- data.frame(x, n, mx, ax, qx, px, lx, dx, Lx, Tx, ex)
  return(return.df)
}

## function to derive e0 from mx (based on previous function)
e0.mx <- function(x, mx, sex = "m", ax = NULL) {
  lt <- lifetable.mx(x, mx, sex, ax)
  return.ex <- lt$ex[1]
  return(return.ex)
}
