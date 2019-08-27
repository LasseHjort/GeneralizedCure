
#Function for calculting the cure rate log-odds ratio (for explicit cure models)
get_lor <- function(object, newdata, var.link = function(x) x, X.cr = X.cr) {
  gamma <- object$coefs
  link.type.cr <- object$link.type.cr
  eta_pi <- as.vector(X.cr %*% gamma)
  pi <- cuRe:::get.link(link.type.cr)(eta_pi)
  odds <- pi / (1 - pi)
  est <- log(odds[1] / odds[2])
  return(est)
}

#Function for calculting the cure rate log-odds ratio (for latent cure models)
get_lor2 <- function(object, newdata, var.link = function(x) x,
                     X = X) {
  link.surv <- object@link
  eta <- as.vector(X %*% object@fullcoef)
  
  pi <- link.surv$ilink(eta)
  odds <- pi / (1 - pi)
  est <- log(odds[1] / odds[2])
  return(est)
}

#Function for calculating standard errors by using the delta method and numerical differentiation
numDeltaMethod2 <- function (object, fun, gd = NULL, ...)
{
  coef <- object@fullcoef
  est <- fun(coef, ...)
  Sigma <- object@vcov
  if (is.null(gd))
    gd <- rstpm2:::grad(fun, coef, ...)
  se.est <- as.vector(sqrt(colSums(gd * (Sigma %*% gd))))
  data.frame(Estimate = est, SE = se.est)
}

#Wrapper function to obtain both LOR estimates and corresponding standard errors
predictnl.default2 <- function (object, fun, newdata = NULL, gd = NULL, ...)
{
  if (is.null(newdata) && !is.null(object$data))
    newdata <- object$data
  localf <- function(coef, ...) {
    object@fullcoef <- coef
    fun(object, ...)
  }
  numDeltaMethod2(object, localf, newdata = newdata, gd = gd,
                  ...)
}



#Function for simulating survival data including a binary covariate effect
sim_surv_cov <- function(pars, age, n, type = "weibull", link = "logit"){
  #Generate general population survival. As implemented, age, gender, and year are fixed.
  sim_ages <- rnorm(n, mean = age, sd = 0)
  sim_gender <- factor(rbinom(n, size = 1, prob = 1), levels = 0:1, labels = c("male", "female"))
  dates <- as.Date("1980-01-01") #+ 1:as.numeric(as.Date("1980-01-01") - as.Date("1980-01-01"))
  diag_dates <- sample(dates, n, replace = T)
  D <- data.frame(age = sim_ages * ayear, sex = sim_gender, diag_date = diag_dates)
  S_exp <- survexp(~1, rmap = list(age = age, sex = sex, year = diag_date), 
                   data = D,
                   times = seq(0, (110 - age) * ayear, length.out = 1000),
                   ratetable = survexp.dk, 
                   scale = ayear)
  
  #Function to draw simulations from
  len <- length(S_exp$time)
  sim_fun_gen <- function(x){
    if(x > 1 - S_exp$surv[len]){
      S_exp$time[len]
    }else{
      S_exp$time[findInterval(x, vec = 1 - S_exp$surv) + 1] 
    }
  }
  
  #Simulate binary covariate and calculate pi
  z <- rbinom(n = n, size = 1, prob = 0.5)
  M <- model.matrix(~z)
  pi <- as.numeric(links[[link]](M %*% pars[1:2]))
  lp <- as.numeric(M %*% pars[3:4])
  
  #Set survival of the uncured and relative survival functions using the Weibull distribution
  if(type == "weibull"){
    unifind <- function(time, i, u){
      surv <- exp(-exp(lp[i]) * time ^ pars[5])
      surv - u
    }
  }
  
  #Function to simulate data
  sim_fun_rel <- function(x){
    res <- rep(Inf, length(x))
    for(i in 1:length(x)){
      #cat(i, "\n")
      if(x[i] >= pi[i]){
        res[i] <- uniroot(unifind, interval = c(0, 150), u = (x[i] - pi[i]) /(1 - pi[i]), i = i)$root
      }
    }
    res
  }
  
  #Simulate uniform variable for general population and disease specific survival
  uni_sim1 <- runif(nrow(D))
  uni_sim2 <- runif(nrow(D))
  
  #Simulate from both distributions
  sim_gen <- sapply(uni_sim1, sim_fun_gen)
  sim_pop <- sim_fun_rel(uni_sim2)
  D$fu <- pmin(sim_gen, sim_pop)
  
  #Simulate from censoring distribution
  #Set parameters
  max <- 15
  sim_cens <- runif(n = nrow(D), min = 0, max = max)
  
  #Generated follow-up as the minimum of survival time and censoring time
  D$FU <- pmin(D$fu, sim_cens)
  D$status <- as.numeric(D$fu <= sim_cens)
  D$FU[D$FU < 1e-3] <- 1e-3
  #Check the imperical survival function
  # plot(survfit(Surv(FU, status) ~ 1, data = D))
  # curve(rel_surv, from = 0, to = 15, add = T, col = 2)
  #Follow-up in days
  D$FU <- D$FU *  ayear
  #Check the imperical relative survival function
  # rsfit <- rs.surv(Surv(FU, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_date), 
  #                  data = D, ratetable = survexp.dk, method = "ederer2")
  # 
  # plot(rsfit)
  # abline(h = 0.5)
  D$FU_years <- D$FU / ayear
  #Get general population hazard
  D$exp_haz <- general.haz(time = "FU", age = "age", sex = "sex", year = "diag_date", 
                           data = D, ratetable = survexp.dk)
  D$z <- z
  list(D = D, pars = pars, pi = pi, lp = lp, link = link, S_exp = S_exp)
}


#Create lists for easy use fo different link functions
links <- list(logit = function(x) exp(x) / (exp(x) + 1), 
              loglog = function(x) exp(-exp(x)))

inv.links <- list(logit = function(x) log(x / (1 - x)), 
                  loglog = function(x) log(-log(x)))

links.math <- list(logit = "$\\log\\left(\\frac{x}{1-x}\\right)$",
                   loglog = "$\\log[-\\log(x)]$")

#Specify cases and save to data frame
age <- 60
pis <- c(0.25, 0.5, 0.6)
link_S <- c("loglog", "loglog", "loglog")
link <- c("logit", "loglog", "logit")
beta_0 <- sapply(1:length(pis), function(i) inv.links[[link[i]]](pis[i]))
beta_1 <- c(0.5, 0.5, 0.5)
exp_theta_0 <- c(1, 0.1, 0.5)
theta_1 <- c(0, 0.5, 0.5) 
theta_2 <- c(1, 1.4, 1)
age <- rep(age, length(pis))
n <- rep(1000, length(pis))

D <- data.frame(pi = pis, beta_0 = beta_0, beta_1 = beta_1, theta_0 = log(exp_theta_0), 
                theta_1 = theta_1, theta_2 = theta_2)

#Specify simulations
cases <- vector("list", nrow(D))
for(i in 1:length(cases)){
  cases[[i]]$pars <- as.numeric(D[i,-1])
  cases[[i]]$link <- link[i]
  cases[[i]]$age <- age[i]
  cases[[i]]$n <- n[i]
  cases[[i]]$type <- "weibull"
}

#Function for calculating the log odds ratio
calcLogOR <- function(case){
  lp <- case$pars[1]
  lp_treat <- case$pars[1] + case$pars[2]
  prob <- links[[case$link]](lp)
  prob_treat <- links[[case$link]](lp_treat)
  log((prob_treat / (1 - prob_treat)) / (prob / (1 - prob)))
}

#Calculate the log true odds ratio in each scenario
trueor <- rep(NA, length(pis))
for(i in 1:length(trueor)){
  trueor[i] <- calcLogOR(cases[[i]])
}

#Create table with parameters of each case
D$link_S <- unlist(links.math[link_S])
D$link <- unlist(links.math[link])
D$LOR <- trueor
names(D) <- c("$\\pi$", "$\\beta_0$", "$\\beta_1$", "$\\theta_0$", "$\\theta_1$", 
              "$\\theta_2$", "$g$", "$g_{\\pi}$", "LOR")
is.num <- sapply(D, is.numeric)
D[is.num] <- lapply(D[is.num], round, 3)
D <- cbind(Scenario = 1:nrow(D), D)

print(xtable(D, align = "lcccccccccc", 
             caption = "Parameter values used for simulations including covariate effects. 
             LOR: true log-odds ratio for the cure proportion.", 
             label = "tab:paramval"), 
      file = file.path(tab.out, "ParamValues.tex"), include.rownames = F,
      sanitize.colnames.function = identity, sanitize.text.function = identity)

#Check relative survival trajectory
# plot_case <- function(case){
#   time <- seq(0, 30, length.out = 100)
#   z <- 0
#   pi <- links[[case$link]](case$pars[1] + case$pars[2] * z)
#   surv <- exp(-exp(case$pars[3] + case$pars[4] * z + case$pars[5] * log(time)))
#   rsurv <- pi + (1 - pi) * surv
#   plot(rsurv ~ time, type = "l", ylim = c(0, 1))
#   z <- 1
#   pi <- links[[case$link]](case$pars[1] + case$pars[2] * z)
#   surv <- exp(-exp(case$pars[3] + case$pars[4] * z + case$pars[5] * log(time)))
#   rsurv <- pi + (1 - pi) * surv
#   lines(rsurv ~ time, col = 2)
# }

# plot_case(cases[[1]])
# plot_case(cases[[2]])
# plot_case(cases[[3]])

#Set number of nodes for the Gauss-Legendre quadrature
gaussxw <- statmod::gauss.quad(100)
weights <- gaussxw$weights

#Set upper limit of integral
upper.lim <- 15
time <- upper.lim / 2 * gaussxw$nodes + upper.lim / 2

#Select relevant values of z
z_values <- c(0, 1)

#Set simulation parameters
n.models <- 5
n.sim <- 500
n.cores <- 48

#Run simulations
filename <- file.path(data.out, "sim_res_cov.RData")
if(file.exists(filename)){
  load(filename)
} else {
  set.seed(20180625, "L'Ecuyer")
  L <- vector("list", length(cases))
  for(i in 1:length(cases)){
    cat(i, "\n")
    #Calculate the true pi, survival of the uncured, and relative survival
    pi_0_fun <- function(z) links[[cases[[i]]$link]](D[i, 3] + D[i, 4] * z)
    survu_0_fun <- function(z) exp(-exp(D[i,5] + D[i, 6] * z + D[i, 7] * log(time)))
    surv_0_fun <- function(z) pi_0_fun(z) + (1 - pi_0_fun(z)) * survu_0_fun(z)
    pi_0 <- sapply(z_values, pi_0_fun)
    survu_0 <- sapply(z_values, survu_0_fun)
    surv_0 <- sapply(z_values, surv_0_fun)
    #logor <- var <- int <- diae <- diae_surv <- matrix(ncol = 5, nrow = n.sim)
    
    res <- mclapply(1:n.sim, function(j){
      logor_j <- var_j <- int_j <- diae_j <- diae_surv_j <- rep(NA, n.models)
      #cat(j, "\n")
      sim_data <- do.call(sim_surv_cov, cases[[i]])
      
      #Check imperical relative survival trajectory
      # fit <- rs.surv(Surv(FU, status) ~ z + ratetable(age = age, sex = sex, year = diag_date),
      #               ratetable = survexp.dk, data = sim_data$D, method = "ederer2")
      # plot(fit, col = 1:2)
      # plot_case(cases[[i]])
      
      #Infinite cure point models
      #Weibull mixture cure model
      fit <- fit.cure.model(Surv(FU_years, status) ~ z, formula.surv = list(~ z, ~ 1), 
                            data = sim_data$D, bhazard = "exp_haz")
      #Calculate odds ratio and CI
      logor_j[1] <- fit$coefs[[1]][2]
      if(!is.null(fit$covariance)){
        var_j[1] <- sqrt(fit$covariance[2,2]) 
      }
      
      #Get absolute error
      pred_z <- predict(fit, newdata = data.frame(z = z_values), type = "curerate", var.type = "n")
      int_j[1] <- mean(abs(do.call(rbind, pred_z)$Estimate - pi_0))
      
      #Get integrated absolute error
      surv_u <- predict(fit, time = time, 
                        newdata = data.frame(z = z_values), type = "survuncured", var.type = "n")
      surv_u <- do.call(cbind, surv_u)
      diae_j[1] <- upper.lim / 2 * sum(weights * rowMeans(abs(surv_u - survu_0)))
      
      #Get integrated absolute error - relative survival
      surv <- do.call(cbind, predict(fit, time = time, 
                                     newdata = data.frame(z = z_values), var.type = "n"))
      diae_surv_j[1] <- upper.lim / 2 * sum(weights * rowMeans(abs(surv - surv_0)))
      
      #Mixture proportional hazards with logit link for pi
      fit1 <- GenFlexCureModel(Surv(FU_years, status) ~ z, 
                               cr.formula = ~ z,
                               data = sim_data$D,
                               bhazard = "exp_haz",
                               df = 4, verbose = F)
      logor_j[2] <- fit1$coefs[2]
      if(!is.null(fit1$covariance)){
        var_j[2] <- sqrt(fit1$covariance[2,2])
      }
      #Get absolute error
      pred <- predict(fit1, newdata = data.frame(z = z_values), type = "curerate")
      pred <- do.call(rbind, pred)$Estimate
      int_j[2] <- mean(abs(pred - pi_0))
      
      #Get integrated absolute error
      surv_u <- predict(fit1, time = time, 
                        newdata = data.frame(z = z_values), type = "survuncured", var.type = "n")
      surv_u <- do.call(cbind, surv_u)
      diae_j[2] <- upper.lim / 2 * sum(weights * rowMeans(abs(surv_u - survu_0)))
      
      #Get integrated absolute error - relative survival
      surv <- do.call(cbind, predict(fit1, time = time, 
                                     newdata = data.frame(z = z_values), var.type = "n"))
      diae_surv_j[2] <- upper.lim / 2 * sum(weights * rowMeans(abs(surv - surv_0)))
      
      
      #Mixture proportional hazards with loglog link for pi
      fit1 <- GenFlexCureModel(Surv(FU_years, status) ~ z, 
                               cr.formula = ~ z,
                               data = sim_data$D,
                               bhazard = "exp_haz", 
                               link.type.cr = "loglog",
                               df = 4, verbose = F)
      
      #Get absolute error
      X.cr <- model.matrix( ~ z, data = data.frame(z = c(1, 0)))
      if(is.null(fit1$covariance)){
        logor_j[3] <- get_lor(fit1, newdata = NULL, X.cr = X.cr)
      } else {
        est <- cuRe:::predictnl.default(fit1, get_lor, X.cr = X.cr)
        logor_j[3] <- est$Estimate 
        var_j[3] <- est$SE
      }
      pred <- predict(fit1, newdata = data.frame(z = z_values), type = "curerate")
      pred <- do.call(rbind, pred)$Estimate
      int_j[3] <- mean(abs(pred - pi_0))
      
      #Get integrated absolute error
      surv_u <- predict(fit1, time = time, 
                        newdata = data.frame(z = z_values), type = "survuncured", var.type = "n")
      surv_u <- do.call(cbind, surv_u)
      diae_j[3] <- upper.lim / 2 * sum(weights * rowMeans(abs(surv_u - survu_0)))
      
      #Get integrated absolute error - relative survival
      surv <- do.call(cbind, predict(fit1, time = time, 
                                     newdata = data.frame(z = z_values), var.type = "n"))
      diae_surv_j[3] <- upper.lim / 2 * sum(weights * rowMeans(abs(surv - surv_0)))
      
      
      #Finite cure point models
      
      #Proportional hazards
      fit <- stpm2(Surv(FU_years, status) ~ z, data = sim_data$D, bhazard = sim_data$D$exp_haz, 
                   df = 5, cure = T, tvc = list(z = 1))
      X <- rstpm2:::lpmatrix.lm(fit@lm, newdata = data.frame(z = c(1, 0), FU_years = max(sim_data$D$FU_years)))
      
      #Get absolute error 
      if(is.null(fit@vcov)){
        logor_j[4] <- get_lor2(fit, newdata = NULL, X = X)
      } else {
        est <- predictnl.default2(fit, fun = get_lor2, X = X, newdata = X)
        logor_j[4] <- est$Estimate 
        var_j[4] <- est$SE
      }
      pred <- predict(fit, newdata = data.frame(z = z_values, FU_years = max(sim_data$D$FU_years)), 
                      keep.attributes = F)
      int_j[4] <- mean(abs(pred - pi_0))
      
      
      #Get integrated absolute error
      get_it <- function(z){
        pi <- predict(fit, newdata = data.frame(z = z, FU_years = 50))
        rsurv <- predict(fit, newdata = data.frame(z = z, FU_years = time))
        (rsurv - pi) / (1 - pi)
      }
      
      surv_u <- sapply(z_values, get_it)
      diae_j[4] <- upper.lim / 2 * sum(weights * rowMeans(abs(surv_u - survu_0)))
      
      #Get integrated absolute error - relative survival
      surv <- sapply(z_values, function(z) predict(fit, newdata = data.frame(z = z, FU_years = time)))
      diae_surv_j[4] <- upper.lim / 2 * sum(weights * rowMeans(abs(surv - surv_0)))
      
      
      #Proportional odds
      fit <- stpm2(Surv(FU_years, status) ~ z, data = sim_data$D, bhazard = sim_data$D$exp_haz, 
                   df = 5, cure = T, tvc = list(z = 1), link.type = "PO")
      
      #Get absolute error
      X <- rstpm2:::lpmatrix.lm(fit@lm, newdata = data.frame(z = c(1, 0), FU_years = max(sim_data$D$FU_years)))
      
      if(is.null(fit@vcov)){
        logor_j[5] <- get_lor2(fit, newdata = NULL, X = X)
      } else {
        est <- predictnl.default2(fit, fun = get_lor2, X = X, newdata = X) 
        logor_j[5] <- est$Estimate 
        var_j[5] <- est$SE
      }
      if(is.na(var_j[5])) break
      pred <- predict(fit, newdata = data.frame(z = z_values, FU_years = max(sim_data$D$FU_years)), 
                      keep.attributes = F)
      int_j[5] <- mean(abs(pred - pi_0))
      
      #Get integrated absolute error
      surv_u <- sapply(z_values, get_it)
      diae_j[5] <- upper.lim / 2 * sum(weights * rowMeans(abs(surv_u - survu_0)))
      
      #Get integrated absolute error - relative survival
      surv <- sapply(z_values, function(z) predict(fit, newdata = data.frame(z = z, FU_years = time)))
      diae_surv_j[5] <- upper.lim / 2 * sum(weights * rowMeans(abs(surv - surv_0)))
      
      #Save dataframe with results
      data.frame(logor = logor_j, var = var_j, int = int_j, 
                 diae = diae_j, diae_surv = diae_surv_j, 
                 model = 1:n.models)
    }, mc.cores = n.cores)
    
    L[[i]] <- do.call(rbind, res)
  }
  save(L, file = filename)
}


#Calculate biases

#Aggregate various results for
meanLOR <- lapply(1:length(L), function(i) aggregate(logor ~ model, data = L[[i]], FUN = mean)[,2])

meanCov <- lapply(1:length(L), function(i) aggregate(var ~ model, data = L[[i]], FUN = median)[,2])

meanNAs <- lapply(1:length(L), function(i) aggregate(var ~ model, data = L[[i]], FUN = function(x) mean(is.na(x)))[,2])

meanInt <- lapply(1:length(L), function(i) aggregate(int ~ model, data = L[[i]], FUN = mean)[,2])

meanDiae <- lapply(1:length(L), function(i) aggregate(diae ~ model, data = L[[i]], FUN = mean)[,2])

meanDiaeS <- lapply(1:length(L), function(i) aggregate(diae_surv ~ model, data = L[[i]], FUN = mean)[,2])

#Calculate percentage of NAs
for(j in 1:5){
  NAs <- lapply(1:length(L), function(i){
    wh <- which(is.na(L[[i]]$var))
    length(which(L[[i]][wh,]$model == j))
  }) 
  print(sum(unlist(NAs)) /( n.sim * length(L)) * 100)
}

#Do some rearrangements and output results to data frame
D <- data.frame(unlist(meanLOR), unlist(meanCov), unlist(meanNAs), 
                unlist(meanInt), unlist(meanDiae), unlist(meanDiaeS))
D <- cbind(Model = c("WMC", "GMC1", "GMC2", "LC1", "LC2"), D, stringsAsFactors = F)
names(D) <- c("Model", "LOR", "Var", "NA's","$A(\\hat\\pi, \\pi_0)$", "IA$(\\hat S_u, S_u)$", "IA$(\\hat R, R)$")
ncol.D <- ncol(D)
col_names <- "Model & LOR & SE(LOR) & NAs & AE$(\\hat\\pi(z), \\pi)$ & IAE$(\\hat S_u, S_u)$ & IAE$(\\hat R, R)$\\\\\n"
digits.D <- c(1, 2, 2, 2, 2, 2, 3, 3, 3)
na.all <- all(D$`NA's` == 0)

if(na.all){
  D <- subset(D, select = -c(`NA's`))
  ncol.D <- ncol(D)
  col_names <- gsub(" NAs &", "", col_names)
  digits.D <- c(1, 2, 2, 2, 3, 3, 3)
}

addtorow <- list()
addtorow$pos <- list(0, 0, 0, 5, 5, 10, 10, 15)
truelor <- sprintf("%.2f", trueor)
addtorow$command <- c(col_names,
                      "\\hline\n", 
                      paste0("\\multicolumn{", ncol.D, "}{l}{Scenario 1 (LOR = ", truelor[1],")}\\\\\n"),
                      "\\hline\n",
                      paste0("\\multicolumn{", ncol.D, "}{l}{Scenario 2 (LOR = ", truelor[2],")}\\\\\n"),
                      "\\hline\n",
                      paste0("\\multicolumn{", ncol.D, "}{l}{Scenario 3 (LOR = ", truelor[3],")}\\\\\n"), 
                      "\\hline\n")
print(xtable(D, align = paste0(rep("c", ncol.D + 1), collapse = ""), label = "tab:rescov", 
             caption = "Empirical mean LOR, AE, IAE, and median SE of the five models in Table \\ref{tab:modelscov}. 
             The IAE was calculated both for the survival of the uncured and the entire relative survival function. 
             LOR: log-odds ratio, SE: standard error, AE: absolute error, IAE: integrated absolute error.", 
             digits = digits.D),
      add.to.row = addtorow, include.colnames = F, include.rownames = F,
      sanitize.text.function = identity, hline.after = c(-1), 
      file = file.path(tab.out, "simulateResultsCov.tex"))#, size="\\fontsize{8pt}{10pt}\\selectfont")

