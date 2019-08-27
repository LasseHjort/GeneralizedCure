###########Simulations for generalized cure paper##########


#Function for simulating survival data
sim_surv <- function(pars, age, n, type = "weibull"){
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
  
  #Set survival of the uncured, relative survival, density, and quantile functions 
  #for the Weibull and generalized gamma case
  if(type == "weibull"){
    surv_can_fun <- function(pars, time) exp(-pars[3] * time ^ pars[2]) 
    dens_can_fun <- function(pars, time) exp(-pars[3] * time ^ pars[2]) * pars[3] * pars[2] * time ^ (pars[2] - 1)
    rel_surv <- function(time) pars[1] + (1 - pars[1]) * exp(-pars[3] * time ^ pars[2])
    qrel_surv <- function(q){
      res <- rep(Inf, length(q))
      wh <- q < 1 - pars[1]
      res[wh] <- qweibull(q[wh] / (1 - pars[1]), 
                          shape = pars[2], 
                          scale = 1 / pars[3] ^ (1 / pars[2]))
      res
    }
  }
  if(type == "gengamma"){
    surv_can_fun <- function(pars, time) pgengamma.stacy(time, scale = pars[2], 
                                                         d = pars[3], k = pars[4], lower.tail = F)
    dens_can_fun <- function(pars, time) dgengamma.stacy(time, scale = pars[2], 
                                                         d = pars[3], k = pars[4])
    rel_surv <- function(time) pars[1] + (1 - pars[1]) * pgengamma.stacy(time, scale = pars[2], 
                                                                         d = pars[3], k = pars[4], lower.tail = F)
    qrel_surv <- function(q){
      res <- rep(Inf, length(q))
      wh <- q < 1 - pars[1]
      res[wh] <- qgengamma.stacy(q[wh] / (1 - pars[1]), 
                                 scale = pars[2], 
                                 d = pars[3], 
                                 k = pars[4])
      res
    }
  }
  
  #Survival of the uncured and relative survival for the polynomial spline case
  if(type == "bs"){
    inknots <- c(1, 7)
    bdrknots <- c(0, 15)
    surv_can_fun <- function(pars, time) exp(-exp(cbind(1, bs(time, knots = inknots, Boundary.knots = bdrknots)) %*% pars[-1]))
    rel_surv <- function(time) pars[1] + (1 - pars[1]) * surv_can_fun(pars = pars, time = time)
  }
  
  #Function for which roots have to found
  unifind <- function(time, u) surv_can_fun(pars = pars, time) - u
  
  #Function for simulating survival times by finding the root of unifind
  sim_fun_rel <- function(x){
    res <- rep(Inf, length(x))
    for(i in 1:length(x)){
      #cat(i, "\n")
      if(x[i] >= pars[1]){
        u.new <- (x[i] - pars[1]) /(1 - pars[1])
        eval <- unifind(1e-15, u.new)
        if(eval < 0){
          res[i] <- .Machine$double.eps
        } else {
          res[i] <- uniroot(unifind, interval = c(1e-15, 150), u = u.new)$root 
        }
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
  list(D = D, pars = pars, rel_surv = rel_surv, S_exp = S_exp, 
       surv_can_fun = surv_can_fun)#, dens_can_fun = dens_can_fun)
}


##Plot relative survival scenarios
#Weibull
cases_wei <- list(list(c(pi = 0.25, shape = 1, scale = 1)),
                  list(c(pi = 0.4, shape = 0.9, scale = 0.7)),
                  list(c(pi = 0.5, shape = 1.3, scale = 0.2)),
                  list(c(pi = 0.75, shape = 1, scale = 0.5)),
                  list(c(pi = 0.9, shape = 1.2, scale = 1.5)),
                  list(c(pi = 0.2, shape = 1.2, scale = 0.1)))

#Generalized gamma
cases_gam <- list(list(c(pi = 0.25, scale = 1.2, d =0.7, k = 0.7)),
                  list(c(pi = 0.25, scale = 5, d = 0.8, k = 0.7)),
                  list(c(pi = 0.5, scale = 3, d = 0.8, k = 0.7)),
                  list(c(pi = 0.75, scale = 1.2, d = 0.6, k = 0.7)),
                  list(c(pi = 0.75, scale = 6, d = 0.9, k = 0.7)),
                  list(c(pi = 0, scale = 5, d = 0.9, k = 1.2)))

#Polynomial spline
cases_gen <- list(list(c(pi = 0.25, -6, 4.5, 7.8, 7.8, 8, 8.5)),
                  list(c(pi = 0.4, -7.5, 6.2, 7.8, 8.5, 9.8, 10.9)),
                  list(c(pi = 0.5, -7.5, 5.3, 7.1, 8.9, 8.9, 9.2)),
                  list(c(pi = 0.75, -6, 4.5, 6.8, 7.8, 7.5, 8)),
                  list(c(pi = 0.9, -8.5, 6.2, 10, 11, 12, 13)),
                  list(c(pi = 0.2, -6.9, 4.8, 5.9, 7, 7.5, 7.7)))

#Check the trajectory of the polynomial splines
# f <- function(t){
#   pars <- cases_gen[[5]][[1]]
#   exp(-exp(cbind(1, bs(x = t, knots = c(1, 7), Boundary.knots = c(0, 15))) %*% pars[-1]))
# }
# 
# curve(f, from = 0, to = 30, ylim = c(0, 1), n = 1000)
# abline(v = 15, lty = 2, h = 0)


#Make table for supplementary containing the parameter values for the simulations
wei_pars <- t(sapply(cases_wei, function(x) x[[1]][-1]))
gen_pars <- t(sapply(cases_gen, function(x) x[[1]][-1]))

D <- as.data.frame(cbind(wei_pars, gen_pars))
D <- cbind(1:nrow(D), D)

addtorow <- list()
addtorow$pos <- list(0, 0)
addtorow$command <- c("& \\multicolumn{2}{|c|}{Weibull} & \\multicolumn{6}{c}{Polynomial splines}\\\\\n", 
                      "Scenario & $\\gamma_0$ & $\\gamma_1$ & $\\beta_0$ & $\\beta_1$ & $\\beta_2$ & $\\beta_3$ & $\\beta_4$ & $\\beta_5$\\\\\n")

print(xtable(D, caption = "Parameter values used for simulating survival data.", label = "tab:params", 
       align = "cc|cc|cccccc", digits = 1), include.rownames = F, include.colnames = F, add.to.row = addtorow, 
      santize.text.function = identity, file = file.path(tab.out, "ParamWithoutCov.tex"))


##Creat plot to show the relative survival trajectories
#Choose time points
time.points <- seq(0, 30, length.out = 100)

#Weibull
wei <- function(pars) pars[1] + (1 - pars[1]) * exp(-pars[3] * time.points ^ pars[2])
L <- lapply(cases_wei, function(pars) wei(pars[[1]]))
D_wei <- data.frame(surv = do.call(c, L), time.points = rep(time.points, length(cases_wei)), 
                    Scenario = rep(1:length(cases_wei), each = length(time.points)), 
                    Model = "Weibull")

#Check the relative survival trajectory
# gam <- function(pars) pars[1] + (1 - pars[1]) * pgengamma.stacy(time.points,
#                                                                 scale = pars[2],
#                                                                 d = pars[3],
#                                                                 k = pars[4],
#                                                                 lower.tail = F)
# 
# L <- lapply(cases_gam, function(pars) gam(pars[[1]]))
# D_gam <- data.frame(surv = do.call(c, L), time.points = rep(time.points, length(cases_gam)), 
#                     Scenario = rep(1:length(cases_gam), each = length(time.points)), 
#                     Model = "Generalized gamma")

#Generalized curves
gen <- function(pars) pars[1] + (1 - pars[1]) * exp(-exp(cbind(1, bs(x = time.points, knots = c(1, 7), Boundary.knots = c(0, 15))) %*% pars[-1]))
L <- lapply(cases_gen, function(pars) gen(pars[[1]]))
D_gen <- data.frame(surv = do.call(c, L), time.points = rep(time.points, length(cases_gen)), 
                    Scenario = rep(1:length(cases_gen), each = length(time.points)), 
                    Model = "Polynomial splines")

#Combine data and plot the trajectories
D <- rbind(D_wei, D_gen)
D$Scenario <- factor(D$Scenario)

p <- ggplot(D, aes(x = time.points, y = surv, colour = Scenario, linetype = Scenario)) + geom_line() +
  ylab("Relative survival") + xlab("Time") + theme_bw() + 
  scale_linetype_manual(values = c("solid", "solid", "dashed", 
                                   "dashed", "twodash", "twodash")) + 
  scale_color_manual(values = rep(c("black", "grey"), 3)) + 
  theme(legend.position = "bottom", 
        legend.text=element_text(size=13), 
        legend.title = element_text(size = 13),
        axis.title=element_text(size=14),
        strip.text = element_text(size=13), 
        axis.text = element_text(size = 13),
        legend.key.size = unit(2,"line")) +
  geom_vline(xintercept = 15, linetype = "dashed") + guides(colour = guide_legend(nrow = 1)) +
  facet_grid(.~ Model) + ylim(0,1)# + scale_color_brewer(palette = "Set2") 


pdf(file.path(fig.out, "Cases.pdf"), width = 8, height = 5)
print(p)
dev.off()


#Run the simulations without covariates
source("Scripts/SimpleSimulations.R")

#Run the simulations without covariates for varying degress of freedom
source("Scripts/DfSimulations.R")

#Run the simulations with covariates
source("Scripts/Simulations_cov.R")

#Run the simulations with different link functions
source("Scripts/Simulations_link.R")

#Run the simulations with different initial values
source("Scripts/Simulations_initialValues.R")

