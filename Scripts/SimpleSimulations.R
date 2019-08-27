#Simulation options
n.sim <- 500
n.obs <- c(1000, 500, 200)
age <- 60
types <- c("weibull", "bs")

#Number of cores for the multi threading
n.cores <- 48

#Set settings for survival of the uncured
gaussxw <- statmod::gauss.quad(100)
weights <- gaussxw$weights
#Set upper limit of integral
upper.lim <- 15
time <- upper.lim / 2 * gaussxw$nodes + upper.lim / 2

#List containing the true survival of the uncured
true_survus <- list(weibull = function(time, pars) exp(-pars[3] * time ^ pars[2]), 
                    bs = function(time, pars) exp(-exp(cbind(1, bs(x = time, 
                                                                   knots = c(1, 7), 
                                                                   Boundary.knots = c(0, 15))) %*% pars[-1]))
)

#Initialize list for results
L.all <- vector("list", length(n.obs))
for(i in 1:length(L.all)){
  L.all[[i]] <- vector("list", length(types))
  names(L.all[[i]]) <- types
}
names(L.all) <- n.obs

#Run simulations
for(k in 1:length(n.obs)){
  cat("k = ", k, "\n")
  for(type in types){
    cat(type, "\n")
    #Create list with parameters for simulation
    cases <- switch(type,
                    weibull = cases_wei,
                    gengamma = cases_gam,
                    bs = cases_gen)
    cases_obs <- lapply(cases, function(x) c(x, age = age, n = n.obs[k], type = type))
    
    #If simulations have already been run, the load these
    filename <- file.path(data.out, paste0("sim_res_", type, "_", n.obs[k],".RData"))
    if(file.exists(filename)){
      load(filename)
    }else{
      set.seed(20190514, "L'Ecuyer")
      L <- vector("list", length(cases_obs))
      for(i in 1:length(cases_obs)){
        #Confidence intervals are not computed for the last scenario
        covariance <- if(i != 6) T else F
        var.type <- if(i != 6) "ci" else "n"
        cat(i, "\n")
        #True survival of uncured and cure proportion
        survu_0 <- function(time) true_survus[[type]](time, pars = cases_obs[[i]][[1]])
        pi_0 <- cases_obs[[i]][[1]]["pi"]
        #Start simulations
        curerate <- mclapply(1:n.sim, function(j){
          #for(j in 1:n.sim){
          #cat(j, "\n")
          sim_data <- do.call(sim_surv, cases_obs[[i]])
          # Check imperical relative survival
          # fit <- rs.surv(Surv(FU, status) ~ ratetable(age = age, sex = sex, year = diag_date),
          #                 ratetable = survexp.dk, data = sim_data$D)
          # fit$time <- fit$time / 365.24
          # plot(fit)
          #f <- function(time) sim_data$rel_surv(time)
          #curve(expr = f, from = 1e-06, to = 15, add = T, col = 2)
          
          #True model - Weibull
          fit <- try(fit.cure.model(Surv(FU_years, status) ~ 1, data = sim_data$D, 
                                    bhazard = "exp_haz", covariance = covariance))
          if(inherits(fit, "try-error")){
            pi2 <- cp2 <- iae2 <- NA
          } else{
            if(is.null(fit$covariance)){
              pi1 <- predict(fit, type = "curerate", var.type = "n")[[1]]$Estimate
              cp1 <- NA
            } else{ 
              pi1 <- predict(fit, type = "curerate", var.type = var.type)[[1]]
              cp1 <- pi_0 <= pi1$upper & pi_0 >= pi1$lower
              pi1 <- pi1$Estimate
            }
            surv_u1 <- predict(fit, time = time, type = "survuncured", var.type = "n")[[1]]$Estimate
            iae1 <- upper.lim / 2 * sum(weights * abs(surv_u1 - survu_0(time)))
          }
          
          #Explicit cure point models
          #With conventional knot selection
          fit1 <- try(GenFlexCureModel(Surv(FU_years, status) ~ 1, 
                                       data = sim_data$D,
                                       bhazard = "exp_haz",
                                       df = 4, verbose = F, covariance = covariance))
          if(inherits(fit1, "try-error")){
            pi2 <- cp2 <- iae2 <- NA
          } else{
            if(is.null(fit1$covariance)){
              pi2 <- predict(fit1, type = "curerate", var.type = "n")[[1]]$Estimate
              cp2 <- NA
            } else{ 
              pi2 <- predict(fit1, type = "curerate", var.type = var.type)[[1]]
              cp2 <- pi_0 <= pi2$upper & pi_0 >= pi2$lower
              pi2 <- pi2$Estimate
            }
            surv_u2 <- predict(fit1, time = time, type = "survuncured", var.type = "n")[[1]]$Estimate
            iae2 <- upper.lim / 2 * sum(weights * abs(surv_u2 - survu_0(time))) 
          }
          
          #With knots only in the beginning of the follow-up
          event_times <- sim_data$D$FU_years[sim_data$D$status == 1]
          knots <- log(c(min(event_times), 0.5, 1, 3, 5))
          bdrknots <- knots[c(1, length(knots))]
          knots <- knots[-c(1, length(knots))]
          fit2 <- try(GenFlexCureModel(Surv(FU_years, status) ~ 1, data = sim_data$D, 
                                       smooth.formula = ~ nsx(log(FU_years), knots = knots, Boundary.knots = bdrknots),
                                       bhazard = "exp_haz", verbose = F, covariance = covariance))      
          if(inherits(fit2, "try-error")){
            pi3 <- cp3 <- iae3 <- NA
          }else{
            if(is.null(fit2$covariance)){
              pi3 <- predict(fit2, type = "curerate", var.type = "n")[[1]]$Estimate
              cp3 <- NA
            } else{ 
              pi3 <- predict(fit2, type = "curerate", var.type = var.type)[[1]]
              cp3 <- pi_0 <= pi3$upper & pi_0 >= pi3$lower
              pi3 <- pi3$Estimate
            }
            surv_u3 <- predict(fit2, time = time, type = "survuncured", var.type = "n")[[1]]$Estimate
            iae3 <- upper.lim / 2 * sum(weights * abs(surv_u3 - survu_0(time))) 
          }
          
          #Latent cure point models
          #With the last knot within the follow-up
          knots <- log(quantile(event_times, probs = c(0, 0.25, 0.5, 0.75, 0.95, 1)))
          fit <- stpm2(Surv(FU_years, status) ~ 1, data = sim_data$D, bhazard = sim_data$D$exp_haz, 
                       smooth.formula = ~ nsx(x = log(FU_years), knots = knots[-c(1, length(knots))], 
                                              Boundary.knots = range(knots), cure = T))
          
          pi4 <- try(predict(fit, newdata = data.frame(age = 50, FU_years = max(exp(knots))), se.fit = T))
          if(inherits(pi4, "try-error")){
            pi4 <- iae4 <- cp4 <- NA
          }else{
            surv4 <- predict(fit, newdata = data.frame(FU_years = time))
            surv_u4 <- (surv4 - pi4$Estimate) / (1 - pi4$Estimate)
            iae4 <- upper.lim / 2 * sum(weights * abs(surv_u4 - survu_0(time)))
            cp4 <- pi_0 <= pi4$upper & pi_0 >= pi4$lower
            pi4 <- pi4$Estimate 
          }
          
          #With the last knot beyond the available follow-up
          knots <- log(c(quantile(event_times, probs = c(0, 0.25, 0.5, 0.75)), 8, 20))
          fit <- stpm2(Surv(FU_years, status) ~ 1, data = sim_data$D, bhazard = sim_data$D$exp_haz, 
                       smooth.formula = ~ nsx(x = log(FU_years), knots = knots[-c(1, length(knots))], 
                                              Boundary.knots = range(knots), cure = T))
          pi5 <- try(predict(fit, newdata = data.frame(age = 50, FU_years = max(exp(knots))), se.fit = T))
          if(inherits(pi5, "try-error")){
            pi5 <- iae5 <- cp5 <- NA
          }else{
            surv5 <- predict(fit, newdata = data.frame(FU_years = time))
            surv_u5 <- (surv5 - pi5$Estimate) / (1 - pi5$Estimate)
            iae5 <- upper.lim / 2 * sum(weights * abs(surv_u5 - survu_0(time)))
            cp5 <- pi_0 <= pi5$upper & pi_0 >= pi5$lower
            pi5 <- pi5$Estimate
          }
          #}
          #Ouput list, CP is not reported in scenario 6
          if(i == 6){
            list(CR = data.frame(pi1, pi2, pi3, pi4, pi5),
                 SU = data.frame(iae1, iae2, iae3, iae4, iae5))
          } else{
            list(CR = data.frame(pi1, pi2, pi3, pi4, pi5),
                 SU = data.frame(iae1, iae2, iae3, iae4, iae5),
                 CP = data.frame(cp1, cp2, cp3, cp4, cp5)) 
          }
          
        }, mc.cores = n.cores)
        
        L1 <- do.call(rbind, lapply(curerate, function(x)return(x$CR)))
        L2 <- do.call(rbind, lapply(curerate, function(x)return(x$SU)))
        L3 <- do.call(rbind, lapply(curerate, function(x)return(x$CP)))
        
        L[[i]] <- list(CR = L1, SU = L2, CP = L3)
        #L_wei[[i]]$M <- M
      }
      save(L, file = filename)
    }
    L.all[[as.character(n.obs[k])]][[type]] <- L
  }
}

#List model names
models <- c("WMC", "GMC1", "GMC2", "LC1", "LC2")


#Function for creating boxplots of the cure proportions
create_plot <- function(n){
  #Reformat weibull results
  for(i in 1:length(L.all[[n]]$weibull)){
    crs <- L.all[[n]]$weibull[[i]]$CR
    model <- rep(1:5, each = n.sim)
    L.all[[n]]$weibull[[i]]$CR <- data.frame(crs = c(as.matrix(crs)), model = model, Scenario = i)
  }
  
  #Reformat splines results
  for(i in 1:length(L.all[[n]]$bs)){
    crs <- L.all[[n]]$bs[[i]]$CR
    model <- rep(1:5, each = n.sim)
    L.all[[n]]$bs[[i]]$CR <- data.frame(crs = c(as.matrix(crs)), model = model, Scenario = i)
  }
  
  #Combine data and plot the cure rate estimates
  res_wei <- do.call(rbind, lapply(L.all[[n]]$weibull, function(x) return(x$CR)))
  res_wei$type <- "Weibull"
  res_bs <- do.call(rbind,  lapply(L.all[[n]]$bs, function(x) return(x$CR)))
  res_bs$type <- "Polynomial spline"
  
  plot_data <- rbind(res_wei, res_bs)
  plot_data$type <- factor(plot_data$type, levels = c("Weibull", "Polynomial spline"))
  plot_data$Model <- factor(models[plot_data$model], levels = models)
  plot_data$Scenario <- factor(plot_data$Scenario)
  
  true_pis <- sapply(cases_wei, function(x) x[[1]][1])
  df <- data.frame(x1 = 1:6 - 0.5, y1 = true_pis, y2 = true_pis, x2 = 1:6 + 0.5, Scenario = 1, Model = "WMC")
  
  p <- ggplot(plot_data, aes(x = Scenario, y = crs, group = Scenario:Model, fill = Model)) +
    geom_boxplot(outlier.size = 0.8) + ylab("Estimated cure probability") + xlab("Scenario") + theme_bw() +
    theme(legend.position = "bottom", 
          legend.text=element_text(size=18), 
          legend.title = element_text(size = 18),
          axis.title=element_text(size=20),
          strip.text = element_text(size=18), 
          axis.text = element_text(size = 17)) + 
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) + facet_grid(~ type) + 
    scale_fill_brewer(palette = "Set2") + 
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df, linetype = "dashed", color = "black")
  print(p)
}


#Create figures for n = 1000, 500, and 200
pdf(file.path(fig.out, "SimulateCureRate1000.pdf"), width = 13, height = 8)
create_plot("1000")
dev.off()

pdf(file.path(fig.out, "SimulateCureRate500.pdf"), width = 13, height = 8)
create_plot("500")
dev.off()

pdf(file.path(fig.out, "SimulateCureRate200.pdf"), width = 13, height = 8)
create_plot("200")
dev.off()

#Create table with survival of uncured and CP
create_table <- function(n, type){
  for(i in 1:length(L.all[[n]][[type]])){
    sus <- L.all[[n]][[type]][[i]]$SU
    model <- rep(1:5, each = n.sim)
    L.all[[n]][[type]][[i]]$SU <- data.frame(sus = c(as.matrix(sus)), model = model, Scenario = i)
  }
  
  L <- lapply(L.all[[n]][[type]], function(x) aggregate(sus ~ model, data = x$SU, FUN = mean))
  D1 <- as.data.frame(sapply(L, function(x) sprintf("%.3f", x[,2])), stringsAsFactors = F)
  #D1 <- round(D1, 3)
  
  
  D2 <- sapply(L.all[[n]][[type]], function(x){
    if(is.null(x$CP)){
      rep(NA, length(models))
    }else{
      df <- x$CP
      cp <- sprintf("%.1f", colMeans(df, na.rm = T) * 100)
      cp.missing <- sprintf("%.1f", colMeans(is.na(df)) * 100)
      paste0(cp, "(", cp.missing, ")")
    }
  })
  
  rbind(D1, D2)
}


D1 <- create_table("1000", type = "weibull")
D2 <- create_table("1000", type = "bs")
D <- rbind(D1, D2)

#Do some rearrangements and output results to data frame
names(D) <- paste0("Scenario ", 1:6)
D <- cbind("Models" = rep(models, 4), D, stringsAsFactors = F)
ncol.D <- ncol(D)

col_names <- paste0(" &", paste0("Scenario ", 1:6, collapse = "&"), "\\\\\n")
digits.D <- c(1, 2, 2, 2, 2, 2, 3, 3, 3)
na.all <- all(D$`NA's` == 0)

addtorow <- list()
addtorow$pos <- list(0, 0, 0, 5, 5, 10, 10, 10, 15, 15, 20)
addtorow$command <- c(col_names, 
                      "\\hline\n", 
                      paste0("\\multicolumn{", ncol.D, "}{l}{IAE$(\\hat S_u, S_u)$ - Weibull}\\\\\n"),
                      "\\hline\n",
                      paste0("\\multicolumn{", ncol.D, "}{l}{CP - Weibull}\\\\\n"),
                      "\\hline\n", 
                      "\\hline\n",
                      paste0("\\multicolumn{", ncol.D, "}{l}{IAE$(\\hat S_u, S_u)$ - polynomial spline}\\\\\n"),
                      "\\hline\n",
                      paste0("\\multicolumn{", ncol.D, "}{l}{CP - polynomial spline}\\\\\n"),
                      "\\hline\n")

print(xtable(D, align = paste0(rep("l", ncol.D + 1), collapse = ""), label = "tab:ressimple", 
             caption = "Empirical mean IAE of $S_u(t)$ and CP for $\\pi$. 
             For the CPs, the proportion of simulations where a confidence interval could not be computed is reported in parentheses.
             As appropriate estimation of the cure proportion is not expected in scenario 6, the CPs were not computed in this scenario.
             IAE: integrated absolute error, CP: coverage probability."),
      add.to.row = addtorow, include.colnames = F, include.rownames = F,
      sanitize.text.function = identity, hline.after = c(-1), 
      file = file.path(tab.out, "simulateResultsSimple.tex"))#, size="\\fontsize{8pt}{10pt}\\selectfont")

