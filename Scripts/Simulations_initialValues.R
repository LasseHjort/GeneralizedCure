#Simulation options
n.sim <- 500
n.obs <- 1000
age <- 60
types <- c("weibull", "bs")

#Number of cores for the multi threading
n.cores <- 48

# #Set settings for survival of the uncured
# gaussxw <- statmod::gauss.quad(100)
# weights <- gaussxw$weights
# #Set upper limit of integral
# upper.lim <- 15
# time <- upper.lim / 2 * gaussxw$nodes + upper.lim / 2
# 
# 
# true_survus <- list(weibull = function(time, pars) exp(-pars[3] * time ^ pars[2]), 
#                     bs = function(time, pars) exp(-exp(cbind(1, bs(x = time, 
#                                                                    knots = c(1, 7), 
#                                                                    Boundary.knots = c(0, 15))) %*% pars[-1]))
# )

#Initialize list for results
L.all <- vector("list", length(n.obs))
for(i in 1:length(L.all)){
  L.all[[i]] <- vector("list", length(types))
  names(L.all[[i]]) <- types
}
names(L.all) <- n.obs

#Run simulations
for(type in types){
  cat(type, "\n")
  cases <- switch(type,
                  weibull = cases_wei,
                  gengamma = cases_gam,
                  bs = cases_gen)
  cases_obs <- lapply(cases, function(x) c(x, age = age, n = n.obs[k], type = type))
  filename <- file.path(data.out, paste0("sim_res_init_", type, "_", n.obs[k],".RData"))
  if(file.exists(filename)){
    load(filename)
  }else{
    set.seed(20190514, "L'Ecuyer")
    L <- vector("list", length(cases_obs))
    for(i in 1:length(cases_obs)){
      #covariance <- if(i != 6) T else F
      #var.type <- if(i != 6) "ci" else "n"
      cat(i, "\n")
      #survu_0 <- function(time) true_survus[[type]](time, pars = cases_obs[[i]][[1]])
      pi_0 <- cases_obs[[i]][[1]]["pi"]
      curerate <- mclapply(1:n.sim, function(j){
        #for(j in 1:n.sim){
        #cat(j, "\n")
        sim_data <- do.call(sim_surv, cases_obs[[i]])
        # 
        # fit <- rs.surv(Surv(FU, status) ~ ratetable(age = age, sex = sex, year = diag_date),
        #                 ratetable = survexp.dk, data = sim_data$D)
        # fit$time <- fit$time / 365.24
        # plot(fit)
        #f <- function(time) sim_data$rel_surv(time)
        #curve(expr = f, from = 1e-06, to = 15, add = T, col = 2)
        
        #Explicit cure point models
        #With conventional knot selection
        fit1 <- try(GenFlexCureModel(Surv(FU_years, status) ~ 1, 
                                     data = sim_data$D,
                                     bhazard = "exp_haz",
                                     df = 4, verbose = F, covariance = F, ini.types = "cure"))
        
        fit2 <- try(GenFlexCureModel(Surv(FU_years, status) ~ 1, 
                                     data = sim_data$D,
                                     bhazard = "exp_haz",
                                     df = 4, verbose = F, covariance = F, ini.types = "flexpara"))
        
        fit3 <- try(GenFlexCureModel(Surv(FU_years, status) ~ 1, 
                                     data = sim_data$D,
                                     bhazard = "exp_haz",
                                     df = 4, verbose = F, covariance = F, ini.types = "cure"))
        
        
        fits <- list(fit1, fit2, fit3)
        
        sapply(fits, function(fit){
          if(inherits(fit, "try-error")){
            NA
          } else{
            # if(is.null(fit1$covariance)){
            #   pi2 <- predict(fit1, type = "curerate", var.type = "n")[[1]]$Estimate
            #   cp2 <- NA
            # } else{ 
            #   pi2 <- predict(fit1, type = "curerate", var.type = var.type)[[1]]
            #   cp2 <- pi_0 <= pi2$upper & pi_0 >= pi2$lower
            #   pi2 <- pi2$Estimate
            # }
            # surv_u2 <- predict(fit1, time = time, type = "survuncured", var.type = "n")[[1]]$Estimate
            # iae2 <- upper.lim / 2 * sum(weights * abs(surv_u2 - survu_0(time))) 
            predict(fit, type = "curerate", var.type = "n")[[1]]$Estimate
          }
        })
      }, mc.cores = n.cores)
      
      L[[i]] <- do.call(rbind, curerate)
    }
    save(L, file = filename)
  }
  L.all[[as.character(n.obs[k])]][[type]] <- L
}


models <- c("Cure", "Flexpara", "Both")

#Reformat the data structure in both types

create_plot <- function(n){
  for(i in 1:length(L.all[[n]]$weibull)){
    crs <- L.all[[n]]$weibull[[i]]
    model <- rep(1:3, each = n.sim)
    L.all[[n]]$weibull[[i]] <- data.frame(crs = c(as.matrix(crs)), model = model, Scenario = i)
  }
  
  
  for(i in 1:length(L.all[[n]]$bs)){
    crs <- L.all[[n]]$bs[[i]]
    model <- rep(1:3, each = n.sim)
    L.all[[n]]$bs[[i]] <- data.frame(crs = c(as.matrix(crs)), model = model, Scenario = i)
  }
  
  #Combine data and plot the cure rate estimates
  res_wei <- do.call(rbind, L.all[[n]]$weibull)
  res_wei$type <- "Weibull"
  res_bs <- do.call(rbind,  L.all[[n]]$bs)
  res_bs$type <- "Polynomial spline"
  
  plot_data <- rbind(res_wei, res_bs)
  plot_data$type <- factor(plot_data$type, levels = c("Weibull", "Polynomial spline"))
  plot_data$Model <- factor(models[plot_data$model], levels = models)
  plot_data$Scenario <- factor(plot_data$Scenario)
  
  true_pis <- sapply(cases_wei, function(x) x[[1]][1])
  df <- data.frame(x1 = 1:6 - 0.5, y1 = true_pis, y2 = true_pis, x2 = 1:6 + 0.5, Scenario = 1, Model = "Both")
  
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



pdf(file.path(fig.out, "SimulateCureRateInit1000.pdf"), width = 13, height = 8)
create_plot("1000")
dev.off()

