#Simulation options
n.sim <- 500
n.obs <- 1000
age <- 60
types <- c("weibull", "bs")

#Number of cores for the multi threading
n.cores <- 48

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
    cases <- switch(type,
                    weibull = cases_wei,
                    gengamma = cases_gam,
                    bs = cases_gen)
    cases_obs <- lapply(cases, function(x) c(x, age = age, n = n.obs[k], type = type))
    filename <- file.path(data.out, paste0("sim_res_link_", type, "_", n.obs[k],".RData"))
    if(file.exists(filename)){
      load(filename)
    }else{
      set.seed(20190514, "L'Ecuyer")
      L <- vector("list", length(cases_obs))
      for(i in 1:length(cases_obs)){
        cat(i, "\n")
        pi_0 <- cases_obs[[i]][[1]]["pi"]
        curerate <- mclapply(1:n.sim, function(j){
          sim_data <- do.call(sim_surv, cases_obs[[i]])
          # Check relative survival trajectory
          # fit <- rs.surv(Surv(FU, status) ~ ratetable(age = age, sex = sex, year = diag_date),
          #                 ratetable = survexp.dk, data = sim_data$D)
          # fit$time <- fit$time / 365.24
          # plot(fit)
          #f <- function(time) sim_data$rel_surv(time)
          #curve(expr = f, from = 1e-06, to = 15, add = T, col = 2)
          
          #Mixture cure model with conventional knot selection
          fit1 <- try(GenFlexCureModel(Surv(FU_years, status) ~ 1, 
                                       data = sim_data$D,
                                       bhazard = "exp_haz",
                                       df = 4, verbose = F, covariance = F, 
                                       link.type = "PH", link.type.cr = "logit"))
          
          fit2 <- try(GenFlexCureModel(Surv(FU_years, status) ~ 1, 
                                       data = sim_data$D,
                                       bhazard = "exp_haz",
                                       df = 4, verbose = F, covariance = F, 
                                       link.type = "PH", link.type.cr = "loglog"))
          
          fit3 <- try(GenFlexCureModel(Surv(FU_years, status) ~ 1, 
                                       data = sim_data$D,
                                       bhazard = "exp_haz",
                                       df = 4, verbose = F, covariance = F, 
                                       link.type = "PO", link.type.cr = "logit"))
          
          fit4 <- try(GenFlexCureModel(Surv(FU_years, status) ~ 1, 
                                       data = sim_data$D,
                                       bhazard = "exp_haz",
                                       df = 4, verbose = F, covariance = F, 
                                       link.type = "PH", link.type.cr = "loglog"))
          
          fits <- list(fit1, fit2, fit3, fit4)
          
          sapply(fits, function(fit){
            if(inherits(fit, "try-error")){
              NA
            } else{
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
}


models <- c("PH-logit", "PH-loglog", "PO-logit", "PO-loglog")

#Reformat the data structure in both types

create_plot <- function(n){
  for(i in 1:length(L.all[[n]]$weibull)){
    crs <- L.all[[n]]$weibull[[i]]
    model <- rep(1:4, each = n.sim)
    L.all[[n]]$weibull[[i]] <- data.frame(crs = c(as.matrix(crs)), model = model, Scenario = i)
  }
  
  
  for(i in 1:length(L.all[[n]]$bs)){
    crs <- L.all[[n]]$bs[[i]]
    model <- rep(1:4, each = n.sim)
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
  df <- data.frame(x1 = 1:6 - 0.5, y1 = true_pis, y2 = true_pis, x2 = 1:6 + 0.5, Scenario = 1, Model = "PH-logit")
  
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


pdf(file.path(fig.out, "SimulateCureRateLink1000.pdf"), width = 13, height = 8)
create_plot("1000")
dev.off()

