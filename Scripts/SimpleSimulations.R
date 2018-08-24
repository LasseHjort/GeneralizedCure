#Simulation options
n.sim <- 500
n.obs <- 1000
age <- 60
types <- c("weibull", "bs")

#Number of cores for the multi threading
n.cores <- 40


#Initialize list for results
L.all <- vector("list", length(types))
names(L.all) <- types

#Run simulations
for(type in types){
  cat(type, "\n")
  cases <- switch(type,
                  weibull = cases_wei,
                  gengamma = cases_gam,
                  bs = cases_gen)
  cases_obs <- lapply(cases, function(x) c(x, age = age, n = n.obs, type = type))
  filename <- file.path(data.out, paste0("sim_res_", type, ".RData"))
  if(file.exists(filename)){
    load(filename)
  }else{
    set.seed(20180625, "L'Ecuyer")
    L <- vector("list", length(cases_obs))
    for(i in 1:length(cases_obs)){
      cat(i, "\n")
      curerate <- mclapply(1:n.sim, function(j){
        #cat(j, "\n")
        sim_data <- do.call(sim_surv, cases_obs[[i]])
        # 
        #fit <- rs.surv(Surv(FU, status) ~ ratetable(age = age, sex = sex, year = diag_date),
        #                ratetable = survexp.dk, data = sim_data$D)
        #fit$time <- fit$time / 365.24
        #plot(fit)
        #f <- function(time) sim_data$rel_surv(time)
        #curve(expr = f, from = 1e-06, to = 15, add = T, col = 2)
        
        #True model - Weibull
        fit <- fit.cure.model(Surv(FU_years, status) ~ 1, data = sim_data$D, 
                              bhazard = "exp_haz", covariance = F)
        pi1 <- predict(fit, type = "curerate")[[1]]$Estimate
        
        
        #Explicit cure point models
        #With conventional knot selection
        fit1 <- GenFlexCureModel(Surv(FU_years, status) ~ 1, 
                                 data = sim_data$D,
                                 bhazard = "exp_haz",
                                 df = 4, verbose = F, covariance = F)
        pi2 <- predict(fit1, type = "curerate", var.type = "n")[[1]]$Estimate
        
        #With knots only in the beginning of the follow-up
        event_times <- sim_data$D$FU_years[sim_data$D$status == 1]
        knots <- log(c(min(event_times), 0.5, 1, 3, 5))
        bdrknots <- knots[c(1, length(knots))]
        knots <- knots[-c(1, length(knots))]
        fit2 <- GenFlexCureModel(Surv(FU_years, status) ~ 1, data = sim_data$D, 
                                 smooth.formula = ~ nsx(log(FU_years), knots = knots, Boundary.knots = bdrknots),
                                 bhazard = "exp_haz", verbose = F, covariance = F)      
        
        pi3 <- predict(fit2, type = "curerate", var.type = "n")[[1]]$Estimate
        
        #Latent cure point models
        #With the last knot within the follow-up
        knots <- log(quantile(event_times, probs = c(0, 0.25, 0.5, 0.75, 0.95, 1)))
        fit <- stpm2(Surv(FU_years, status) ~ 1, data = sim_data$D, bhazard = sim_data$D$exp_haz, 
                     smooth.formula = ~ nsx(x = log(FU_years), knots = knots[-c(1, length(knots))], 
                                            Boundary.knots = range(knots), cure = T))
        pi4 <- predict(fit, newdata = data.frame(age = 50, FU_years = max(exp(knots))))
        
        #With the last knot beyond the available follow-up
        knots <- log(c(quantile(event_times, probs = c(0, 0.25, 0.5, 0.75)), 8, 20))
        fit <- stpm2(Surv(FU_years, status) ~ 1, data = sim_data$D, bhazard = sim_data$D$exp_haz, 
                     smooth.formula = ~ nsx(x = log(FU_years), knots = knots[-c(1, length(knots))], 
                                            Boundary.knots = range(knots), cure = T))
        pi5 <- predict(fit, newdata = data.frame(age = 50, FU_years = max(exp(knots))))
        
        data.frame(pi1, pi2, pi3, pi4, pi5)
        
      }, mc.cores = n.cores)
      
      
      L[[i]] <- do.call(rbind, curerate)
      #L_wei[[i]]$M <- M
    }
    save(L, file = filename)
  }
  L.all[[type]] <- L
}

#Reformat the data structure in both types
for(i in 1:length(L.all$weibull)){
  crs <- L.all$weibull[[i]]
  model <- rep(1:5, each = n.sim)
  L.all$weibull[[i]] <- data.frame(crs = c(as.matrix(crs)), model = model, Scenario = i)
}


for(i in 1:length(L.all$bs)){
  crs <- L.all$bs[[i]]
  model <- rep(1:5, each = n.sim)
  L.all$bs[[i]] <- data.frame(crs = c(as.matrix(crs)), model = model, Scenario = i)
}

#Combine data and plot the cure rate estimates
res_wei <- do.call(rbind, L.all$weibull)
res_wei$type <- "Weibull"
res_bs <- do.call(rbind, L.all$bs)
res_bs$type <- "Polynomial spline"

plot_data <- rbind(res_wei, res_bs)
plot_data$type <- factor(plot_data$type, levels = c("Weibull", "Polynomial spline"))
plot_data$Model <- factor(LETTERS[plot_data$model])
plot_data$Scenario <- factor(plot_data$Scenario)


p <- ggplot(plot_data, aes(x = Scenario, y = crs, group = Scenario:Model, fill = Model)) +
  geom_boxplot(outlier.size = 0.8) + ylab("Estimated cure rates") + xlab("Scenario") + theme_bw() +
  theme(legend.position = "bottom", 
        legend.text=element_text(size=18), 
        legend.title = element_text(size = 18),
        axis.title=element_text(size=20),
        strip.text = element_text(size=18), 
        axis.text = element_text(size = 17)) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) + facet_grid(~ type) + 
  scale_fill_brewer(palette = "Set2") 


pdf(file.path(fig.out, "SimulateCureRate.pdf"), width = 13, height = 8)
print(p)
dev.off()
