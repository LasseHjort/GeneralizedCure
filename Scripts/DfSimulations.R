###Simulations for degrees of freedom

#Simulation options
n.sim <- 500
n.obs <- 1000
age <- 60
degrees <- c(2, 4, 6, 8, 10, 12)
types <- c("weibull", "bs")

#Number of cores for the multi threading
n.cores <- 48

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
  filename <- file.path(data.out, paste0("sim_res_degree_", type, ".RData"))
  if(file.exists(filename)){
    load(filename)
  }else{
    set.seed(20190514, "L'Ecuyer")
    L <- vector("list", length(cases_obs))
    for(i in 1:length(cases_obs)){
      cat(i, "\n")
      curerate <- mclapply(1:n.sim, function(j){
        #cat(j, "\n")
        sim_data <- do.call(sim_surv, cases_obs[[i]])
        
        #Looping over the degrees of freedom in degrees
        res <- matrix(nrow = 1, ncol = length(degrees))
        for(k in 1:length(degrees)){
          fit <- try(GenFlexCureModel(Surv(FU_years, status) ~ 1, 
                                  data = sim_data$D,
                                  bhazard = "exp_haz",
                                  df = degrees[k] - 1 , verbose = F, covariance = F))
          if(inherits(fit, "try-error")){
            res[1, k] <- NA
          }else{
            res[1, k] <- predict(fit, type = "curerate", var.type = "n")[[1]]$Estimate 
          }
        }
        as.data.frame(res)
        
      }, mc.cores = n.cores)
      
      
      L[[i]] <- do.call(rbind, curerate)
    }
    save(L, file = filename)
  }
  L.all[[type]] <- L
}

#Reformat the data structure in both types
for(i in 1:length(L.all$weibull)){
  crs <- L.all$weibull[[i]]
  model <- rep(1:length(degrees), each = n.sim)
  L.all$weibull[[i]] <- data.frame(crs = c(as.matrix(crs)), model = model, Scenario = i)
}


for(i in 1:length(L.all$bs)){
  crs <- L.all$bs[[i]]
  model <- rep(1:length(degrees), each = n.sim)
  L.all$bs[[i]] <- data.frame(crs = c(as.matrix(crs)), model = model, Scenario = i)
}

#Combine data and plot the cure rate estimates
res_wei <- do.call(rbind, L.all$weibull)
res_wei$type <- "Weibull"
res_bs <- do.call(rbind, L.all$bs)
res_bs$type <- "Polynomial spline"


plot_data <- rbind(res_wei, res_bs)
plot_data$type <- factor(plot_data$type, levels = c("Weibull", "Polynomial spline"))
plot_data$Scenario <- factor(plot_data$Scenario)
plot_data$model <- factor(plot_data$model, labels = degrees)
plot_data$Scenario <- factor(plot_data$Scenario)

true_pis <- sapply(cases_wei, function(x) x[[1]][1])
df <- data.frame(x1 = 1:6 - 0.5, y1 = true_pis, y2 = true_pis, x2 = 1:6 + 0.5, Scenario = 1, model = "2")

p <- ggplot(plot_data, aes(x = Scenario, y = crs, group = Scenario:model, fill = model)) +
  geom_boxplot(outlier.size = 0.8) + ylab("Estimated cure probability") + xlab("Scenario") + theme_bw() +
  scale_fill_brewer(name = "Number of parameters", palette = "Set2") + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) + facet_grid(~ type) + 
  theme(legend.position = "bottom", 
        legend.text=element_text(size=18), 
        legend.title = element_text(size = 18),
        axis.title=element_text(size=20),
        strip.text = element_text(size=18), 
        axis.text = element_text(size = 17)) + 
  guides(fill = guide_legend(nrow = 1)) + 
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df, linetype = "dashed", color = "black")

pdf(file.path(fig.out, "SimulateCureRateDegrees.pdf"), width = 13, height = 8)
print(p)
dev.off()

