#Create age group variable
Colon$age_groups <- cut(Colon$age_years, breaks = c(0, 55, 65, 75, 120), labels = c("-55", "55-65", "65-75", "75-"))

#Compute non-parametric and parametric relative survival estimates
levs <- levels(Colon$age_groups)
L <- L_conf <- vector("list", length(levs))
for(i in 1:length(levs)){
  these <- Colon$age_groups == levs[i]
  #Non parametric
  fit_nonp <- rs.surv(Surv(FU, status) ~ ratetable(age = age, sex = sex, year = diag_date), 
                      data = Colon[these,],
                      ratetable = survexp.dk, 
                      method = "ederer1") 
  D1 <- data.frame(surv = fit_nonp$surv, time = fit_nonp$time / ayear, method = "Ederer I")
  L_conf[[i]] <- data.frame(lower = fit_nonp$lower, upper = fit_nonp$upper, 
                            time = fit_nonp$time / ayear, lev = levs[i], method = "Ederer I")
  
  #Stpm2 model
  fitp <- stpm2(Surv(FU_years, status) ~ 1, data = Colon[these, ], 
                bhazard = Colon$exp_haz[these], df = 6)
  
  # knots <- quantile(Colon$FU_years[Colon$status == 1], probs = seq(0, 1, length.out = 6))
  # knots <- log(sort(c(knots, 10)))
  # bdr.knots <- range(knots)
  # inner.knots <- knots[-c(1, length(knots))]
  # fitp <- stpm2(Surv(FU_years, status) ~ 1, data = Colon[these, ], 
  #               bhazard = Colon$exp_haz[these], smooth.formula = ~nsx(log(FU_years), knots = inner.knots, Boundary.knots = bdr.knots))
  
  pred <- predict(fitp, newdata = data.frame(FU_years = D1$time), type = "surv")
  D2 <- data.frame(surv = pred, time = D1$time, method = "stpm2")
  
  #Stpm2 - cure model
  fitpc <- stpm2(Surv(FU_years, status) ~ 1, data = Colon[these, ], 
                 bhazard = Colon$exp_haz[these], df = 6, cure = T)
  
  predc <- predict(fitpc, newdata = data.frame(FU_years = D1$time), type = "surv")
  D3 <- data.frame(surv = predc, time = D1$time, method = "stpm2-cure")
  
  #Flexible mixture cure model
  fitmix <- GenFlexCureModel(Surv(FU_years, status) ~ 1, data = Colon[these,], 
                             bhazard = Colon$exp_haz[these], df = 5, verbose = F)
  
  predmix <- predict(fitmix, time = D1$time)
  D4 <- data.frame(surv = predmix[[1]]$Estimate, time = D1$time, method = "Mixture cure")
  D <- rbind(D1, D2, D3, D4)
  D$lev <- levs[i]
  L[[i]] <- D
}

#Plot relative survival estimates for each age group and output to file
plot_data <- do.call(rbind, L)
plot_data_conf <- do.call(rbind, L_conf)
pal <- c("black", brewer.pal(3, "Set2"))

p <- ggplot(plot_data, aes(x = time, y = surv, group = method, colour = method)) + 
  geom_line(size = 0.6) + facet_wrap(~lev, ncol = 2) + ylab("Relative survival") + 
  xlab("Time since diagnosis (years)") +  coord_cartesian(ylim = c(0, 1)) + 
  geom_step(data = plot_data_conf, aes(x = time, y = lower), linetype = "dashed", size = 0.4) + 
  geom_step(data = plot_data_conf, aes(x = time, y = upper), linetype = "dashed", size = 0.4) +
  scale_colour_manual(values = c("Ederer I" = pal[1], 
                                 "stpm2" = pal[2], 
                                 "stpm2-cure" = pal[3],
                                 "Mixture cure" = pal[4])) + 
  theme_bw() + 
  theme(legend.title = element_blank(), legend.position = "bottom", 
        legend.text=element_text(size=15), 
        axis.title=element_text(size=17),
        strip.text = element_text(size=15), 
        axis.text = element_text(size = 13),
        legend.key.size = unit(2,"line"))


pdf(file.path(fig.out, "ColonRelSurv.pdf"), width = 10, height = 7.3)
print(p)
dev.off()


##Covariate assessment
#Without time-varying effects
#Select variables
vars <- c("Age", "Gender", "Charlson", "Metastases", "Stage")

#Fit the Esteve model from package "relsurv"
#Select knots
deaths <- Colon$FU_years[Colon$status == 1]
knots <- c(0, quantile(deaths, probs = c(0.2, 0.4, 0.6, 0.8, 1)))
knots <- quantile(deaths, probs = seq(0, 1, length.out = 9))
knots[which.min(knots)] <- 0
#Fit model
fitnp <- rsadd(Surv(FU, status) ~ Age + Gender + Charlson + Metastases + Stage + 
                 ratetable(age = age, sex = sex, year = diag_date), 
               data = Colon, ratetable = survexp.dk, int = knots)
#Tried to fit piecewise constant excess hazard model. Does not work due to the lack of differentiability
# x <- Colon$FU_years
# get_b <- function(knots, x){
#   indices <- findInterval(x, knots)
#   basis <- matrix(0, ncol = length(knots), nrow = length(x))
#   for(i in 1:nrow(basis)){
#     basis[i, indices[i]] <- 1
#   }
#   basis <- basis[,-ncol(basis)]
#   attributes(basis) <- c(attributes(basis), list(knots = knots))
#   basis
# }
# a <- get_b(knots = knots, x = Colon$FU_years)
# 
# 
# fit <- stpm2(Surv(FU_years, status) ~ -1 + Age + Gender + Charlson + Metastases + Stage, data = Colon, 
#             bhazard = Colon$exp_haz, smooth.formula = ~get_b(x = FU_years, knots = knots))

#Fit stpm2 and stpm2-cure models
fitp <- stpm2(Surv(FU_years, status) ~ Age + Gender + Charlson + Metastases + Stage, data = Colon, 
             bhazard = Colon$exp_haz, df = 6)

fitpc <- stpm2(Surv(FU_years, status) ~ Age + Gender + Charlson + Metastases + Stage, data = Colon, 
             bhazard = Colon$exp_haz, df = 6, cure = T)

# plot(fit, newdata = data.frame(Age = 0, Gender = 0, Charlson = 0, Metastases = 0, Stage = 0), ylim = c(0, 1))
# 
# a <- rs.surv.rsadd(fitnp, newdata = data.frame(Age = 0, Gender = 0, Charlson = 0, 
#                                           Metastases = 0, Stage = 0, sex = "male", diag_date = 2010, age = 50))
#plot(a)

#Function for computing estimates and corresponding confidence intervals
get_sum <- function(fit){
  sum <- as.data.frame(summary(fit)@coef[vars,])
  sum$lower <- sum$Estimate - qnorm(0.975) * sum$`Std. Error`
  sum$upper <- sum$Estimate + qnorm(0.975) * sum$`Std. Error`
  sum <- round(sum, 3)
  est <- paste0(sum$Estimate, "(", sum$lower, ";", sum$upper, ")") 
  data.frame(est, pval = sprintf(tab.format, sum[, "Pr(z)"]), stringsAsFactors = F)
}

#Combine results into dataframe and output to file
D <- cbind(get_sum(fitp), get_sum(fitpc))
coefs <- summary(fitnp)$coefficients
coefs <- coefs[!grepl("fu ", rownames(coefs)),]
lower <- sprintf(tab.format, coefs[,1] - coefs[,2] * qnorm(0.975))
upper <- sprintf(tab.format, coefs[,1] + coefs[,2] * qnorm(0.975))
est <- paste0(sprintf(tab.format, coefs[,1]), "(", lower, ";", upper, ")")
df <- data.frame(est = est, pval = sprintf(tab.format, coefs[,4]), stringsAsFactors = F)

D <- cbind(df, D)
names(D) <- NULL

D <- rbind(rep(c("$\\beta$", "p-value"), 3), D)
rownames(D) <- c("", gsub("_", "\\\\_", vars))
rownames(D) <- c("","Age", "Male gender", "CCI $\\geq$2", "Metastases", "UICC Stage III-IV")



addtorow <- list()
addtorow$pos <- list(0)
addtorow$command <- c("& \\multicolumn{2}{c|}{Est\\`{e}ve ASM} & \\multicolumn{2}{c|}{stpm2} & \\multicolumn{2}{c}{stpm2-cure}\\\\\n")
print(xtable(D, align = c("lcc|cc|cc"), label = "tab:covcolon", 
             caption = "Effect estimates and their 95\\% confidence intervals estimated 
             from proportional excess hazards models. 
             ASM: additive survival model, CCI: Charlson Comorbidity Index, 
             UICC: Union for International Cancer Control."), 
      add.to.row = addtorow, include.colnames = FALSE, 
      sanitize.text.function = identity,
      file = file.path(tab.out, "coefscomp.tex"))

#Save coefficients for usage in latter plot
coefs <- fitp@coef[2:6]
coefs.cure <- fitpc@coef[2:6]


##With time-varying coefficients
#Fit stpm2
fit <- stpm2(Surv(FU_years, status) ~ Age + Gender + Charlson + Metastases + Stage, data = Colon, 
             tvc = list(Age = 3, Gender = 3, Charlson = 3, Metastases = 3, Stage = 3),
             bhazard = Colon$exp_haz, df = 6)
#Fit stpm2 - cure
fitc <- stpm2(Surv(FU_years, status) ~ Age + Gender + Charlson + Metastases + Stage, data = Colon, 
              tvc = list(Age = 3, Gender = 3, Charlson = 3, Metastases = 3, Stage = 3),
              bhazard = Colon$exp_haz, df = 6, cure = T)

for(age in c(50, 60, 70, 80)){
  plot(fitc, newdata = data.frame(FU_years = seq(0.001, 17), Charlson = 1, Stage = 1, Metastases = 0, 
                                  Gender = 0, Age = age), ylim = c(0, 1), main = age) 
}


#Check for the proportional hazards assumption of each variable by doing a likelihood ratio test.
#In each loop, one time-varying effect is removed to assess the assumption
res1 <- res2 <- rep(NA, length(vars))
for(i in 1:length(vars)){
  cat(vars[i], "\n")
  tvc.list <- list(3, 3, 3, 3)
  names(tvc.list) <- vars[-i]
  fit2 <- stpm2(Surv(FU_years, status) ~ Age + Gender + Charlson + Metastases + Stage, data = Colon, 
                tvc = tvc.list,
                bhazard = Colon$exp_haz, df = 6)
  
  an <- anova(fit2, fit)
  res1[i] <- an[2, "Pr(>Chisq)"]
  
  fitc2 <- stpm2(Surv(FU_years, status) ~ Age + Gender + Charlson + Metastases + Stage, data = Colon, 
                 tvc = tvc.list,
                 bhazard = Colon$exp_haz, df = 6, cure = T)
  
  an <- anova(fitc2, fitc)
  res2[i] <- an[2, "Pr(>Chisq)"]
}

#Display the LRT p-values
res1
res2


#Plot the time-varying coefficients with confidence intervals
#Set time points and standard newdata matrix, where all covariates are equal to 1
time <- seq(5 / ayear, 16, length.out = 100)
D <- data.frame(Age = 1, Gender = 1, Charlson = 1, Metastases = 1, Stage = 1,
                FU_years = time)


#The time-varying coefficients in the stpm2 model with CIs
X <- rstpm2:::lpmatrix.lm(fit@lm, newdata = D)
pred <- lapply(vars, function(var){
  th.coefs <- grepl(var, names(fit@coef))
  th.bases <- grepl(var, names(fit@coef))
  res <- X[,th.bases] %*% fit@coef[th.coefs]
  vcov <- fit@vcov[th.bases, th.bases]
  vars <- rep(NA, length(time))
  for(i in 1:length(time)){
    vars[i] <- X[i, th.bases] %*% vcov %*% X[i, th.bases] 
  }
  lower <- res - sqrt(vars) * qnorm(0.975)
  upper <- res + sqrt(vars) * qnorm(0.975)
  
  data.frame(Estimate = res, lower = lower, upper = upper)
})

pred <- do.call(rbind, pred)
pred$variable <- rep(vars, each = length(time))
pred$model <- "stpm2"

#The time-varying coefficients in the stpm2 model with CIs
X <- rstpm2:::lpmatrix.lm(fitc@lm, newdata = D)
predc <- lapply(vars, function(var){
  th.coefs <- grepl(var, names(fitc@coef))
  th.bases <- grepl(var, names(fitc@coef))
  res <- X[,th.bases] %*% fitc@coef[th.coefs]
  vcov <- fitc@vcov[th.bases, th.bases]
  vars <- rep(NA, length(time))
  for(i in 1:length(time)){
    vars[i] <- X[i, th.bases] %*% vcov %*% X[i, th.bases] 
  }
  lower <- res - sqrt(vars) * qnorm(0.975)
  upper <- res + sqrt(vars) * qnorm(0.975)
  
  data.frame(Estimate = res, lower = lower, upper = upper)
})

predc <- do.call(rbind, predc)
predc$variable <- rep(vars, each = length(time))
predc$model <- "stpm2-cure"


#The results are combined and the estimates are plotted
plot_data <- rbind(pred, predc)
plot_data$time <- rep(time, length(vars) * 2)
plot_data$variable[plot_data$variable == "Charlson"] <- "CCI \u2265 2"
plot_data$variable[plot_data$variable == "Stage"] <- "UICC stage III-IV"
plot_data$variable[plot_data$variable == "Gender"] <- "Male gender"

add.data <- data.frame(Estimate = rep(c(coefs, coefs.cure), each = length(time)), 
                       model = rep(c("stpm2", "stpm2-cure"), each = length(time) * length(coefs)))

add.data <- data.frame(Estimate = c(coefs, coefs.cure), 
                       model = rep(c("stpm2", "stpm2-cure"), each =length(coefs)))
add.data$variable <- rep(vars, 2)

add.data$variable[add.data$variable == "Charlson"] <- "CCI \u2265 2"
add.data$variable[add.data$variable == "Stage"] <- "UICC stage III-IV"
add.data$variable[add.data$variable == "Gender"] <- "Male gender"

p <- ggplot(plot_data, aes(x = time, y = Estimate, colour = model, group = model)) + geom_line(size = 0.8) + 
  facet_wrap(~variable, ncol = 2, scales = "free_y") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = model, colour = NULL, x = time), 
              alpha = 0.3, show.legend = F) + ylab("Time-varying coefficient") + 
  scale_x_continuous(breaks = seq(0, 16, by = 2)) + 
  xlab("Time since diagnosis (years)") +  
  scale_color_manual(values = pal[2:3], breaks = c("stpm2", "stpm2-cure")) + 
  scale_fill_manual(values = pal[2:3], breaks = c("stpm2", "stpm2-cure")) + 
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank(), 
                     legend.text=element_text(size=15), 
                     axis.title=element_text(size=17),
                     strip.text = element_text(size=15), 
                     axis.text = element_text(size = 13),
                     legend.key.size = unit(2,"line")) + 
  geom_hline(data = add.data, aes(yintercept = Estimate, colour = model), size = 0.8, linetype = "dashed")


#Cairo_pdf is used to get the unichar correctly
cairo_pdf(file.path(fig.out, "coveffectstime.pdf"), width = 9.6, height = 8.4)
print(p)
dev.off()



############Cure rate estimation#################

#Fit the parametric mixture cure model with time-varying covariate effects (may take a long time)
filename <- file.path(data.out, "ColonMixModel3.RData")
if(file.exists(filename)){
  load(filename)
} else {
  fitmix <- GenFlexCureModel(Surv(FU_years, status) ~ ns(Age, 2) + Gender + Charlson + Metastases + Stage, data = Colon, 
                             bhazard = Colon$exp_haz, 
                             cr.formula = ~ ns(Age, 2) + Gender + Charlson + Metastases + Stage, 
                             tvc = list(Gender = 2, Charlson = 2, Metastases = 2, Stage = 2), 
                             smooth.formula = ~ns(Age, 2):nsx(log(FU_years), df = 2) + nsx(log(FU_years), df = 5),
                             control = list(maxit = 100000)) 
  save(fitmix, file = filename)
}

#Fit the stpm2-cure model with a spline on the age-variable
fitc <- stpm2(Surv(FU_years, status) ~ ns(Age, df = 2) + Gender + Charlson + Metastases + Stage, data = Colon, 
              tvc = list(Gender = 3, Charlson = 3, Metastases = 3, Stage = 3), 
              smooth.formula = ~ns(Age, df = 2):nsx(FU_years, df = 3, cure = T) + nsx(log(FU_years), df = 6, cure = T),
              bhazard = Colon$exp_haz)

#Set newdata from which to do predictions
Age <- c(50, 60, 70, 80)
Charlson <- c(0, 1)
Metastases <- 0
Gender <- 0
Stage <- c(0, 1)
newdata <- expand.grid(Age = Age, Charlson = Charlson, Metastases = Metastases, Gender = Gender, Stage = Stage)

#Do predictions from the two cure models
pred <- predict(fitmix, newdata = newdata, type = "curerate")
pred <- do.call(rbind, pred)

plot_data1 <- cbind(pred, newdata)
plot_data1$model <- "Mixture cure model"

pred2 <- predict(fitc, newdata = cbind(newdata, FU_years = 20), se.fit = T)
plot_data2 <- cbind(pred2, newdata)
plot_data2$model <- "stpm2-cure"

#Combine results
plot_data <- rbind(plot_data1, plot_data2)
plot_data$y <- rep(1:length(Age), nrow(plot_data) / length(Age))
plot_data$y[plot_data$model == "stpm2-cure"] <- plot_data$y[plot_data$model == "stpm2-cure"] - 0.1
plot_data$y[plot_data$model == "Mixture cure model"] <- plot_data$y[plot_data$model == "Mixture cure model"] + 0.1
plot_data$Metastases <- factor(plot_data$Metastases, levels = c(0, 1), labels = c("No metastases", "Metastases"))
plot_data$Stage <- factor(plot_data$Stage, levels = c(0, 1), labels = c("UICC Stage I-II", "UICC Stage III-IV"))
plot_data$Charlson <- factor(plot_data$Charlson, levels = c(0, 1), labels = c("CCI < 2", "CCI \u2265 2"))


xlabs <- as.character(Age)
names(xlabs) <- 1:length(Age)

#Plot the estimated cure rates
p <- ggplot(plot_data, aes(y = Estimate, x = y, colour = model)) + 
  geom_point() + facet_grid(. ~ Charlson + Stage) + 
  xlab("Age at diagnosis") + ylab("Cure proportion") + 
  geom_segment(aes(y = lower, yend = upper, x = y, xend = y), linetype = "dashed") + theme_bw() + 
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text=element_text(size=15), 
        axis.title=element_text(size=17),
        strip.text = element_text(size=15), 
        axis.text = element_text(size = 13), 
        legend.key.size = unit(3,"line")) + 
  scale_x_continuous(breaks = 1:length(Age), labels = xlabs) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0.2, 1)) +  
  scale_colour_grey()


#Cairo_pdf is used to get the unichar correctly
cairo_pdf(file.path(fig.out, "CureRateColon.pdf"), width = 10, height = 7)
print(p)
dev.off()



#Compare 5-year relative survival probabilities


#Set newdata from which to do predictions
Age <- c(50, 60, 70, 80)
Charlson <- c(0, 1)
Metastases <- 0
Gender <- 0
Stage <- c(0, 1)
newdata <- expand.grid(Age = Age, Charlson = Charlson, Metastases = Metastases, Gender = Gender, Stage = Stage)

#Do predictions from the two cure models
pred <- predict(fitmix, newdata = newdata, type = "surv", time = 5)
pred <- do.call(rbind, pred)

plot_data1 <- cbind(pred, newdata)
plot_data1$model <- "Mixture cure model"

pred2 <- predict(fitc, newdata = cbind(newdata, FU_years = 5), se.fit = T)
plot_data2 <- cbind(pred2, newdata)
plot_data2$model <- "stpm2-cure"

#Combine results
plot_data <- rbind(plot_data1, plot_data2)
plot_data$y <- rep(1:length(Age), nrow(plot_data) / length(Age))
plot_data$y[plot_data$model == "stpm2-cure"] <- plot_data$y[plot_data$model == "stpm2-cure"] - 0.1
plot_data$y[plot_data$model == "Mixture cure model"] <- plot_data$y[plot_data$model == "Mixture cure model"] + 0.1
plot_data$Metastases <- factor(plot_data$Metastases, levels = c(0, 1), labels = c("No metastases", "Metastases"))
plot_data$Stage <- factor(plot_data$Stage, levels = c(0, 1), labels = c("UICC Stage I-II", "UICC Stage III-IV"))
plot_data$Charlson <- factor(plot_data$Charlson, levels = c(0, 1), labels = c("CCI < 2", "CCI \u2265 2"))


xlabs <- as.character(Age)
names(xlabs) <- 1:length(Age)

#Plot the estimated cure rates
p <- ggplot(plot_data, aes(y = Estimate, x = y, colour = model)) + 
  geom_point() + facet_grid(. ~ Charlson + Stage) + 
  xlab("Age at diagnosis") + ylab("Cure proportion") + 
  geom_segment(aes(y = lower, yend = upper, x = y, xend = y), linetype = "dashed") + theme_bw() + 
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text=element_text(size=15), 
        axis.title=element_text(size=17),
        strip.text = element_text(size=15), 
        axis.text = element_text(size = 13), 
        legend.key.size = unit(3,"line")) + 
  scale_x_continuous(breaks = 1:length(Age), labels = xlabs) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0.2, 1)) +  
  scale_colour_grey()


#Cairo_pdf is used to get the unichar correctly
cairo_pdf(file.path(fig.out, "RelSurv5Colon.pdf"), width = 10, height = 7)
print(p)
dev.off()
