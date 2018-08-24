load("K:/AUH-Databse_Projekter/009_DLBCL_vs_normal_followup/GeneratedData/surv_data.RData")

surv_data$status2 <- as.numeric(surv_data$death_cause_group != "Lymphoma") + 1
surv_data$status2[is.na(surv_data$status2)] <- 0

fit <- rs.surv(Surv(FU_days, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_dato), 
               data = surv_data, ratetable = survexp.dk)
plot(fit)

levels <- c("Censored", "Lymphoma", "Cardiovascular disease", "Other cancers", 
            "Other causes", "No information/\nno contact to hospital")


surv_data$status3 <- as.numeric(factor(surv_data$death_cause_group_new, levels = levels)) - 1

cif <- etmCIF(Surv(FU, status3 != 0) ~ 1, data = surv_data, 
              etype = status3, failcode = 1)

plot(cif, which.cif = 1:5)
sum <- summary(cif)

nrows <- sapply(sum[[1]], nrow)
res.cif <- do.call(rbind, sum[[1]])
rownames(res.cif) <- NULL
res.cif$cause <- factor(rep(levels[-1], nrows), levels = levels)

pdf(file.path(fig.out, "CIFAllPatients.pdf"), width = 8, height = 6)
ggplot(res.cif, aes(x = time, y = P, colour = cause)) + geom_step() + ylim(0, 0.2) +
  ylab("Cumulative incidence") + xlab("Post-treatment follow-up (years)") + theme_bw() + 
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text=element_text(size=12), 
        axis.title=element_text(size=14),
        strip.text = element_text(size=12), 
        axis.text = element_text(size = 12))
dev.off()



surv_data$exp_haz <- general.haz(time = "FU_days", age = "age", sex = "sex", year = "diag_dato", data = surv_data)
surv_data$age_group <- cut(surv_data$age_diagnosis, breaks = c(0, 50, 60, 110), 
                           labels = c("-50", "50-60", "60-"))


#time <- seq(0, 12.5, length.out = 100)


# get_data <- function(data){
#   fit <- stpm2(Surv(FU, status) ~ 1, data = data, bhazard = data$exp_haz, df = 3, cure = F)
#   fitc <- stpm2(Surv(FU, status) ~ 1, data = data, bhazard = data$exp_haz, df = 3, cure = T)
#   cif <- etmCIF(Surv(FU, status2 != 0) ~ 1, data = data, etype = status2, failcode = 1)
#   
#   res <- calc.Crude(fit, time = time, rmap = list(year = diag_dato), var.type = "n")
#   res2 <- calc.Crude(fitc, time = time, rmap = list(year = diag_dato), var.type = "n")
#   
#   
#   ci_cancer <- summary(cif)[[1]]$`CIF 1`[, c(1,2,4,5)]
#   ci_cancer <- rbind(c(0, 0, 0, 0), ci_cancer)
# 
#   probs_cancer <- c(res[[1]]$Estimate, res2[[1]]$Estimate)
#   model.types <- c("Relative survival", "Relative survival - cure", "Death cause")
#   
#   D1 <- data.frame(P = probs_cancer, time = rep(time, 2), lower = NA, upper = NA)
#   D1 <- rbind(D1, ci_cancer)
#   D1$cause <- "Cancer"
#   D1$method <- rep(model.types, c(length(time), length(time), nrow(ci_cancer)))
# 
#   
#   res <- calc.Crude(fit, time = time, rmap = list(year = diag_dato), var.type = "n", type = "other")
#   res2 <- calc.Crude(fitc, time = time, rmap = list(year = diag_dato), var.type = "n", type = "other")
#   ci_other <- summary(cif)[[1]]$`CIF 2`[, c(1,2,4,5)]
#   ci_other <- rbind(c(0, 0, 0, 0), ci_other)
# 
#   probs_other <- c(res[[1]]$Estimate, res2[[1]]$Estimate)
#   D2 <- data.frame(P = probs_other, time = rep(time, 2), lower = NA, upper = NA)
#   D2 <- rbind(D2, ci_other)
#   D2$cause <- "Other causes"
#   D2$method <- rep(model.types, c(length(time), length(time), nrow(ci_other)))
#   
#   rbind(D1, D2)
# }
# 
# 
# res1 <- get_data(surv_data[surv_data$age_group == "-50",])
# res2 <- get_data(surv_data[surv_data$age_group == "50-60",])
# res3 <- get_data(surv_data[surv_data$age_group == "60-",])
# 
# 
# res <- rbind(res1, res2, res3)
# res$age_group <- rep(levels(surv_data$age_group), c(nrow(res1), nrow(res2), nrow(res3)))
# 
# p <- ggplot(res[res$method != "Death cause",], aes(x = time, y = P, colour = method)) + geom_line(size = 1) + 
#   facet_grid(cause ~ age_group) + ylim(0, 0.5) + xlim(0, 10) + theme_bw() + 
#   ylab("Cumulative incidence") + xlab("Post-treatment follow-up (years)") + 
#   geom_step(data = res[res$method == "Death cause",], aes(x = time, y = P, colour = method)) + 
#   geom_step(data = res[res$method == "Death cause",], aes(x = time, y = lower, colour = method), linetype = "dashed", size = 0.4) + 
#   geom_step(data = res[res$method == "Death cause",], aes(x = time, y = upper, colour = method), linetype = "dashed", size = 0.4) + 
#   scale_colour_manual(values = c("Death cause" = "black", 
#                                  "Relative survival" = "chocolate2", 
#                                  "Relative survival - cure" = "darkolivegreen3")) + 
#   theme(legend.position = "bottom", legend.title = element_blank())
# 
# pdf(file.path(fig.out, "ExcessDeathCause.pdf"), width = 10, height = 6)
# print(p)
# dev.off()





##Covariate analysis

# surv_data_comp <- surv_data[complete.cases(surv_data[, c("ECOG")]),]
# fit1 <- stpm2(Surv(FU, status) ~ ECOG, data = surv_data_comp, bhazard = surv_data_comp$exp_haz, tvc = list(ECOG = 3))
# fit2 <- stpm2(Surv(FU, status) ~ ECOG, data = surv_data_comp, 
#               bhazard = surv_data_comp$exp_haz, tvc = list(ECOG = 3), cure = T)
# 
# 
# fit <- rs.surv(Surv(FU_days, status) ~ ECOG + ratetable(age = age, sex = sex, year = diag_dato), 
#                data = surv_data, ratetable = survexp.dk)
# fit$time <- fit$time / ayear
# plot(fit)
# plot(fit1, newdata = data.frame(ECOG = 0), line.col = 2, add = T)
# plot(fit1, newdata = data.frame(ECOG = 1), line.col = 2, add = T)
# plot(fit2, newdata = data.frame(ECOG = 0), line.col = 3, add = T)
# plot(fit2, newdata = data.frame(ECOG = 1), line.col = 3, add = T)
# 
# surv_data_comp <- surv_data[complete.cases(surv_data[, c("age_diagnosis", "ECOG", "stage", "IPI")]),]
# fit1 <- stpm2(Surv(FU, status) ~ age_diagnosis + ECOG + stage + IPI, 
#               data = surv_data_comp, bhazard = surv_data_comp$exp_haz)
# fit2 <- stpm2(Surv(FU, status) ~ age_diagnosis + ECOG + stage + IPI, 
#               data = surv_data_comp, bhazard = surv_data_comp$exp_haz, cure = T)
# 
# sum1 <- summary(fit1)
# sum2 <- summary(fit2)
# 
# cbind(sum1@coef[, c(1,4)], sum2@coef[, c(1,4)])
# 
# 
# fit <- stpm2(Surv(FU, status) ~ age_diagnosis, data = surv_data[which(surv_data$ECOG == 1),], 
#              bhazard = surv_data$exp_haz[which(surv_data$ECOG == 1)], tvc = list(age_diagnosis = 3))
# fitc <- stpm2(Surv(FU, status) ~ age_diagnosis, data = surv_data[which(surv_data$ECOG == 1),], 
#              bhazard = surv_data$exp_haz[which(surv_data$ECOG == 1)], cure = T, tvc = list(age_diagnosis = 3))
# 
# 
# time <- seq(min(surv_data$FU), 10, length.out = 100)
# D <- data.frame(age_diagnosis = 50, FU = time)
# pred1 <- predict(fit, newdata = D, type = "hr", var = "age_diagnosis", se.fit = T)
# pred2 <- predict(fitc, newdata = D, type = "hr", var = "age_diagnosis", se.fit = T)
# 
# head(pred1)
# plot(Estimate ~ time, data = pred1, type = "l", ylim = c(min(pred2$Estimate), max(pred1$Estimate)))
# lines(Estimate ~ time, data = pred2, col = 2)
# 
# summary(fit)
# summary(fitc)

###Investigate hazards#######
#We model the whole age effect with covariates because there is not enough events in the young group.
#This model is unstable in this case and the standard errors become very large or NA. 
#See last analysis for confirmation

surv_data$gender <- ifelse(surv_data$sex == "female", 0, 1)
surv_data$caltime <- as.numeric(surv_data$diag_dato - min(surv_data$diag_dato)) / ayear

time <- seq(0.01, 12, length.out = 300)
time.df <- 2
df.smooth <- 4

fitdc <- stpm2(Surv(FU, status2 == 1) ~ ns(age_diagnosis, df = 2) + gender + ns(caltime, df = 2), data = surv_data, df = df.smooth, 
               tvc.formula = ~ns(age_diagnosis, df = 2):ns(log(FU), df = time.df) + gender:ns(log(FU), df = time.df) + 
                 ns(caltime, df = 2):ns(log(FU), df = time.df))

fitdc.cure <- stpm2(Surv(FU, status2 == 1) ~ ns(age_diagnosis, df = 2) + gender + ns(caltime, df = 2), data = surv_data, 
                    tvc.formula = ~ns(age_diagnosis, df = 2):nsx(log(FU), df = time.df, cure = T) + gender:nsx(log(FU), df = time.df, cure = T) + 
                      ns(caltime, df = 2):nsx(log(FU), df = time.df, cure = T), 
                    smooth.formula = ~nsx(log(FU), df = df.smooth, cure = T))

fitdc2 <- stpm2(Surv(FU, status2 == 2) ~ ns(age_diagnosis, df = 2) + gender + ns(caltime, df = 2), data = surv_data, df = df.smooth, 
                tvc.formula = ~ns(age_diagnosis, df = 2):nsx(log(FU), df = time.df) + gender:nsx(log(FU), df = time.df) + 
                  ns(caltime, df = 2):nsx(log(FU), df = time.df))

fitex <- stpm2(Surv(FU, status) ~ ns(age_diagnosis, df = 2) + gender + ns(caltime, df = 2), 
               data = surv_data, df = df.smooth, bhazard = surv_data$exp_haz, 
               tvc.formula = ~ns(age_diagnosis, df = 2):nsx(log(FU), df = time.df) + gender:nsx(log(FU), df = time.df) + 
                 ns(caltime, df = 2):nsx(log(FU), df = time.df))

fitex.cure <- stpm2(Surv(FU, status) ~ ns(age_diagnosis, df = 2) + gender + ns(caltime, df = 2), 
                    data = surv_data, df = df.smooth, bhazard = surv_data$exp_haz, 
                    tvc.formula = ~ns(age_diagnosis, df = 2):nsx(log(FU), df = time.df, cure = T) + gender:nsx(log(FU), df = time.df, cure = T) + 
                      ns(caltime, df = 2):nsx(log(FU), df = time.df, cure = T), cure = T)

# fitex <- stpm2(Surv(FU, status) ~ ns(age_diagnosis, df = 2), data = surv_data, df = 4, bhazard = surv_data$exp_haz, 
#                tvc.formula = ~ns(age_diagnosis, df = 2):nsx(FU, df = 3, cure = F), cure = F)
# 
# fitex.cure <- stpm2(Surv(FU, status) ~ ns(age_diagnosis, df = 2), data = surv_data, df = 4, bhazard = surv_data$exp_haz, 
#                     tvc.formula = ~ns(age_diagnosis, df = 2):nsx(FU, df = 3, cure = T), cure = T)


models <- c("Death cause", "Death cause - cure", "Relative survival", "Relative survival - cure")
date.eval <- as.Date("2008-01-01")
cal.eval <- as.numeric(date.eval - min(surv_data$diag_dato)) / ayear

get_haz <- function(age){
  newdata <- data.frame(age_diagnosis = age, gender = 0, caltime = cal.eval, FU = time)
  preddc <- predict(fitdc, newdata = newdata, type = "hazard", se.fit = T)
  preddc.cure <- predict(fitdc.cure, newdata = newdata, type = "hazard", se.fit = T)
  predex <- predict(fitex, newdata = newdata, type = "hazard", se.fit = T)
  predex.cure <- predict(fitex.cure, newdata = newdata, type = "hazard", se.fit = T)
  pred <- rbind(preddc, preddc.cure, predex, predex.cure)
  pred$method <- rep(models, each = length(time))
  pred$time <- rep(time, length(models))
  pred
}

ages <- c(50, 60, 70)
preds <- lapply(ages, get_haz)
preds <- do.call(rbind, preds)
preds$age <- factor(rep(paste0(ages, " years of age"), each = length(time) * length(models)))

#Colours for plotting
breaks <- c("Relative survival - cure", "Death cause", "Death cause - cure","Relative survival")
cols <- c("black", "dodgerblue3","chocolate2","darkolivegreen3")
names(cols) <- breaks

p <- ggplot(preds, aes(x = time, y = Estimate, colour = method)) + geom_line(size = 1) + 
  facet_grid(~ age) + coord_cartesian(ylim=c(0, 0.1)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = method, colour = NULL, x = time), 
              alpha = 0.3, show.legend = F) + theme_bw() + 
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text=element_text(size=15), 
        axis.title=element_text(size=17),
        strip.text = element_text(size=15), 
        axis.text = element_text(size = 13)) + 
  ylab("Lymphoma-spcific hazard") + xlab("Post-treatment follow-up (years)") + 
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols)

pdf(file.path(fig.out, "DLBCLhazards.pdf"), width = 10, height = 6)
print(p)
dev.off()

preds2 <- preds[preds$method == "Death cause",]
preds2 <- preds[preds$method == "Relative survival - cure",]

p <- ggplot(preds, aes(x = time, y = Estimate, colour = method)) + geom_line(size = 1) + 
  facet_grid(~ age) + coord_cartesian(ylim=c(0, 0.12)) + 
  geom_ribbon(data = preds2, aes(ymin = lower, ymax = upper, fill = method, colour = NULL, x = time), 
              alpha = 0.3, show.legend = F) + theme_bw() + 
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text=element_text(size=15), 
        axis.title=element_text(size=17),
        strip.text = element_text(size=15), 
        axis.text = element_text(size = 13), 
        legend.key.size = unit(2,"line")) + 
  ylab("Lymphoma-spcific hazard") + xlab("Post-treatment follow-up (years)") + 
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) + 
  scale_x_continuous(breaks = seq(0, 12, 3))

pdf(file.path(fig.out, "DLBCLhazards.pdf"), width = 11, height = 7)
print(p)
dev.off()

calc_cumhaz <- function(time, haz_cs, surv1, surv2, gaussxw){
  scale <- time / 2
  eval <- rep(NA, length(time))
  for(i in 1:length(time)){
    if(time[i] == 0){
      eval[i] <- 0
    } else {
      points <- scale[i] * (gaussxw$nodes + 1)
      eval_gen <- surv1(points)
      eval_rel <- surv2(points)
      eval_haz <- haz_cs(points)
      eval[i] <- sum(gaussxw$weights * (eval_gen * eval_rel * eval_haz))
    }
  }
  prob <- scale * eval
  prob[time == 0] <- 0
  log(-log(prob))
}


get_cumhaz <- function(age, time){
  newdata <- data.frame(age_diagnosis = age, age = age * ayear, 
                        sex = "female", gender = 0, 
                        diag_dato = date.eval, caltime = cal.eval)
  res1 <- calc.Crude(fitex, newdata = newdata, time = time, rmap = list(year = diag_dato), var.type = "n")
  res2 <- calc.Crude(fitex.cure, newdata = newdata, time = time, rmap = list(year = diag_dato), var.type = "n")
  
  gaussxw <- statmod::gauss.quad(100)
  haz_cs <- function(time) predict(fitdc, newdata = data.frame(age_diagnosis = age, 
                                                               gender = 0, caltime = cal.eval, 
                                                               FU = time), type = "hazard")
  surv1 <- function(time) predict(fitdc, newdata = data.frame(age_diagnosis = age, FU = time, 
                                                              gender = 0, caltime = cal.eval))
  haz_cs.cure <- function(time) predict(fitdc.cure, newdata = data.frame(age_diagnosis = age, FU = time, 
                                                                         gender = 0, caltime = cal.eval), type = "hazard")
  surv1.cure <- function(time) predict(fitdc.cure, newdata = data.frame(age_diagnosis = age, FU = time, 
                                                                        gender = 0, caltime = cal.eval))
  surv2 <- function(time) predict(fitdc2, newdata = data.frame(age_diagnosis = age, FU = time, 
                                                               gender = 0, caltime = cal.eval))
  haz_cs.cure(time[time != 0])
  res <- calc_cumhaz(time = time, haz_cs, surv1, surv2, gaussxw = gaussxw)
  res <- exp(-exp(res))
  res.cure <- calc_cumhaz(time = time, haz_cs.cure, surv1.cure, surv2, gaussxw = gaussxw)
  res.cure <- exp(-exp(res.cure))
  probs_cancer <- c(res1[[1]]$Estimate, res2[[1]]$Estimate, res, res.cure)
  
  #D1 <- data.frame(P = probs_cancer, time = rep(time, 3), lower = NA, upper = NA)
  #D1 <- rbind(D1, ci_cancer)
  #D1$cause <- "Cancer"
  #D1$method <- rep(model.types, c(length(time), length(time), nrow(ci_cancer)))
  
  
  res1 <- calc.Crude(fitex, newdata = newdata, time = time, rmap = list(year = diag_dato), var.type = "n", type = "other")
  res2 <- calc.Crude(fitex.cure, newdata = newdata, time = time, rmap = list(year = diag_dato), var.type = "n", type = "other")
  haz_cs <- function(time) predict(fitdc2, newdata = data.frame(age_diagnosis = age, FU = time, 
                                                                gender = 0, caltime = cal.eval), 
                                   type = "hazard")
  res <- calc_cumhaz(time = time, haz_cs, surv1, surv2, gaussxw = gaussxw)
  res <- exp(-exp(res))
  res.cure <- calc_cumhaz(time = time, haz_cs = haz_cs, surv1 = surv1.cure, surv2 = surv2, gaussxw = gaussxw)
  res.cure <- exp(-exp(res.cure))
  probs_other <- c(res1[[1]]$Estimate, res2[[1]]$Estimate, res, res.cure)
  model.types <- c("Relative survival", "Relative survival - cure", "Death cause", "Death cause - cure")
  
  data.frame(probs = c(probs_cancer, probs_other), 
             cause = rep(c("Lymphoma", "Other causes"), each = length(probs_cancer)), 
             method = rep(rep(model.types, each = length(time)), 2), 
             age = age, time = rep(time, 2))
}

time <- seq(0, 12, length.out = 100)
res1 <- get_cumhaz(50, time = time)
res2 <- get_cumhaz(60, time = time)
res3 <- get_cumhaz(70, time = time)


res <- rbind(res1, res2, res3)
res$age <- paste0(res$age, " years of age")

p <- ggplot(res, aes(x = time, y = probs, colour = method)) + geom_line(size = 1) + 
  facet_grid(cause ~ age) + ylim(0, 0.3) + theme_bw() + 
  ylab("Cumulative incidence") + xlab("Post-treatment follow-up (years)") + 
  scale_colour_manual(values = cols) + 
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text=element_text(size=15), 
        axis.title=element_text(size=17),
        strip.text = element_text(size=15), 
        axis.text = element_text(size = 13), 
        legend.key.size = unit(2,"line")) + 
  scale_x_continuous(breaks = seq(0, 12, 3))


pdf(file.path(fig.out, "ExcessDeathCausePara.pdf"), width = 11, height = 7)
print(p)
dev.off()







# time <- seq(0.05, 10, length.out = 100)
# fitdc <- stpm2(Surv(FU, status2 == 1) ~ 1, data = surv_data, df = 5)
# fitex <- stpm2(Surv(FU, status) ~ 1, data = surv_data, df = 5, bhazard = surv_data$exp_haz, cure = F)
# 
# newdata <- data.frame(age_diagnosis = 50, FU = time)
# preddc <- predict(fitdc, newdata = newdata, type = "hazard", se.fit = T)
# predex <- predict(fitex, newdata = newdata, type = "hazard", se.fit = T)
#   
# pred <- rbind(preddc, predex)
# pred$method <- rep(c("Death cause", "Relative survival"), each = length(time))
# pred$time <- rep(time, 2)
# 
# 
# p <- ggplot(pred, aes(x = time, y = Estimate, colour = method)) + geom_line() + 
#   geom_ribbon(aes(ymin = lower, ymax = upper, fill = method, colour = NULL, x = time), 
#               alpha = 0.3, show.legend = F) + theme_bw() + 
#   theme(legend.position = "bottom", legend.title = element_blank()) + 
#   ylab("Hazard") + xlab("Post-treatment follow-up (years)") + 
#   scale_color_manual(values = c("black", "darkolivegreen3"), breaks = c("Death cause", "Relative survival")) + 
#   scale_fill_manual(values = c("black", "darkolivegreen3"), breaks = c("Death cause", "Relative survival"))
# 
# pdf(file.path(fig.out, "DLBCLhazardAll.pdf"), width = 8, height = 6)
# print(p)
# dev.off()



###Stratified analyses

time <- seq(0.05, 10, length.out = 100)
fit_models <- function(surv_data){
  fitdc <- stpm2(Surv(FU, status2 == 1) ~ 1, data = surv_data, df = 3)
  fitex <- stpm2(Surv(FU, status) ~ 1, data = surv_data, df = 3, bhazard = surv_data$exp_haz)
  preddc <- predict(fitdc, newdata = data.frame(x = 1, FU = time), type = "hazard", se.fit = T)
  predex <- predict(fitex, newdata = data.frame(x = 1, FU = time), type = "hazard", se.fit = T)
  preddc$method <- "Death cause"
  predex$method <- "Relative survival"
  pred <- rbind(preddc, predex)
}


res1 <- fit_models(surv_data[surv_data$age_group == "-50",])
res2 <- fit_models(surv_data[surv_data$age_group == "50-60",])
res3 <- fit_models(surv_data[surv_data$age_group == "60-",])


preds <- rbind(res1, res2, res3)
lvs <- levels(surv_data$age_group)
preds$age_group <- factor(rep(lvs, each = length(time)), levels = lvs)
preds$time <- rep(time, 2 * length(lvs))



p <- ggplot(preds, aes(x = time, y = Estimate, colour = method)) + geom_line() + 
  facet_grid(~ age_group) + coord_cartesian(ylim=c(0, 0.1)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = method, colour = NULL, x = time), 
              alpha = 0.3, show.legend = F) + theme_bw() + 
  theme(legend.position = "bottom", legend.title = element_blank()) + 
  ylab("Hazard") + xlab("Post-treatment follow-up (years)") + 
  scale_color_manual(values = c("black", "darkolivegreen3"), breaks = c("Death cause", "Relative survival")) + 
  scale_fill_manual(values = c("black", "darkolivegreen3"), breaks = c("Death cause", "Relative survival"))
p

