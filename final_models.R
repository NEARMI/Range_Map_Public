library(ggplot2)
library(mgcv)
library(rgdal)
library(dplyr)
library(raster)
library(voxel)
library(gridExtra)
library(terra)
library(tidyverse)

'%notin%' <- Negate('%in%')
### THESE ARE NOT REAL LOCATIONS >>> THIS WILL NOT RECOVER RESULTS FROM PAPER###
sh.plot <- read.csv("non_geo_data.csv")
cov <- read.csv("covariate_non_geo.csv")

source("Graphing_Set_up.R")


## change aspect into radians and then into Northing (0 at due E or due W so cos) and Easting
## standardize HLI
## add 0/1s instead of counts
sh.plot <- sh.plot %>% mutate(asp_rad = pi*aspect/180) %>%
  mutate(north = cos(asp_rad), east = sin(asp_rad)) %>%
  mutate(std_HLI = scale(HLI)[,1], std_IMI = scale(IMI)[,1])

###########################    Set up cross validation   #########################
chy <- (max(sh.plot$y) -  min(sh.plot$y))/3
chx <- (max(sh.plot$x) -  min(sh.plot$x))/3

wiggles <- c(25, 50, 200)
scen <- c(1,2,3)

mod_out <- vector("list",3)
mod_out[[1]] <- vector("list", 3)
mod_out[[1]][[1]] <- vector("list", 3)
mod_out[[1]][[2]] <- vector("list", 3)
mod_out[[1]][[3]] <- vector("list", 3)

mod_out[[2]] <- vector("list", 3)
mod_out[[2]][[1]] <- vector("list", 3)
mod_out[[2]][[2]] <- vector("list", 3)
mod_out[[2]][[3]] <- vector("list", 3)


mod_out[[3]] <- vector("list", 3)
mod_out[[3]][[1]] <- vector("list", 3)
mod_out[[3]][[2]] <- vector("list", 3)
mod_out[[3]][[3]] <- vector("list", 3)

###########################   Cross Validation  #########################
## takes about 1.5 hours to fit ##


system.time(
for(t in 1: length(scen)){
  print("tt")
for(j in 1:length(scen)){
  print("jj")
  ## leave out one- ninth of the area 
  ## summarize each site to mean temp and mean ppt -- other values are the same by line
    ## this is because I want to predict by site not by visit
  sh.plot1 <- sh.plot %>%
    filter(x >= min(sh.plot$x) + (j-1)*chx & x < min(sh.plot$x) + j*chx  &
             y >= min(sh.plot$y) + (t-1)*chy & y < min(sh.plot$y) + t*chy) %>%
    group_by(ID) %>%
    summarize(temp10 = mean(temp10, na.rm =T), ppt = mean(ppt, na.rm =T), std_HLI = mean(std_HLI), elev = mean(elev),
              x = mean(x), y = mean(y), north = mean(north), east = mean(east), std_IMI= mean(std_IMI),
              shc = sum(shencount), cic = sum(cincount)) %>%
    ungroup()  %>% mutate(stt = ifelse(shc == 0, 0, 1), ctt = ifelse(cic == 0, 0, 1))
    
    
  sh.plot2 <- sh.plot %>% filter(ID %notin% sh.plot1$ID)
  
for(i in 1:length(wiggles)){
  if(wiggles[i] < length(unique(sh.plot2$ID)) ){
 shen_gam = gam(shencount ~-1+ s(elev, k = 6) 
                + s(std_IMI, k=6) 
                + s(std_HLI, k = 6) 
                + s(ppt, temp10, bs='ds', k =14, m=2)
                + s(x, y, bs='ds', k =wiggles[i], m=2) 
                + s(ID,bs = "re")
               ,
               data=sh.plot2, family = "ziP",offset = log(area/10000), method = 'REML',
               drop.unused.levels = FALSE)
#summary(shen_gam)
mod_out[[i]][[t]][[j]] <-  shen_gam
#shen_gam <- mod_out[[i]][[t]][[j]]


## predict for left out sites
shen.est <- predict(shen_gam, newdata=sh.plot1, se.fit = T )


## calculate root mean square error 
rmse = sqrt(sum((sh.plot1$shc - exp(shen.est$fit))^2/nrow(sh.plot1)))

## calculate coverage
Upper <- exp(shen.est$fit + 2*shen.est$se.fit)
Lower <- exp(shen.est$fit - 2*shen.est$se.fit)
ccc <- ifelse(sh.plot1$shc > Lower & sh.plot1$shc < Upper, 1, 0)

cov_prop <- mean(ccc, na.rm =T)

## store 

if( t == 1 & j == 1 & i == 1){
  out <- cbind(rmse, cov_prop, wiggles[i], j, t, length(unique(sh.plot1$ID)))
} else {
  new <- cbind(rmse, cov_prop, wiggles[i], j,t, length(unique(sh.plot1$ID)))
  out <- rbind(out, new)
} 
}else {
  mod_out[[i]][[t]][[j]] < 0
  rmse = NA
  cov_prop = NA
  new <- cbind(rmse, cov_prop, wiggles[i], j, t, length(unique(sh.plot1$ID)))
  out <- rbind(out, new)
}
  

## compare to data 1000 times
#sitet <- seq(1,length(unique(sh.plot1$ID)))

#for(p in 1:length(sitet)){
#temp <- rnorm(1000, shen.est$fit[sitet[p]], shen.est$se.fit[sitet[p]])  %>% plogis() %>% rbinom(length(.), 1, .)
# prop_correct = length(which(temp == sh.plot1$stt[p]))/1000


 
# if(p == 1 & t == 1 & j == 1 & i == 1){
 #  out <- cbind(prop_correct,  wiggles[i], j, t, length(unique(sh.plot1$ID)))
# } else {
 #  new <- cbind(prop_correct,  wiggles[i], j,t, length(unique(sh.plot1$ID)))
 #  out <- rbind(out, new)
# }
 
 
 
#} ## close p

  
  } ## close i
} ## close j
} ## close t

)



saveRDS(out, file = "prop_correct_ZiPIm_aug9.RDS")
saveRDS(mod_out, file = "model_output_ZiPIm_aug9.RDS") 


out <- as.data.frame(out)
names(out) <- c("rmse","cov_prop","wiggles", "j", "t", "rm_sites")


outpred <-out %>% ungroup() %>% filter(rm_sites > 0) %>% group_by(wiggles, rm_sites) %>% summarize(tt = mean(rmse), cc = mean(cov_prop))

out %>% ungroup() %>% filter(rm_sites > 0) %>% group_by(wiggles) %>% summarize(tt = mean(rmse, na.rm =T), cc = mean(cov_prop, na.rm=T))

## remove largest RMSE 
out %>% ungroup() %>% filter(rm_sites > 0) %>% group_by(wiggles) %>% filter(rmse < max(rmse)) %>%
  summarize(tt = mean(rmse, na.rm =T), cc = mean(cov_prop, na.rm=T))



rm(out)
gc()


## repeat for P. cinereus


## loop takes about 15 min
system.time(
for(t in 1: length(scen)){
  print("tt")
  for(j in 1:length(scen)){
    print("jj")
    ## leave out one- ninth of the area 
    ## summarize each site to mean temp and mean ppt -- other values are the same by line
    ## this is because I want to predict by site not by visit
    sh.plot1 <- sh.plot %>%
      filter(x >= min(sh.plot$x) + (j-1)*chx & x < min(sh.plot$x) + j*chx  &
               y >= min(sh.plot$y) + (t-1)*chy & y < min(sh.plot$y) + t*chy) %>%
      group_by(ID) %>%
      summarize(temp10 = mean(temp10, na.rm =T), ppt = mean(ppt, na.rm =T), std_HLI = mean(std_HLI), elev = mean(elev),
                x = mean(x), y = mean(y), north = mean(north), east = mean(east),std_IMI= mean(std_IMI),
                shc = sum(shencount), cic = sum(cincount)) %>%
      ungroup() %>% mutate(stt = ifelse(shc == 0, 0, 1), ctt = ifelse(cic == 0, 0, 1))
    
    
    sh.plot2 <- sh.plot %>% filter(ID %notin% sh.plot1$ID)
    
    for(i in 1:length(wiggles)){
        cin_gam = gam(cincount ~ -1+ s(elev, k = 6) 
                    + north
                    +east 
                    +s(std_HLI, k = 6) 
                    + s(ppt, temp10, bs='ds', k =14, m=2)
                    + s(x, y, bs='ds', k =wiggles[i], m=2) + s(ID,bs = "re"),
                    data=sh.plot2, family = "ziP", method = 'REML', offset = log(area/10000),
                    drop.unused.levels = FALSE)
      

      mod_out[[i]][[t]][[j]] <- cin_gam
      
      ## predict for left out sites
      cin.est <- predict(cin_gam, newdata=sh.plot1, se.fit = T, type = "response")
      
      ## calculate root mean square error 
      rmse = sqrt(sum((sh.plot1$cic - (cin.est$fit))^2/nrow(sh.plot1), na.rm = T))
      
      ## calculate coverage
      Upper <- cin.est$fit + 2*cin.est$se.fit
      Lower <- cin.est$fit - 2*cin.est$se.fit
      ccc <- ifelse(sh.plot1$cic > Lower & sh.plot1$cic < Upper, 1, 0)
      
      cov_prop <- mean(ccc, na.rm =T)
      
      ## store 
      
      if( t == 1 & j == 1 & i == 1){
        out <- cbind(rmse, cov_prop, wiggles[i], j, t, length(unique(sh.plot1$ID)))
      } else {
        new <- cbind(rmse, cov_prop,  wiggles[i], j,t, length(unique(sh.plot1$ID)))
        out <- rbind(out, new)
      }
      
      
      
    } ## close i
  } ## close j
} ## close t

)






out <- as.data.frame(out)
names(out) <- c("rmse","cov_prop", "wiggles", "j", "t", "rm_sites")



out %>% filter(rm_sites > 0) %>% group_by(wiggles) %>% summarize(tt = mean(rmse, na.rm =T), cc = mean(cov_prop, na.rm=T))


## remove largest RMSE 
out %>% ungroup() %>% filter(rm_sites > 0) %>% group_by(wiggles) %>% filter(rmse < max(rmse)) %>%
        summarize(tt = mean(rmse, na.rm =T), cc = mean(cov_prop, na.rm=T))





saveRDS(out, file = "prop_correct_cinzipM_aug9.RDS")
saveRDS(mod_out, file = "model_output_cinzipM_aug9.RDS") 

rm(out)
rm(mod_out)
gc()

###########################    FIT GAMS   #########################
#### run full models without any sites removed for final covariate estimates

shen_gam25 = gam(shencount ~-1+ s(elev, k = 6)
              +s(std_IMI, k =6)
               + s(std_HLI, k = 6) 
               + s(ppt, temp10, bs='ds', k =14, m=2)
               + s(x, y, bs='ds', k =25, m=2) 
               + s(ID,bs = "re")
               ,
               data=sh.plot, family = "ziP", method = 'REML', offset = log(area/10000),
               drop.unused.levels = FALSE)


shen_gam50 = gam(shencount ~ -1+ s(elev, k = 6)
                 +s(std_IMI, k =6)
                 + s(std_HLI, k = 6) 
                 + s(ppt, temp10, bs='ds', k =14, m=2)
                 + s(x, y, bs='ds', k =50, m=2) 
                 + s(ID,bs = "re")
                 ,
                 data=sh.plot, family = "ziP", method = 'REML', offset = log(area/10000),
                 drop.unused.levels = FALSE)

## takes about 6 minutes to run
shen_gam200 = gam(shencount ~ -1+s(elev, k = 6)
                 +s(std_IMI, k =6)
                 + s(std_HLI, k = 6) 
                 + s(ppt, temp10, bs='ds', k =14, m=2)
                 + s(x, y, bs='ds', k =200, m=2) 
                 + s(ID,bs = "re")
                 ,
                 data=sh.plot, family = "ziP", method = 'REML', offset = log(area/10000),
                 drop.unused.levels = FALSE)


fin_mod_out <- vector("list",3)
fin_mod_out[[1]] <- shen_gam25
fin_mod_out[[2]] <- shen_gam50
fin_mod_out[[3]] <- shen_gam200
names(fin_mod_out) <- c("psK25", "psK50", "psK200")
saveRDS(fin_mod_out, "Pshen_mod_outzIP_aug9.RDS")

cin_gam25 = gam(cincount ~ -1+s(elev, k = 6)
               
                 + s(std_HLI, k = 6) 
                +s(std_IMI, k =6)
                 + s(ppt, temp10, bs='ds', k =14, m=2)
                 + s(x, y, bs='ds', k =25, m=2) 
                 + s(ID,bs = "re")
                 ,
                 data=sh.plot, family = "ziP", method = 'REML', offset = log(area/10000),
                 drop.unused.levels = FALSE)


cin_gam50 = gam(cincount ~ s(elev, k = 6)
             
                 + s(std_HLI, k = 6) 
                +s(std_IMI, k =6)
                 + s(ppt, temp10, bs='ds', k =14, m=2)
                 + s(x, y, bs='ds', k =50, m=2) 
                 + s(ID,bs = "re")
                 ,
                 data=sh.plot, family = "ziP", method = 'REML', offset = log(area/10000),
                 drop.unused.levels = FALSE)

## takes about 6 minutes to run
cin_gam200 = gam(cincount ~ s(elev, k = 6)
                 # + north
                  #+ east
                  + s(std_HLI, k = 6) 
                 +s(std_IMI, k =6)
                  + s(ppt, temp10, bs='ds', k =14, m=2)
                  + s(x, y, bs='ds', k =200, m=2) 
                  + s(ID,bs = "re")
                  ,
                  data=sh.plot, family = "ziP", method = 'REML', offset = log(area/10000),
                  drop.unused.levels = FALSE)

fin_mod_out <- vector("list",3)
fin_mod_out[[1]] <- cin_gam25
fin_mod_out[[2]] <- cin_gam50
fin_mod_out[[3]] <- cin_gam200
names(fin_mod_out) <- c("pcK25", "pcK50", "pcK200")
saveRDS(fin_mod_out, "cin_mod_out_Zip_aug9.RDS")





###########################   COVARIATE FIGURES  #########################
fin_mod_out <- readRDS("Pshen_mod_outzIP_aug9.RDS")
fin_mod_outC <- readRDS("cin_mod_out_Zip_aug9.RDS")

shen_gam200 <- fin_mod_out[[3]]
cin_gam200 <- fin_mod_outC[[3]]

### equivalent elev effect size plot in ggplot
p_objS <- plot(shen_gam200, residuals = TRUE)
p_obj1 <- p_objS[[1]] # just one smooth so select the first component-- Elevation
sm_df <- as.data.frame(p_obj1[c("x", "se", "fit")])
data_df <- as.data.frame(p_obj1[c("raw", "p.resid")])

## unscale elevation for graphing 

elv <- cov %>% dplyr::select(elev, elev1)


sm2 <- sm_df 
sm2 <- sm2 %>% 
  mutate(unscaled = plyr::mapvalues(round(x,2), from = round(elv$elev,2), to = elv$elev1))
sm2 <- sm2 %>% filter(unscaled >200) ## 3 values not matching

(elev <- ggplot(sm2, aes(x = unscaled, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Elevation (m)", y = "Effect Size") +ggtitle("A. P. Shenandoah") )#+ ylim(-40, 5) )

p_obj2 <- p_objS[[2]] # just one smooth so select the second component-- IMI
sm_df <- as.data.frame(p_obj2[c("x", "se", "fit")])
data_df <- as.data.frame(p_obj2[c("raw", "p.resid")])


(imi <- ggplot(sm_df, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized IMI", y = "Effect Size") +ggtitle("B.")
    + 
    ylim(-8, 2.5) 
    
    )



p_obj3 <- p_objS[[3]] # just one smooth so select the third component-- HLI
sm_df <- as.data.frame(p_obj3[c("x", "se", "fit")])
data_df <- as.data.frame(p_obj3[c("raw", "p.resid")])


(hli <- ggplot(sm_df, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized HLI", y = "Effect Size") +ggtitle("C.") + 
    ylim(-3, 2))


p_objC <- plot(cin_gam200, residuals = TRUE)
p_obj1C <- p_objC[[1]] # just one smooth so select the first component-- Elevation
sm_dfC <- as.data.frame(p_obj1C[c("x", "se", "fit")])
data_dfC <- as.data.frame(p_obj1[c("raw", "p.resid")])


sm2c <- sm_dfC 
sm2c <- sm2c %>% 
  mutate(unscaled = plyr::mapvalues(round(x,2), from = round(elv$elev,2), to = elv$elev1))
sm2c <- sm2c %>% filter(unscaled >200) ## 3 values not matching

(elevC <- ggplot(sm2c, aes(x = unscaled, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Elevation (m)", y = "Effect Size") +ggtitle("D. P. cinereus")  +ylim(-40,5))

p_obj2C <- p_objC[[2]] # just one smooth so select the second component-- HLI
sm_dfC <- as.data.frame(p_obj2C[c("x", "se", "fit")])
data_dfC <- as.data.frame(p_obj2C[c("raw", "p.resid")])


(hliC <- ggplot(sm_dfC, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized HLI", y = "Effect Size") +ggtitle("F.") + 
    ylim(-3, 2) )



p_obj3C <- p_objC[[3]] # just one smooth so select imi
sm_dfC <- as.data.frame(p_obj3C[c("x", "se", "fit")])
data_dfC <- as.data.frame(p_obj3C[c("raw", "p.resid")])


(imiC <- ggplot(sm_dfC, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized IMI", y = "Effect Size") +ggtitle("E.") + 
    ylim(-8, 2.5)
  )


grid.arrange(elev, imi, hli,elevC, imiC, hliC, nrow = 2)



pl4 <- plot(shen_gam200, select = 4, scale = 0, 
            ylab = "Temp", shade = TRUE, shade.col = "gray80", lwd = 2)
ps4 <- pl4[[4]]


ps4q <- as.data.frame(ps4[c("x",  "y")])


tempo <- expand.grid(ps4q$x, ps4q$y) 

temp2 <- cbind(tempo, ps4$fit) %>% rename(x= Var1, y = Var2, fit = "ps4$fit" ) %>%
  filter(!is.na(fit))

#mycol <- c( "blueviolet","blue3" ,"turquoise3", "springgreen4" ,  "springgreen2" ,"aquamarine1",
           # "gray73", "goldenrod4" , "gold3", "goldenrod2" , "orange" ,"orangered")

mycol <- c( "blueviolet","blue3" ,"turquoise3" ,"aquamarine1",
            "gray73", "yellow" , "gold1", "orange" , "orangered")

(stp <- ggplot(temp2, aes(x=x, y=y, z = fit)) + 
    geom_contour_filled(breaks = c(-1, -0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1,1.2)) + 
  xlab("3-day precipitation") + ylab("Temperature") + ggtitle("A. P. shenandoah") +
    scale_fill_manual("Effect Size", values = mycol))
(stp <- stp + theme(legend.position = (c(0.7,0.25)), legend.direction = "horizontal", legend.background = element_blank()) + 
  annotate("text", x= 22, y = 7, label= "Lower count") +
  annotate("text", x= 38, y = 7, label = "Higher count"))



cl4 <- plot(cin_gam200, select = 4, scale = 0, 
            ylab = "Temp", shade = TRUE, shade.col = "gray80", lwd = 2)
cs4 <- cl4[[4]]


cs4q <- as.data.frame(cs4[c("x",  "y")])


tempoc <- expand.grid(cs4q$x, cs4q$y) 

temp2c <- cbind(tempoc, cs4$fit) %>% rename(x= Var1, y = Var2, fit = "cs4$fit" ) %>%
  filter(!is.na(fit))

#geom_contour_filled(breaks = c(-1.4, -1.2,-1,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6))

(ctp <- ggplot(temp2c, aes(x=x, y=y, z = fit)) +  geom_contour_filled(breaks = c(-1,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1,1.2)) + 
    xlab("3-day precipitation") + ylab("Temperature") + ggtitle("B. P. cinereus") +
    scale_fill_manual("Effect Size", values = mycol)+ theme(legend.position = "none"))
 

gridExtra::grid.arrange(stp, ctp)



###########################   MAKE MAPS   #########################

qq <- readRDS("Pshen_mod_outzIP_aug9.RDS")
zz <-readRDS("cin_mod_out_Zip_aug9.RDS")

cin_gam200 <- zz[[3]]

shen_gam200 <- qq[[3]]
rm(zz, qq)

cov$area <- log(10000)

cov <- subset(cov, !is.na(IMI))

cov <- cov %>% mutate(asp_rad = pi*aspect/180) %>%
  mutate(north = cos(asp_rad), east = sin(asp_rad)) %>%
  mutate(std_HLI = scale(HLI)[,1], std_IMI = scale(IMI)[,1])

pred_weather <- read.csv("weath_cov.csv")

cov <- cbind(cov, pred_weather[,c("ppt3", "temp10")])

cov <- cov %>%
  rename("ppt" = "ppt3")



pred_gamS <- predict(shen_gam200, data.frame(cov, area = rep(10000,nrow(cov))), se.fit = TRUE)



pg2 <-  data.frame(predS = pred_gamS$fit,seS = pred_gamS$se.fit, 
                        lat = cov$lat, lon = cov$lon,x = cov$UTMx,
                        y = cov$UTMy, elev <- cov$elev, aspect = cov$aspect, 
                        IMI = cov$IMI, std_HLI = cov$std_HLI, left = cov$left_, right = cov$right_,
                        top = cov$top_, bottom = cov$bottom, predS1 = pred_gamS$fit)



pg2$LI <- exp(pg2$predS- 2*pg2$seS)
pg2$UI <- exp(pg2$predS+ 2*pg2$seS)
pg2$status <- "Uncertain"
pg2$status[pg2$LI > 1 | exp(pg2$predS) > 20] <- "Highly Likely"

pg2$status[exp(pg2$predS) <= 1 ] <- "Unlikely"
pg3 <- pg2 %>% mutate(predS = ifelse(predS < -1, -1, predS))


pg3$dummy <- ifelse(exp(pg3$predS) > 1000, "A", ifelse(
  exp(pg3$predS) < 1000 & exp(pg3$predS) > 500, "B", ifelse(
    exp(pg3$predS) < 500 & exp(pg3$predS) > 250, "C", ifelse(
      exp(pg3$predS) < 250 & exp(pg3$predS) > 125, "D", ifelse(
        exp(pg3$predS) < 125 & exp(pg3$predS) > 75, "E", ifelse(
          exp(pg3$predS) < 75 & exp(pg3$predS) > 35, "F", ifelse(
            exp(pg3$predS) < 35 & exp(pg3$predS) > 18, "G", ifelse(
              exp(pg3$predS) < 18 & exp(pg3$predS) > 9,   "H", ifelse(
                exp(pg3$predS) < 9 & exp(pg3$predS) > 4,   "H", ifelse(
                  exp(pg3$predS) < 4 & exp(pg3$predS) > 1,   "I",ifelse(
                    exp(pg3$predS) < 1 & exp(pg3$predS) > 0 , "J", "K"
                  )
                )
              )
            )     
          )
        )
      )  
    )
  )
))




## scale ses to 0,1 for alpha

pg3 <- pg3 %>% mutate(se2 = (seS - min(seS))/ (max(seS)-min(seS)))


pg3$erDum <- ifelse(exp(pg3$seS) > 1000, 0.2, ifelse(
  exp(pg3$seS) < 1000 & exp(pg3$seS) > 100, 0.4, ifelse(
    exp(pg3$seS) < 100 & exp(pg3$seS) > 20, 0.6, 1
  )
))




(presg<- ggplot() +  
    geom_tile(data = pg3, aes(x = x, y = y, fill = dummy, alpha=erDum)) + 
    geom_point(data = sh.plot, aes(x = UTMx, y = UTMy, color = shencount>0 ))+
    scale_color_manual("P. shenandoah \n sampling",
                       labels = c("None found", "At least 1 found"), values = c("grey30","thistle3")) +
    scale_fill_manual("P. shenandoah \n predicted counts",
                      
                      values = cols.d$colo,
                      labels = c(1000, 500, 250, 125,75,35,18,9,4,1)) +
   
    
    ylab("Latitude") + xlab("Longitude") + ggtitle("")  +
    theme(legend.position = c(0.2, 0.8), legend.background = element_blank(), 
          plot.margin = unit(c(0.5,0.6,0.5,0.4), 'cm'),
          legend.box = "horizontal") + scale_alpha(guide = 'none') + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank()) + 
    guides(color = guide_legend(order = 2), fill = guide_legend(order = 1)) )

 


#### P. cinereus 



pred_gamC <- predict(cin_gam200, data.frame(cov, area = rep(10000,nrow(cov))), se.fit = TRUE)



pg2C <-  data.frame(predC = pred_gamC$fit,seC = pred_gamC$se.fit, 
                         lat = cov$lat, lon = cov$lon,x = cov$UTMx,
                         y = cov$UTMy, elev <- cov$elev, aspect = cov$aspect, 
                         IMI = cov$IMI, std_HLI = cov$std_HLI, left = cov$left_, right = cov$right_,
                         top = cov$top_, bottom = cov$bottom, predC1 = pred_gamC$fit)



pg2C$LI <- exp(pg2C$predC- 2*pg2C$seC)
pg2C$UI <- exp(pg2C$predC+ 2*pg2C$seC)
pg2C$status <- "Uncertain"
pg2C$status[pg2C$LI > 1 | exp(pg2C$predC) > 20] <- "Highly Likely"
pg2C$status[exp(pg2C$predC) <= 1 ] <- "Unlikely"
pg3c <- pg2C %>% mutate(predC = ifelse(predC < -1, -1, predC))


pg3c$dummy <- ifelse(exp(pg3c$predC) > 1000, "A", ifelse(
  exp(pg3c$predC) < 1000 & exp(pg3c$predC) > 500, "B", ifelse(
    exp(pg3c$predC) < 500 & exp(pg3c$predC) > 250, "C", ifelse(
     exp(pg3c$predC) < 250 & exp(pg3c$predC) > 125, "D", ifelse(
        exp(pg3c$predC) < 125 & exp(pg3c$predC) > 75, "E", ifelse(
          exp(pg3c$predC) < 75 & exp(pg3c$predC) > 35, "F", ifelse(
           exp(pg3c$predC) < 35 & exp(pg3c$predC) > 18, "G", ifelse(
              exp(pg3c$predC) < 18 & exp(pg3c$predC) > 9,   "H", ifelse(
               exp(pg3c$predC) < 9 & exp(pg3c$predC) > 4,   "H", ifelse(
                  exp(pg3c$predC) < 4 & exp(pg3c$predC) > 1,   "I",ifelse(
                   exp(pg3c$predC) < 1 & exp(pg3c$predC) > 0 , "J", "K"
                  )
               )
              )
            )     
          )
        )
      )  
    )
  )
))



## scale ses to 0,1 for alpha

pg3c <- pg3c %>% mutate(se2 = (seC - min(seC))/ (max(seC)-min(seC)))


pg3c$erDum <- ifelse(exp(pg3c$seC) > 1000, 0.2, ifelse(
  exp(pg3c$seC) < 1000 & exp(pg3c$seC) > 100, 0.4, ifelse(
    exp(pg3c$seC) < 100 & exp(pg3c$seC) > 20, 0.6, 1
  )
))

(presgC<- ggplot() +  
    geom_tile(data = pg3c, aes(x = x, y = y, fill = dummy, alpha=erDum)) + 
    geom_point(data = sh.plot, aes (x = UTMx, y = UTMy, color = cincount>0 ))+
    scale_color_manual("P. cinereus \n sampling",
                       labels = c("None found", "At least 1 found"), values = c("grey30","thistle3")) +
    scale_fill_manual("P. cinereus \n predicted counts",
                      
                      values = cols.d$colo,
                      labels = c(1000, 500, 250, 125,75,35,18,9,4,1)) +
    
    
    ylab("Latitude") + xlab("Longitude") + ggtitle("")  +
    theme(legend.position = c(0.2, 0.8), legend.background = element_blank(), 
          plot.margin = unit(c(0.5,0.6,0.5,0.4), 'cm'),
          legend.box = "horizontal") + scale_alpha(guide = 'none') + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank()) + 
    guides(color = guide_legend(order = 2), fill = guide_legend(order = 1)) )


## combine pshen and pcin onto one map to see overlap

ndat <- pg2 %>% filter(status == "Highly Likely") %>% 
  mutate(spec = "P. shenandoah") %>% 
  rename("count" = "predS") %>% dplyr::select(-c(predS1,seS))
nd <- pg2C %>% filter(status == "Highly Likely") %>%
  mutate(spec = "P. cinereus")%>% 
  rename("count" = "predC")%>% dplyr::select(-c(predC1,seC))
ndat <- rbind(ndat, nd)


ndat$cd <- ifelse(exp(ndat$count) > 1000, "A", ifelse(
  exp(ndat$count) < 1000 & exp(ndat$count) > 100, "B", ifelse(
    exp(ndat$count) < 100 & exp(ndat$count) > 20, "C", "D"
  )
))

pg3sa <- pg3 %>% dplyr::select(lat,lon, predS, x, y)
pg3ca <- pg3c %>% dplyr::select(lat, lon, predC, x, y)

join_dat <- left_join(pg3ca, pg3sa)
join_dat <- join_dat %>% mutate(
  dummy1000 = ifelse(exp(predS) >1000 & exp(predC) > 1000, "Both > 1000", ifelse(
    exp(predS) < 1000  & exp(predC) > 1000, "PCIN Only > 1000",ifelse(
      exp(predS) > 1000  & exp(predC) < 1000,"PSHEN Only > 1000", "Both < 1000"))),
  
  dummy500 = ifelse(exp(predS) >500 & exp(predC) > 500, "Both > 500", ifelse(
    exp(predS) < 500  & exp(predC) > 500, "PCIN Only > 500",ifelse(
      exp(predS) > 500  & exp(predC) < 500,"PSHEN Only > 500", "Both < 500"))), 
  
  dummy250 = ifelse(exp(predS) >250 & exp(predC) > 250, "Both > 250", ifelse(
    exp(predS) < 250  & exp(predC) > 250, "PCIN Only > 250",ifelse(
      exp(predS) > 250  & exp(predC) < 250,"PSHEN Only > 250", "Both < 250"))), 
  
  dummy1 = ifelse(exp(predS) > 1 & exp(predC) >1, "Both Present", ifelse(
      exp(predS) < 1 & exp(predC) >1, "PCIN Only", ifelse(
        exp(predS) > 1 & exp(predC) < 1, "PSHEN Only", "Both Absent"
      )
    )) )
        
 
ndat$dummy <- interaction(ndat$cd, ndat$spec)

ggplot() + geom_tile(data = ndat, aes(x=, y=y, color= spec, fill = dummy, alpha = spec), key_glyph = "dotplot") +
  scale_color_manual("species", values = c("deepskyblue4", "black") ) + 
  scale_alpha_manual("species", values = c(0.5, 1)) +
  scale_fill_manual("Predicted counts",
                    labels = c(">1000 P. cinereus ",
                                "> 100 P. cinereus",
                                 "> 20 P.cinereus", 
                                "< 20 P. cinereus",
                               ">1000 P. shenandoah ",
                               "> 100 P. shenandoah",
                               "> 20 P.shenandoah", 
                               "< 20 P. shenandoah"),
                    values = c("blue4", "deepskyblue4", "deepskyblue2", "cadetblue1",
                               "coral4", "coral3", "coral2", "coral")) +
  xlab("Longitude") + ylab("Latitude")+ 
  guides(color = guide_legend(override.aes = list(fill = c("deepskyblue4", "black")) )) 
#plum2




(highct <- ggplot() + geom_tile(data = join_dat, aes(x=x, y=y,  fill = dummy1000, ), key_glyph = "dotplot") +
  xlab("Longitude") + ylab("Latitude") + 
  scale_fill_manual("Expected \n counts ", values = c("slategray3", "lightcoral", "aquamarine" )) +
    ggtitle("A.") +
  theme(legend.position = c(0.25,0.83), legend.background = element_blank() ))

(ct500 <- ggplot() + geom_tile(data = join_dat, aes(x=x, y=y,  fill = dummy500, ), key_glyph = "dotplot") +
  xlab("Longitude") + ylab("Latitude") + 
  scale_fill_manual("Expected counts ", values = c("slategray3", "darkorchid1", "lightcoral", "aquamarine" )) +
  ggtitle("B.") +
  theme(legend.position = c(0.25,0.83), legend.background = element_blank()))


(ct250 <- ggplot() + geom_tile(data = join_dat, aes(x=x, y=y,  fill = dummy250, ), key_glyph = "dotplot") +
  xlab("Longitude") + ylab("Latitude") + 
  scale_fill_manual("Expected Counts ", values = c("slategray3", "darkorchid1", "lightcoral", "aquamarine" )) +
  ggtitle("C.")+
    theme(legend.position = c(0.25,0.83), legend.background = element_blank()))

(pres <- ggplot() + geom_tile(data = join_dat, aes(x=x, y=y,  fill = dummy1, ), key_glyph = "dotplot") +
  xlab("Longitude") + ylab("Latitude") + 
  scale_fill_manual("Expected counts ", values = c("slategray3", "darkorchid1", "lightcoral", "aquamarine")) + ggtitle("D.")+
    theme(legend.position = c(0.25,0.83), legend.background = element_blank()))

grid.arrange(highct, ct500, ct250, pres, nrow = 2)

###########################   DATA SUMMARIES   #########################

## how big each status is:  ## divide by 100 = sq 
# add both present + species of interest for total area with predicted mean >1
join_dat %>% group_by(dummy1) %>% summarize(ct = n())

## how much area where LI is greater than 1
pg2 %>% mutate(new_stat = ifelse(exp(predS) >=1 & LI > 1, "pres", "abs")) %>%
  group_by(new_stat) %>% summarize(ct = n())

pg2C %>% mutate(new_stat = ifelse(exp(predC) >=1 & LI > 1, "pres", "abs")) %>%
  group_by(new_stat) %>% summarize(ct = n())
km





###########################    Supplemental Covariate Figures  #########################
fin_mod_out <- readRDS("Pshen_mod_outzIP_aug9.RDS")
fin_mod_outC <- readRDS("cin_mod_out_Zip_aug9.RDS")

shen_gam25 <- fin_mod_out[[1]]
p_objS25 <- plot(shen_gam25, residuals = TRUE)
p_obj125 <- p_objS25[[1]] # just one smooth so select the first component-- Elevation
sm_df25 <- as.data.frame(p_obj125[c("x", "se", "fit")])
data_df25 <- as.data.frame(p_obj125[c("raw", "p.resid")])


sm25 <- sm_df25 
sm25 <- sm25 %>% 
  mutate(unscaled = plyr::mapvalues(round(x,2), from = round(elv$elev,2), to = elv$elev1))
sm25 <- sm25 %>% filter(unscaled >200) ## 3 values not matching



(elev25 <- ggplot(sm25, aes(x = unscaled, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Elevation (m)", y = "Effect Size") +ggtitle("A. P. Shenandoah") + ylim(c(-22,2)) )


p_obj225 <- p_objS25[[2]] 
sm_df25 <- as.data.frame(p_obj225[c("x", "se", "fit")])
data_df25 <- as.data.frame(p_obj225[c("raw", "p.resid")])


(imi25 <- ggplot(sm_df25, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized IMI", y = "Effect Size") +ggtitle("B.") + ylim(c(-15,5)) )

p_obj325 <- p_objS25[[3]]
sm_df25 <- as.data.frame(p_obj325[c("x", "se", "fit")])
data_df25 <- as.data.frame(p_obj225[c("raw", "p.resid")])


(hli25 <- ggplot(sm_df25, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized HLI", y = "Effect Size") +ggtitle("C.") + ylim(c(-5,1.5) ))


shen_gam50 <- fin_mod_out[[2]]
p_objS50 <- plot(shen_gam50, residuals = TRUE)
p_obj150 <- p_objS50[[1]] # just one smooth so select the first component-- Elevation
sm_df50 <- as.data.frame(p_obj150[c("x", "se", "fit")])
data_df50 <- as.data.frame(p_obj150[c("raw", "p.resid")])

sm50 <- sm_df50
sm50 <- sm50 %>% 
  mutate(unscaled = plyr::mapvalues(round(x,2), from = round(elv$elev,2), to = elv$elev1))
sm50 <- sm50 %>% filter(unscaled >200) ## 3 values not matching


(elev50 <- ggplot(sm50, aes(x = unscaled, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Elevation (m)", y = "Effect Size") +ggtitle("A. P. Shenandoah") + ylim(c(-22,2)) )


p_obj250 <- p_objS50[[2]] # just one smooth so select the first component-- Elevation
sm_df50 <- as.data.frame(p_obj250[c("x", "se", "fit")])
data_df50 <- as.data.frame(p_obj250[c("raw", "p.resid")])


(imi50 <- ggplot(sm_df50, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized IMI", y = "Effect Size") +ggtitle("B.") )

p_obj350 <- p_objS50[[3]] # just one smooth so select the first component-- Elevation
sm_df50 <- as.data.frame(p_obj350[c("x", "se", "fit")])
data_df50 <- as.data.frame(p_obj250[c("raw", "p.resid")])


(hli50 <- ggplot(sm_df50, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized HLI", y = "Effect Size") +ggtitle("C.")+ ylim(c(-2.5,1)) )


### cinereus
cin_gam25 <- fin_mod_outC[[1]]
p_objC25 <- plot(cin_gam25, residuals = TRUE)
p_objC125 <- p_objC25[[1]] # just one smooth so select the first component-- Elevation
sm_dfC25 <- as.data.frame(p_objC125[c("x", "se", "fit")])
data_dfC25 <- as.data.frame(p_objC125[c("raw", "p.resid")])

sm25c <- sm_dfC25
sm25c <- sm25c%>% 
  mutate(unscaled = plyr::mapvalues(round(x,2), from = round(elv$elev,2), to = elv$elev1))
sm25c <- sm25c %>% filter(unscaled >200) ## 3 values not matching 

(elevC25 <- ggplot(sm25c, aes(x = unscaled, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Elevation (m)", y = "Effect Size") +ggtitle("D. P. cinereus") + ylim(c(-22,2)) )

p_obj2C25 <- p_objC25[[2]] # just one smooth so select the first component-- Elevation
sm_dfC25 <- as.data.frame(p_obj2C25[c("x", "se", "fit")])
data_dfC25 <- as.data.frame(p_obj2C25[c("raw", "p.resid")])


(hliC25 <- ggplot(sm_dfC25, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized HLI", y = "Effect Size") +ggtitle("F.") + ylim(c(-5,1.5) ))



p_obj2C25 <- p_objC25[[3]] 
sm_dfC25 <- as.data.frame(p_obj2C25[c("x", "se", "fit")])
data_dfC25 <- as.data.frame(p_obj2C25[c("raw", "p.resid")])


(imiC25 <- ggplot(sm_dfC25, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized IMI", y = "Effect Size") +ggtitle("E.") + ylim(c(-15,5)) )


cin_gam50 <- fin_mod_outC[[2]]
p_objC50 <- plot(cin_gam50, residuals = TRUE)
p_obj1C50 <- p_objC50[[1]] # just one smooth so select the first component-- Elevation
sm_dfC50 <- as.data.frame(p_obj1C50[c("x", "se", "fit")])
data_dfC50 <- as.data.frame(p_obj1C50[c("raw", "p.resid")])

sm50c <- sm_dfC50
sm50c <- sm50c%>% 
  mutate(unscaled = plyr::mapvalues(round(x,2), from = round(elv$elev,2), to = elv$elev1))
sm50c <- sm50c %>% filter(unscaled >200) ## 3 values not matching 


(elevC50 <- ggplot(sm50c, aes(x = unscaled, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Elevation (m)", y = "Effect Size") +ggtitle("D. P. cinereus") + ylim(c(-22,2) ))

p_obj2C50 <- p_objC50[[2]] # just one smooth so select the first component-- Elevation
sm_dfC50 <- as.data.frame(p_obj2C50[c("x", "se", "fit")])
data_dfC50 <- as.data.frame(p_obj2C50[c("raw", "p.resid")])


(hliC50 <- ggplot(sm_dfC50, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized HLI", y = "Effect Size") +ggtitle("F.") + ylim(c(-2.5,1)) )



p_obj2C50 <- p_objC50[[3]] # just one smooth so select the first component-- Elevation
sm_dfC50 <- as.data.frame(p_obj2C50[c("x", "se", "fit")])
data_dfC50 <- as.data.frame(p_obj2C50[c("raw", "p.resid")])


(imiC50 <- ggplot(sm_dfC50, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized IMI", y = "Effect Size") +ggtitle("E.") )



### 25
grid.arrange(elev25,imi25, hli25, elevC25, imiC25, hliC25, nrow = 2)

### 50 
grid.arrange(elev50,imi50, hli50, elevC50, imiC50, hliC50, nrow = 2)



##### ppt and temp covariate figures ---- 
## 25 
pl4 <- plot(shen_gam25, select = 4, scale = 0, 
            ylab = "Temp", shade = TRUE, shade.col = "gray80", lwd = 2)
ps4 <- pl4[[4]]


ps4q <- as.data.frame(ps4[c("x",  "y")])


tempo <- expand.grid(ps4q$x, ps4q$y) 

temp2 <- cbind(tempo, ps4$fit) %>% rename(x= Var1, y = Var2, fit = "ps4$fit" ) %>%
  filter(!is.na(fit))

mycol <- c( "blueviolet","blue3" ,"turquoise3" ,"aquamarine1",
            "gray73", "lightyellow", "yellow" , "gold1", "orange" , "orangered")

(stp <- ggplot(temp2, aes(x=x, y=y, z = fit)) + 
    geom_contour_filled(breaks =  c(-1,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1,1.2)) + 
    xlab("3-day precipitation") + ylab("Temperature") + ggtitle("A. P. shenandoah") +
    scale_fill_manual("Effect Size", values = mycol))
(stp <- stp + theme(legend.position = (c(0.7,0.25)), legend.direction = "horizontal", legend.background = element_blank()) + 
    annotate("text", x= 27, y = 7, label= "Lower count") +
    annotate("text", x= 37, y = 7, label = "Higher count"))



cl4 <- plot(cin_gam25, select = 4, scale = 0, 
            ylab = "Temp", shade = TRUE, shade.col = "gray80", lwd = 2)
cs4 <- cl4[[4]]


cs4q <- as.data.frame(cs4[c("x",  "y")])


tempoc <- expand.grid(cs4q$x, cs4q$y) 

temp2c <- cbind(tempoc, cs4$fit) %>% rename(x= Var1, y = Var2, fit = "cs4$fit" ) %>%
  filter(!is.na(fit))



(ctp <- ggplot(temp2c, aes(x=x, y=y, z = fit)) + geom_contour_filled(breaks = c(-1,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1,1.2)) + 
    xlab("3-day precipitation") + ylab("Temperature") + ggtitle("B. P. cinereus") +
    scale_fill_manual("Effect Size", values = mycol) + theme(legend.position = "none"))


gridExtra::grid.arrange(stp, ctp)

##  50 
pl4 <- plot(shen_gam50, select = 4, scale = 0, 
            ylab = "Temp", shade = TRUE, shade.col = "gray80", lwd = 2)
ps4 <- pl4[[4]]


ps4q <- as.data.frame(ps4[c("x",  "y")])


tempo <- expand.grid(ps4q$x, ps4q$y) 

temp2 <- cbind(tempo, ps4$fit) %>% rename(x= Var1, y = Var2, fit = "ps4$fit" ) %>%
  filter(!is.na(fit))



(stp <- ggplot(temp2, aes(x=x, y=y, z = fit)) + 
    geom_contour_filled(breaks = c(-1,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6)) + 
    xlab("3-day precipitation") + ylab("Temperature") + ggtitle("A. P. shenandoah") +
    scale_fill_manual("Effect Size", values = mycol))
(stp <- stp + theme(legend.position = (c(0.7,0.25)), legend.direction = "horizontal", legend.background = element_blank()) + 
    annotate("text", x= 27, y = 7, label= "Lower count") +
    annotate("text", x= 37, y = 7, label = "Higher count"))



cl4 <- plot(cin_gam50, select = 4, scale = 0, 
            ylab = "Temp", shade = TRUE, shade.col = "gray80", lwd = 2)
cs4 <- cl4[[4]]


cs4q <- as.data.frame(cs4[c("x",  "y")])


tempoc <- expand.grid(cs4q$x, cs4q$y) 

temp2c <- cbind(tempoc, cs4$fit) %>% rename(x= Var1, y = Var2, fit = "cs4$fit" ) %>%
  filter(!is.na(fit))



(ctp <- ggplot(temp2c, aes(x=x, y=y, z = fit)) + geom_contour_filled(breaks = c( -1,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6)) + 
    xlab("3-day precipitation") + ylab("Temperature") + ggtitle("B. P. cinereus") +
    scale_fill_manual("Effect Size", values = mycol)
  + theme(legend.position = "none")
   )


gridExtra::grid.arrange(stp, ctp)




###########################    Supplemental Map #########################
fin_mod_out <- readRDS("Pshen_mod_outzIP_aug9.RDS")
fin_mod_outC <- readRDS("cin_mod_out_Zip_aug9.RDS")

cov$area <- log(10000)

#cov <- subset(cov, !is.na(IMI))

cov <- cov %>% mutate(asp_rad = pi*aspect/180) %>%
  mutate(north = cos(asp_rad), east = sin(asp_rad)) %>%
  mutate(std_HLI = scale(HLI)[,1], std_IMI = scale(IMI)[,1])

cov <- cbind(cov, pred_weather[,c("ppt3", "temp10")])

cov <- cov %>%
  rename("ppt" = "ppt3")

pred_gamS25 <- predict(fin_mod_out[[1]], data.frame(cov, area = rep(10000,nrow(cov))), se.fit = TRUE)



pg225 <- pg <- data.frame(predS = pred_gamS25$fit,seS = pred_gamS25$se.fit, 
                          lat = cov$lat, lon = cov$lon,x = cov$UTMx,
                          y = cov$UTMy, elev <- cov$elev, aspect = cov$aspect, 
                          IMI = cov$IMI, std_HLI = cov$std_HLI, left = cov$left_, right = cov$right_,
                          top = cov$top_, bottom = cov$bottom, predS1 = pred_gamS25$fit)



pg225$LI <- exp(pg225$predS- 2*pg225$seS)

pg225$status <- ifelse(exp(pg225$predS) >  1, "mean greater than 1", "B")

#### P. cinereus 
pred_gamC25 <- predict(fin_mod_outC[[1]], data.frame(cov, area = rep(10000,nrow(cov))), se.fit = TRUE)



pg2C25 <- pg <- data.frame(predC = pred_gamC25$fit,seC = pred_gamC25$se.fit, 
                           lat = cov$lat, lon = cov$lon,x = cov$UTMx,
                           y = cov$UTMy, elev <- cov$elev, aspect = cov$aspect, 
                           IMI = cov$IMI, std_HLI = cov$std_HLI, left = cov$left_, right = cov$right_,
                           top = cov$top_, bottom = cov$bottom, predC1 = pred_gamC25$fit)


pg2C25$status <- ifelse(exp(pg2C25$predC) >  1, "mean greater than 1", "B")
pg2C25$LI <- exp(pg2C25$predC- 2*pg2C25$seC)
#pg2C25$UI <- exp(pg2C25$predC+ 2*pg2C25$seC)
#pg2C25$status <- "Uncertain"
#pg2C25$status[pg2C25$LI > 1 | exp(pg2C25$predC) > 20] <- "Highly Likely"

#pg2C25$status[exp(pg2C25$predC) <= 1] <- "Unlikely"
#pg3C25 <- subset(pg2C25, exp(predC) > 0.1)


## how big each status is
pg225 %>% group_by(status) %>% summarize(ct = n())
pg225 %>% mutate(new_stat = ifelse(exp(predS1) >=1 & LI > 1, "pres", "abs")) %>%
  group_by(new_stat) %>% summarize(ct = n())
pg2C25 %>% group_by(status) %>% summarize(ct = n())
pg2C25 %>% mutate(new_stat = ifelse(exp(predC) >=1 & LI > 1, "pres", "abs")) %>%
  group_by(new_stat) %>% summarize(ct = n())



## 50

pred_gamS50 <- predict(fin_mod_out[[2]], data.frame(cov, area = rep(10000,nrow(cov))), se.fit = TRUE)



pg250 <- pg <- data.frame(predS = pred_gamS50$fit,seS = pred_gamS50$se.fit, 
                          lat = cov$lat, lon = cov$lon,x = cov$UTMx,
                          y = cov$UTMy, elev <- cov$elev, aspect = cov$aspect, 
                          IMI = cov$IMI, std_HLI = cov$std_HLI, left = cov$left_, right = cov$right_,
                          top = cov$top_, bottom = cov$bottom, predS1 = pred_gamS50$fit)





pg250$status <- ifelse(exp(pg250$predS) >1, "mean greater than 1", "B")
pg250$LI <- exp(pg250$predS- 2*pg250$seS)


#### P. cinereus 
pred_gamC50 <- predict(fin_mod_outC[[2]], data.frame(cov, area = rep(10000,nrow(cov))), se.fit = TRUE)



pg2C50 <- pg <- data.frame(predC = pred_gamC50$fit,seC = pred_gamC50$se.fit, 
                           lat = cov$lat, lon = cov$lon,x = cov$UTMx,
                           y = cov$UTMy, elev <- cov$elev, aspect = cov$aspect, 
                           IMI = cov$IMI, std_HLI = cov$std_HLI, left = cov$left_, right = cov$right_,
                           top = cov$top_, bottom = cov$bottom, predC1 = pred_gamC50$fit)




pg2C50$status <- ifelse(exp(pg2C50$predC) >1, "mean greater than 1", "B")
pg2C50$LI <- exp(pg2C50$predC-2*pg2C50$seC)

## how big each status is
pg250 %>% group_by(status) %>% summarize(ct = n())
pg250 %>% mutate(new_stat = ifelse(exp(predS1) >=1 & LI > 1, "pres", "abs")) %>%
  group_by(new_stat) %>% summarize(ct = n())
pg2C50 %>% group_by(status) %>% summarize(ct = n())
pg2C50 %>% mutate(new_stat = ifelse(exp(predC) >=1 & LI > 1, "pres", "abs")) %>%
  group_by(new_stat) %>% summarize(ct = n())


## 200

pred_gamS200 <- predict(fin_mod_out[[3]], data.frame(cov, area = rep(10000,nrow(cov))), se.fit = TRUE)



pg2200 <- pg <- data.frame(predS = pred_gamS200$fit,seS = pred_gamS200$se.fit, 
                           lat = cov$lat, lon = cov$lon,x = cov$UTMx,
                           y = cov$UTMy, elev <- cov$elev, aspect = cov$aspect, 
                           IMI = cov$IMI, std_HLI = cov$std_HLI, left = cov$left_, right = cov$right_,
                           top = cov$top_, bottom = cov$bottom, predS1 = pred_gamS200$fit)



pg2200$status <- ifelse(exp(pg2200$predS) >1, "mean greater than 1", "B")


#### P. cinereus 
pred_gamC200 <- predict(fin_mod_outC[[3]], data.frame(cov, area = rep(10000,nrow(cov))), se.fit = TRUE)



pg2C200 <- pg <- data.frame(predC = pred_gamC200$fit,seC = pred_gamC200$se.fit, 
                            lat = cov$lat, lon = cov$lon,x = cov$UTMx,
                            y = cov$UTMy, elev <- cov$elev, aspect = cov$aspect, 
                            IMI = cov$IMI, std_HLI = cov$std_HLI, left = cov$left_, right = cov$right_,
                            top = cov$top_, bottom = cov$bottom, predC1 = pred_gamC200$fit)



pg2C200$status <- ifelse(exp(pg2C200$predC) >1, "mean greater than 1", "B")

## how big each status is
pg2200 %>% group_by(status) %>% summarize(ct = n())
pg2C200 %>% group_by(status) %>% summarize(ct = n())
