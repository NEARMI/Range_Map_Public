
ddat <- read.csv("non_geo_dist.csv")

wiggles <- c(15,35,75)

## takes ~30 min out put saved in  "subset_dist_road_wiggles.RDS"
system.time(
  for(t in 1: length(scen)){
    print("tt")
    for(j in 1:length(scen)){
      print("jj")
      ## leave out one- ninth of the area 
      ## summarize each site to mean temp and mean ppt -- other values are the same by line
      ## this is because I want to predict by site not by visit
      sh.plot1 <- ddat %>%
        filter(x >= min(ddat$x) + (j-1)*chx & x < min(ddat$x) + j*chx  &
                 y >= min(ddat$y) + (t-1)*chy & y < min(ddat$y) + t*chy) %>%
        group_by(ID) %>%
        summarize(temp10 = mean(temp10, na.rm =T), ppt = mean(ppt, na.rm =T), std_HLI = mean(std_HLI), elev = mean(elev),
                  x = mean(x), y = mean(y), std_IMI= mean(std_IMI),
                  shc = sum(shencount), cic = sum(cincount)) %>%
        ungroup()  %>% mutate(stt = ifelse(shc == 0, 0, 1), ctt = ifelse(cic == 0, 0, 1))
      
      
      sh.plot2 <- ddat %>% filter(ID %notin% sh.plot1$ID)
      
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
        
      } ## close i
    } ## close j
  } ## close t
  
) 

out <- as.data.frame(out)
names(out) <- c("rmse","cov_prop","wiggles", "j", "t", "rm_sites")
#saveRDS(out, "subset_dist_road_wiggles.RDS")

outpred <- out %>% 
  filter(rm_sites > 0, rm_sites != 64) %>%  ## extremely high rmse when 64 sites removed don't know why- doesn't really change overall pattern
  group_by(wiggles, rm_sites) %>% 
  summarize(tt = mean(rmse, na.rm=T), cc = mean(cov_prop, na.rm=T))

outpred %>% group_by(wiggles) %>% summarize(rm = mean(tt,na.rm=T), cm = mean(cc, na.rm= T))


shen_gam_dis <-  gam(shencount ~-1+ s(elev, k = 6) 
                      + s(std_IMI, k=6) 
                      + s(std_HLI, k = 6) 
                      + s(ppt, temp10, bs='ds', k =14, m=2)
                      + s(x, y, bs='ds', k =15, m=2) 
                      + s(ID,bs = "re")
                      ,
                      data=ddat, family = "ziP",offset = log(area/10000), method = 'REML',
                      drop.unused.levels = FALSE)


shen_gam_dis35 <-  gam(shencount ~-1+ s(elev, k = 6) 
                        + s(std_IMI, k=6) 
                        + s(std_HLI, k = 6) 
                        + s(ppt, temp10, bs='ds', k =14, m=2)
                        + s(x, y, bs='ds', k =35, m=2) 
                        + s(ID,bs = "re")
                        ,
                        data=ddat, family = "ziP",offset = log(area/10000), method = 'REML',
                        drop.unused.levels = FALSE)


shen_gam_dis75 <-  gam(shencount ~-1+ s(elev, k = 6) 
                        + s(std_IMI, k=6) 
                        + s(std_HLI, k = 6) 
                        + s(ppt, temp10, bs='ds', k =14, m=2)
                        + s(x, y, bs='ds', k =75, m=2) 
                        + s(ID,bs = "re")
                        ,
                        data=ddat, family = "ziP",offset = log(area/10000), method = 'REML',
                        drop.unused.levels = FALSE)


## covariates with reduced data
par(mfrow=c(3,3))
p_objS <- plot(shen_gam_dis75, residuals = TRUE)
p_obj1 <- p_objS[[1]] # just one smooth so select the first component-- Elevation
sm_df <- as.data.frame(p_obj1[c("x", "se", "fit")])
data_df <- as.data.frame(p_obj1[c("raw", "p.resid")])



sm25dist <- sm_df
sm25dist <- sm25dist%>% 
  mutate(unscaled = plyr::mapvalues(round(x,2), from = round(elv$elev,2), to = elv$elev1))
sm25dist <- sm25dist %>% filter(unscaled >200)

(elevD <- ggplot(sm25dist, aes(x = unscaled, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Elevation (m)", y = " ") +ggtitle("E. Subset by Distance")) #+ ylim(-40,20) )

p_obj2 <- p_objS[[2]] # just one smooth so select the second component-- IMI
sm_df <- as.data.frame(p_obj2[c("x", "se", "fit")])
data_df <- as.data.frame(p_obj2[c("raw", "p.resid")])


(imiD<- ggplot(sm_df, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized IMI", y = "") +ggtitle("E. Subset by Distance") + ylim(-50,15) )


p_obj2 <- p_objS[[3]] 
sm_df <- as.data.frame(p_obj2[c("x", "se", "fit")])
data_df <- as.data.frame(p_obj2[c("raw", "p.resid")])


(hliD<- ggplot(sm_df, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized HLI", y = " ") +ggtitle("E. Subset by Distance") + ylim(-15,6))



roadppt <- plot(shen_gam_dis75, select = 4, scale = 0, 
                ylab = "Temp", shade = TRUE, shade.col = "gray80", lwd = 2)
road4 <- roadppt[[4]]


road4q <- as.data.frame(road4[c("x",  "y")])


temporoad <- expand.grid(road4q$x, road4q$y) 

temp2road <- cbind(temporoad, road4$fit) %>% rename(x= Var1, y = Var2, fit = "road4$fit" ) %>%
  filter(!is.na(fit))

mycol <- c( "blueviolet","blue3" ,"turquoise3" ,"aquamarine1",
            "gray73", "yellow" , "gold1", "orange" , "orangered")



(roadpt <- ggplot(temp2road, aes(x=x, y=y, z = fit)) + 
    geom_contour_filled(breaks = c(-1, -0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1)) + 
    xlab("3-day precipitation") + ylab("Temperature") + ggtitle("E. Subset by Distance") +
    scale_fill_manual("Effect Size", values = mycol) + theme(legend.position = "none"))



##### map with reduced data


pred_gamS_dis <- predict(shen_gam_dis75, (cov %>% mutate(area = 10000)), se.fit = TRUE)


y = pred_gamS_dis$fit
x = pred_gamS_dis$se.fit
ZZ = "E. Subset by distance"

elmD <- makemap(y=y, x=x, ZZ=ZZ)
#elmerrD <- makemaperror(y=y, x=x, ZZ=ZZ)

pg2dis <- pg <- data.frame(predS = pred_gamS_dis$fit,seS = pred_gamS_dis$se.fit, 
                            lat = cov$lat, lon = cov$lon,x = cov$UTMx,
                            y = cov$UTMy, elev <- cov$elev, aspect = cov$aspect, 
                            IMI = cov$IMI, std_HLI = cov$std_HLI, left = cov$left_, right = cov$right_,
                            top = cov$top_, bottom = cov$bottom)



pg2dis$status <- ifelse(exp(pg2dis$predS) >1 , "greater than 1", "B")

#pg3 <- subset(pg2ehi, plogis(predS) > 0.1)


pg2dis %>% group_by(status) %>% summarize(ct = n())


pg2dis$LI <- exp(pg2dis$predS-2*pg2dis$seS)

## how big each status is

pg2dis %>% mutate(new_stat = ifelse(exp(predS) >=1 & LI > 1, "pres", "abs")) %>%
  group_by(new_stat) %>% summarize(ct = n())
