
sh.plot2 <- sh.plot %>% group_by(UTMx, UTMy, year) %>% slice(1)


wiggles <- c(25, 50, 200)


system.time(
  for(t in 1: length(scen)){
    print("tt")
    for(j in 1:length(scen)){
      print("jj")
      ## leave out one- ninth of the area 
      ## summarize each site to mean temp and mean ppt -- other values are the same by line
      ## this is because I want to predict by site not by visit
      sh.plot1 <- sh.plot2 %>%
        filter(x >= min(sit$x) + (j-1)*chx & x < min(sit$x) + j*chx  &
                 y >= min(sit$y) + (t-1)*chy & y < min(sit$y) + t*chy) %>%
        group_by(ID) %>%
        summarize(temp10 = mean(temp10, na.rm =T), ppt = mean(ppt, na.rm =T), std_HLI = mean(std_HLI, na.rm=T), elev = mean(elev),
                  x = mean(x), y = mean(y), north = mean(north), east = mean(east), std_IMI= mean(std_IMI,na.rm=T),
                  shc = sum(shencount), cic = sum(cincount)) %>%
        ungroup()  %>% mutate(stt = ifelse(shc == 0, 0, 1), ctt = ifelse(cic == 0, 0, 1))
      
      
      sh.plot3 <- sh.plot2 %>% filter(ID %notin% sh.plot1$ID)
      
      for(i in 1:length(wiggles)){
        if(wiggles[i] < length(unique(sh.plot2$ID)) ){
          shen_gam = gam(shencount ~-1+ s(elev, k = 6) 
                         + s(std_IMI, k=6) 
                         + s(std_HLI, k = 6) 
                         + s(ppt, temp10, bs='ds', k =14, m=2)
                         + s(x, y, bs='ds', k =wiggles[i], m=2) 
                         + s(ID,bs = "re")
                         ,
                         data=sh.plot3, family = "ziP",offset = log(area/10000), method = 'REML',
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


outpred <- out %>% 
  filter(rm_sites > 0, rm_sites != 35) %>% 
  group_by(wiggles, rm_sites) %>% 
  summarize(tt = mean(rmse, na.rm=T), cc = mean(cov_prop, na.rm=T))

outpred %>% group_by(wiggles) %>% summarize(rm = mean(tt,na.rm=T), cm = mean(cc, na.rm= T))






shen_gam_vis = gam(shencount~ s(elev, k = 6)
                   +s(std_IMI, k =6)
                   + s(std_HLI, k = 6) 
                   + s(ppt, temp10, bs='ds', k =14, m=2)
                   + s(x, y, bs='ds', k =50, m=2) 
                   + s(ID,bs = "re")
                   ,
                   data= sh.plot2, family = "ziP", method = 'REML', offset = log(area/10000),
                   drop.unused.levels = FALSE)

shen_gam_vis200 = gam(shencount~ s(elev, k = 6)
                      +s(std_IMI, k =6)
                      + s(std_HLI, k = 6) 
                      + s(ppt, temp10, bs='ds', k =14, m=2)
                      + s(x, y, bs='ds', k =200, m=2) 
                      + s(ID,bs = "re")
                      ,
                      data= sh.plot2, family = "ziP", method = 'REML', offset = log(area/10000),
                      drop.unused.levels = FALSE)



shen_gam_vis25 = gam(shencount~ s(elev, k = 6)
                     +s(std_IMI, k =6)
                     + s(std_HLI, k = 6) 
                     + s(ppt, temp10, bs='ds', k =14, m=2)
                     + s(x, y, bs='ds', k =25, m=2) 
                     + s(ID,bs = "re")
                     ,
                     data= sh.plot2, family = "ziP", method = 'REML', offset = log(area/10000),
                     drop.unused.levels = FALSE)


## covariates with reduced data
par(mfrow=c(3,3))
p_objSab <- plot(shen_gam_vis200, residuals = TRUE)
p_obj1ab <- p_objSab[[1]] # just one smooth so select the first component-- Elevation
sm_dfab <- as.data.frame(p_obj1ab[c("x", "se", "fit")])
data_dfab <- as.data.frame(p_obj1ab[c("raw", "p.resid")])


sm25ann <- sm_dfab
sm25ann <- sm25ann %>% 
  mutate(unscaled = plyr::mapvalues(round(x,2), from = round(elv$elev,2), to = elv$elev1))
sm25ann <- sm25ann %>% filter(unscaled >200)

(elev2 <- ggplot(sm25ann, aes(x = unscaled, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Elevation (m)", y = "Effect Size") +
    ggtitle("D. Single vist by year"))  #+ ylim(-40,20) )

p_obj2ab <- p_objSab[[2]] # 
sm_dfab <- as.data.frame(p_obj2ab[c("x", "se", "fit")])
data_dfab <- as.data.frame(p_obj2ab[c("raw", "p.resid")])


(imi2 <- ggplot(sm_dfab, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized IMI ", y = "Effect Size") +ggtitle("D. Single vist by year") + ylim(-50,15) )




p_obj2ab <- p_objSab[[3]] # 
sm_dfab <- as.data.frame(p_obj2ab[c("x", "se", "fit")])
data_dfab <- as.data.frame(p_obj2ab[c("raw", "p.resid")])


(hli2 <- ggplot(sm_dfab, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized HLI ", y = "Effect Size") +ggtitle("D. Single vist by year")  + ylim(-15,6))





plot(shen_gam_vis200, select = 4, scale = 0, 
     ylab = "Temperature", xlab = "Precipitation", 
     shade = TRUE, shade.col = "gray80", lwd = 2 )

annppt <- plot(shen_gam_vis200, select = 4, scale = 0, 
               ylab = "Temp", shade = TRUE, shade.col = "gray80", lwd = 2)
ann4 <- annppt[[4]]


ann4q <- as.data.frame(ann4[c("x",  "y")])


tempoann <- expand.grid(ann4q$x, ann4q$y) 

temp2ann <- cbind(tempoann, ann4$fit) %>% rename(x= Var1, y = Var2, fit = "ann4$fit" ) %>%
  filter(!is.na(fit))



mycol <- c( "blueviolet","blue3" ,"turquoise3" ,"aquamarine1",
            "gray73", "yellow" , "gold1", "orange" , "orangered")

(annpt <- ggplot(temp2ann, aes(x=x, y=y, z = fit)) + 
    geom_contour_filled(breaks = c(-1, -0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1)) + 
    xlab("3-day precipitation") + ylab("Temperature") + ggtitle("D. Single vist by year") +
    scale_fill_manual("Effect Size", values = mycol) + theme(legend.position = "none"))






pred_gamSVis200 <- predict(shen_gam_vis200, (cov %>% mutate(area = 10000)), se.fit = TRUE)



y = pred_gamSVis200$fit
x = pred_gamSVis200$se.fit
ZZ = "D. Single Visit by year"

annual <- makemap(y=y, x=x, ZZ=ZZ)

#annualerr <- makemaperror(y=y, x=x, ZZ=ZZ)

pg2sas <-   data.frame(predS = pred_gamSVis200$fit,seS = pred_gamSVis200$se.fit, 
                       lat = cov$lat, lon = cov$lon,x = cov$UTMx,
                       y = cov$UTMy, elev <- cov$elev, aspect = cov$aspect, 
                       IMI = cov$IMI, std_HLI = cov$std_HLI, left = cov$left_, right = cov$right_,
                       top = cov$top_, bottom = cov$bottom)


pg2sas$status <- ifelse(exp(pg2sas$predS) >1 , "greater than 1", "B")


pg2sas %>% group_by(status) %>% summarize(ct = n())


pg2sas$LI <- exp(pg2sas$predS-2*pg2sas$seS)

## how big each status is

pg2sas %>% mutate(new_stat = ifelse(exp(predS) >=1 & LI > 1, "pres", "abs")) %>%
  group_by(new_stat) %>% summarize(ct = n())
