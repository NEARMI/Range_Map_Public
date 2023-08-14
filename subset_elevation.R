elevs <- sh.plot %>% 
  ungroup() %>%
  distinct(elev) %>%
  arrange(elev) %>%
  mutate(id_col = seq(n())) %>%
  filter((id_col / 3) %% 1 == 0)


el_dat <- sh.plot %>% filter(elev %in% elevs$elev)  ## 127 sites (unique(el_dat$ID)) ~ 1/3

wiggles <- c(15,35,75)

## takes ~10 min
system.time(
  for(t in 1: length(scen)){
    print("tt")
    for(j in 1:length(scen)){
      print("jj")
      ## leave out one- ninth of the area 
      ## summarize each site to mean temp and mean ppt -- other values are the same by line
      ## this is because I want to predict by site not by visit
      sh.plot1 <- el_dat %>%
        filter(x >= min(el_dat$x) + (j-1)*chx & x < min(el_dat$x) + j*chx  &
                 y >= min(el_dat$y) + (t-1)*chy & y < min(el_dat$y) + t*chy) %>%
        group_by(ID) %>%
        summarize(temp10 = mean(temp10, na.rm =T), ppt = mean(ppt, na.rm =T), std_HLI = mean(std_HLI), elev = mean(elev),
                  x = mean(x), y = mean(y), north = mean(north), east = mean(east), std_IMI= mean(std_IMI),
                  shc = sum(shencount), cic = sum(cincount)) %>%
        ungroup()  %>% mutate(stt = ifelse(shc == 0, 0, 1), ctt = ifelse(cic == 0, 0, 1))
      
      
      sh.plot2 <- el_dat %>% filter(ID %notin% sh.plot1$ID)
      
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


outpred <- out %>% 
  filter(rm_sites > 0, rm_sites != 35) %>% 
  group_by(wiggles, rm_sites) %>% 
  summarize(tt = mean(rmse, na.rm=T), cc = mean(cov_prop, na.rm=T))

outpred %>% group_by(wiggles) %>% summarize(rm = mean(tt,na.rm=T), cm = mean(cc, na.rm= T))


shen_gam_elev <-  gam(shencount ~-1+ s(elev, k = 6) 
                      + s(std_IMI, k=6) 
                      + s(std_HLI, k = 6) 
                      + s(ppt, temp10, bs='ds', k =14, m=2)
                      + s(x, y, bs='ds', k =75, m=2) 
                      + s(ID,bs = "re")
                      ,
                      data=el_dat, family = "ziP",offset = log(area/10000), method = 'REML',
                      drop.unused.levels = FALSE)


shen_gam_elev15 <-  gam(shencount ~-1+ s(elev, k = 6) 
                        + s(std_IMI, k=6) 
                        + s(std_HLI, k = 6) 
                        + s(ppt, temp10, bs='ds', k =14, m=2)
                        + s(x, y, bs='ds', k =15, m=2) 
                        + s(ID,bs = "re")
                        ,
                        data=el_dat, family = "ziP",offset = log(area/10000), method = 'REML',
                        drop.unused.levels = FALSE)


shen_gam_elev35 <-  gam(shencount ~-1+ s(elev, k = 6) 
                        + s(std_IMI, k=6) 
                        + s(std_HLI, k = 6) 
                        + s(ppt, temp10, bs='ds', k =14, m=2)
                        + s(x, y, bs='ds', k =35, m=2) 
                        + s(ID,bs = "re")
                        ,
                        data=el_dat, family = "ziP",offset = log(area/10000), method = 'REML',
                        drop.unused.levels = FALSE)


## covariates with reduced data
par(mfrow=c(3,3))
p_objS <- plot(shen_gam_elev, residuals = TRUE)
p_obj1 <- p_objS[[1]] # just one smooth so select the first component-- Elevation
sm_df <- as.data.frame(p_obj1[c("x", "se", "fit")])
data_df <- as.data.frame(p_obj1[c("raw", "p.resid")])


sm25elev <- sm_df
sm25elev <- sm25elev%>% 
  mutate(unscaled = plyr::mapvalues(round(x,2), from = round(elv$elev,2), to = elv$elev1))
sm25elev <- sm25elev %>% filter(unscaled >200)

(elevS <- ggplot(sm25elev, aes(x = unscaled, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = " ", y = " ") +ggtitle("B. Subset by elevation")) #+ ylim(-40,20) )

p_obj2 <- p_objS[[2]] # just one smooth so select the second component-- IMI
sm_df <- as.data.frame(p_obj2[c("x", "se", "fit")])
data_df <- as.data.frame(p_obj2[c("raw", "p.resid")])


(imi<- ggplot(sm_df, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = " ", y = " ") +ggtitle("B. Subset by elevation") + ylim(-50,15) )


p_obj2 <- p_objS[[3]] 
sm_df <- as.data.frame(p_obj2[c("x", "se", "fit")])
data_df <- as.data.frame(p_obj2[c("raw", "p.resid")])


(hli<- ggplot(sm_df, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = " ", y = " ") +ggtitle("B. Subset by elevation") + ylim(-15,6))



plot(shen_gam_elev, select = 4, scale = 0, 
     ylab = "Temperature", xlab = "Precipitation", 
     shade = TRUE, shade.col = "gray80", lwd = 2 )

elevppt <- plot(shen_gam_elev, select = 4, scale = 0, 
            ylab = "Temp", shade = TRUE, shade.col = "gray80", lwd = 2)
elev4 <- elevppt[[4]]


elev4q <- as.data.frame(elev4[c("x",  "y")])


tempoelev <- expand.grid(elev4q$x, elev4q$y) 

temp2elev <- cbind(tempoelev, elev4$fit) %>% rename(x= Var1, y = Var2, fit = "elev4$fit" ) %>%
  filter(!is.na(fit))

#mycol <- c( "blueviolet","blue3" ,"turquoise3", "springgreen4" ,  "springgreen2" ,"aquamarine1",
# "gray73", "goldenrod4" , "gold3", "goldenrod2" , "orange" ,"orangered")

mycol <- c( "blueviolet","blue3" ,"turquoise3" ,"aquamarine1",
            "gray73", "yellow" , "gold1", "orange" , "orangered")

(elevpt <- ggplot(temp2elev, aes(x=x, y=y, z = fit)) + 
    geom_contour_filled(breaks = c(-1, -0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1,1.2)) + 
    xlab("3-day precipitation") + ylab("Temperature") + ggtitle("B. Subset by elevation") +
    scale_fill_manual("Effect Size", values = mycol) + theme(legend.position = "none"))





##### map with reduced data


pred_gamS_elev <- predict(shen_gam_elev, data.frame(cov, area = rep(10000,nrow(cov))), se.fit = TRUE)

y = pred_gamS_elev$fit
x = pred_gamS_elev$se.fit
ZZ = "B. Subset by elevation"

elm <- makemap(y=y, x=x, ZZ=ZZ)
#elmerr <- makemaperror(y=y, x=x, ZZ=ZZ)

pg2elev <- pg <- data.frame(predS = pred_gamS_elev$fit,seS = pred_gamS_elev$se.fit, 
                        lat = cov$lat, lon = cov$lon,x = cov$UTMx,
                        y = cov$UTMy, elev <- cov$elev, aspect = cov$aspect, 
                        IMI = cov$IMI, std_HLI = cov$std_HLI, left = cov$left_, right = cov$right_,
                        top = cov$top_, bottom = cov$bottom, predS1 =pred_gamS_elev$fit)





pg2elev$status <- ifelse(exp(pg2elev$predS) > 1, "greater than 1","B")
#pg3 <- subset(pg2, plogis(predS) > 0.1)
pg2elev$LI <- exp(pg2elev$predS1-2*pg2elev$seS)

## how big each status is

pg2elev %>% mutate(new_stat = ifelse(exp(predS1) >=1 & LI > 1, "pres", "abs")) %>%
  group_by(new_stat) %>% summarize(ct = n())

pg2elev %>% group_by(status) %>% summarize(ct = n())
