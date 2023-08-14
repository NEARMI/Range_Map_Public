
## distance formula
dd <- function(x1, x2,y1, y2,z1, z2){
  
  sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
}


reps = 10000

sites <- matrix(data = NA, nrow = 109, ncol = reps)
temp_vec <- vector(mode= "numeric", length = 109)
out <- vector(mode = "numeric", length = reps)
wiggles <- c(15,35,75)
## randomly select 109 sites


#### takes 21 minutes with reps at 10 000 
system.time(
  for(i in 1:ncol(sites)){
    print(i)
    rand_site <- sample(unique(sh.plot$ID), 109, replace = FALSE)
    sites[,i] <- rand_site
    
    temp_mat <- sh.plot %>% filter(ID %in% rand_site)
    
    ## calculate within selection distance (intracluster distance) ## average diameter distance
    
    for(k in 1:nrow(temp_mat)){
      dist_temp <- dd(x1 = temp_mat$elev[k],   x2 = temp_mat$elev[k+1],
                      y1 = temp_mat$HLI[k],  y2 = temp_mat$HLI[k+1],
                      z1 = temp_mat$IMI[k], z2 = temp_mat$IMI[k+1])  
      
      temp_vec[k] <- dist_temp
      
      out[i] <- mean(temp_vec, na.rm=T)
      
    } 
  })

## choose set with largest within cluster distance (largest spread)  
tt <- which(out == max(out))
site_list <- sites[,tt]
saveRDS(site_list, "max_temp_hli_subset_list.RDS")
site_list <- readRDS("max_temp_hli_subset_list.RDS")
sit <- sh.plot %>% filter(ID %in% site_list) ## subset of sites

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




system.time(
  for(t in 1: length(scen)){
    print("tt")
    for(j in 1:length(scen)){
      print("jj")
      ## leave out one- ninth of the area 
      ## summarize each site to mean temp and mean ppt -- other values are the same by line
      ## this is because I want to predict by site not by visit
      sh.plot1 <- sit %>%
        filter(x >= min(sit$x) + (j-1)*chx & x < min(sit$x) + j*chx  &
                 y >= min(sit$y) + (t-1)*chy & y < min(sit$y) + t*chy) %>%
        group_by(ID) %>%
        summarize(temp10 = mean(temp10, na.rm =T), ppt = mean(ppt, na.rm =T), std_HLI = mean(std_HLI), elev = mean(elev),
                  x = mean(x), y = mean(y), north = mean(north), east = mean(east), std_IMI= mean(std_IMI),
                  shc = sum(shencount), cic = sum(cincount)) %>%
        ungroup()  %>% mutate(stt = ifelse(shc == 0, 0, 1), ctt = ifelse(cic == 0, 0, 1))
      
      
      sh.plot2 <- sit %>% filter(ID %notin% sh.plot1$ID)
      
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

saveRDS(out, "elev_hli_imi_leaveoneout.RDS")


shen_gam_eh = gam(shencount ~ s(elev, k = 6)
                  +s(std_IMI, k = 6)
                  + s(std_HLI, k = 6) 
                  + s(ppt, temp10, bs='ds', k =14, m=2)
                  + s(x, y, bs='ds', k =75, m=2) 
                  + s(ID,bs = "re")
                  ,
                  data=sit, family = "ziP", method = 'REML', offset = log(area/10000),
                  drop.unused.levels = FALSE)


shen_gam_eh15 = gam(shencount ~ s(elev, k = 6)
                    +s(std_IMI, k = 6)
                    + s(std_HLI, k = 6) 
                    + s(ppt, temp10, bs='ds', k =14, m=2)
                    + s(x, y, bs='ds', k =15, m=2) 
                    + s(ID,bs = "re")
                    ,
                    data=sit, family = "ziP", method = 'REML', offset = log(area/10000),
                    drop.unused.levels = FALSE)


shen_gam_eh35 = gam(shencount ~ s(elev, k = 6)
                    +s(std_IMI, k = 6)
                    + s(std_HLI, k = 6) 
                    + s(ppt, temp10, bs='ds', k =14, m=2)
                    + s(x, y, bs='ds', k =35, m=2) 
                    + s(ID,bs = "re")
                    ,
                    data=sit, family = "ziP", method = 'REML', offset = log(area/10000),
                    drop.unused.levels = FALSE)

## covariates with reduced data
par(mfrow=c(3,3))
p_objSa <- plot(shen_gam_eh, residuals = TRUE)
p_obj1a <- p_objSa[[1]] # just one smooth so select the first component-- Elevation
sm_dfa <- as.data.frame(p_obj1a[c("x", "se", "fit")])
data_dfa <- as.data.frame(p_obj1a[c("raw", "p.resid")])


sm25c <- sm_dfa
sm25c <- sm25c%>% 
  mutate(unscaled = plyr::mapvalues(round(x,2), from = round(elv$elev,2), to = elv$elev1))
sm25c <- sm25c %>% filter(unscaled >200)

(elev1 <- ggplot(sm25c, aes(x = unscaled, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Elevation (m)", y = "") +
    ggtitle("C. Subset by elevation, HLI and IMI"))

p_obj2a <- p_objSa[[2]] # just one smooth so select the second component-- IMI
sm_dfa <- as.data.frame(p_obj2a[c("x", "se", "fit")])
data_dfa <- as.data.frame(p_obj2a[c("raw", "p.resid")])


(imi1 <- ggplot(sm_dfa, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized IMI", y = " ") +ggtitle("C. Subset by elevation, HLI and IMI")+ ylim(-40,4) )



p_obj2a <- p_objSa[[3]] 
sm_dfa <- as.data.frame(p_obj2a[c("x", "se", "fit")])
data_dfa <- as.data.frame(p_obj2a[c("raw", "p.resid")])


(hli1 <- ggplot(sm_dfa, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "Standardized HLI", y = " ") +
    ggtitle("C. Subset by elevation, HLI and IMI") + ylim(-15,6))


plot(shen_gam_eh, select = 4, scale = 0, 
     ylab = "Temperature", xlab = "Precipitation", 
     shade = TRUE, shade.col = "gray80", lwd = 2 )

ehppt <- plot(shen_gam_eh, select = 4, scale = 0, 
              ylab = "Temp", shade = TRUE, shade.col = "gray80", lwd = 2)
eh4 <- ehppt[[4]]


eh4q <- as.data.frame(eh4[c("x",  "y")])


tempoeh <- expand.grid(eh4q$x, eh4q$y) 

temp2eh <- cbind(tempoeh, eh4$fit) %>% rename(x= Var1, y = Var2, fit = "eh4$fit" ) %>%
  filter(!is.na(fit))

#mycol <- c( "blueviolet","blue3" ,"turquoise3", "springgreen4" ,  "springgreen2" ,"aquamarine1",
# "gray73", "goldenrod4" , "gold3", "goldenrod2" , "orange" ,"orangered")

mycol <- c( "blueviolet","blue3" ,"turquoise3" ,"aquamarine1",
            "gray73", "yellow" , "gold1", "orange" , "orangered")

(ehpt <- ggplot(temp2eh, aes(x=x, y=y, z = fit)) + 
    geom_contour_filled(breaks = c(-1, -0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1,1.2)) + 
    xlab("3-day precipitation") + ylab("Temperature") + ggtitle("C. Subset by elevation, HLI, and IMI") +
    scale_fill_manual("Effect Size", values = mycol) + theme(legend.position = "none"))


##### map with reduced data


pred_gamSh <- predict(shen_gam_eh, data.frame(cov, area = rep(10000,nrow(cov))), se.fit = TRUE)


y = pred_gamSh$fit
x = pred_gamSh$se.fit
ZZ = "C. Subset by elevation, HLI and IMI"

elhi <- makemap(y=y, x=x, ZZ=ZZ)
#elhierr <- makemaperror(y=y, x=x, ZZ=ZZ)
pg2ehi <-   data.frame(predS = pred_gamSh$fit,seS = pred_gamSh$se.fit, 
                            lat = cov$lat, lon = cov$lon,x = cov$UTMx,
                            y = cov$UTMy, elev <- cov$elev, aspect = cov$aspect, 
                            IMI = cov$IMI, std_HLI = cov$std_HLI, left = cov$left_, right = cov$right_,
                            top = cov$top_, bottom = cov$bottom, predS1 =pred_gamS_elev$fit)


pg2ehi$status <- ifelse(exp(pg2ehi$predS) >1 , "greater than 1", "B")

#pg3 <- subset(pg2ehi, plogis(predS) > 0.1)
pg2ehi$LI <- exp(pg2ehi$predS1-2*pg2ehi$seS)

## how big each status is

pg2ehi %>% mutate(new_stat = ifelse(exp(predS1) >=1 & LI > 1, "pres", "abs")) %>%
  group_by(new_stat) %>% summarize(ct = n())

pg2ehi %>% group_by(status) %>% summarize(ct = n())
