library(ggplot2)
library(mgcv)
library(rgdal)
library(dplyr)
library(raster)
library(voxel)
library(gridExtra)
library(terra)

pred_weather <- read.csv("weath_cov.csv")
sh.plot <- read.csv("non_geo_data.csv") ## this is not real location data, this will not replicate paper results
cov <- read.csv("covariate_non_geo.csv")
source("Graphing_Set_Up.R")
source("map_functions.R")
'%notin%' <- Negate('%in%')
cov$area <- log(10000)

cov <- cov %>% mutate(asp_rad = pi*aspect/180) %>%
  mutate(north = cos(asp_rad), east = sin(asp_rad)) %>%
  mutate(std_HLI = scale(HLI)[,1], std_IMI = scale(IMI)[,1])

cov <- cbind(cov, pred_weather[,c("ppt3", "temp10")])

cov <- cov %>%
  rename("ppt" = "ppt3")

elv <- cov %>% dplyr::select(elev, elev1)


################# set up ##############################
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

############## Subset across elev ---------------
### subsets -- assume only can sample 1/3 of all sites 
## across elevation (transects)

source("subset_elevation.R")

rm(out,outpred, shen_gam_elev15, shen_gam_elev35, p_objS,p_obj1, sm_df, data_df, p_obj2)
gc()


########## Subset across Elevation, HLI and IMI ------------- 

source("subset_elev_hli_imi.R")

rm(out,outpred, shen_gam_eh15, shen_gam_eh35, p_objSa,p_obj1a, sm_dfa, data_dfa, p_obj2a)
gc()




##### random site selection ----
#source("subset_random.R")
#rm(out,outpred, shen_gam_rand15, shen_gam_rand35, p_objSabr ,p_obj1abr, sm_dfabr, data_dfabr, p_obj2abr)
#gc()





######### only one visit per site per year -----------------

source("subset_annual_sampling.R")
rm(out,outpred, shen_gam_vis200, shen_gam_vis25, p_objSab ,p_obj1ab, sm_dfab, data_dfab, p_obj2ab)
gc()



####### subset within 250m (ie. 15 min walk) of trail or road ----------------

source("subset_dist_road_trail.R")



############ Arrange graphs -----------------
# read in full dataset model fits so can graph together for easy comparison
fin_mod_out <- readRDS("Pshen_mod_outzIP_aug9.RDS")
fds <- fin_mod_out[[3]]
rm(fin_mod_out)
gc()
 

p_objS <- plot(fds, residuals = TRUE)
p_obj1 <- p_objS[[1]] # just one smooth so select the first component-- Elevation
sm_df <- as.data.frame(p_obj1[c("x", "se", "fit")])
data_df <- as.data.frame(p_obj1[c("raw", "p.resid")])


sm_df <- sm_df%>% 
  mutate(unscaled = plyr::mapvalues(round(x,2), from = round(elv$elev,2), to = elv$elev1))%>% 
  filter(unscaled >200)


(elevF <- ggplot(sm_df, aes(x = unscaled, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "", y = "Effect Size") +
    ggtitle("A. Full dataset")) #+ ylim(-40, 20) )

p_obj2 <- p_objS[[2]] # just one smooth so select the second component-- IMI
sm_df <- as.data.frame(p_obj2[c("x", "se", "fit")])
data_df <- as.data.frame(p_obj2[c("raw", "p.resid")])


(imiF <- ggplot(sm_df, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = "", y = "Effect Size") +ggtitle("A. Full dataset") + ylim(-50, 15) )



p_obj3 <- p_objS[[3]] # just one smooth so select the third component-- HLI
sm_df <- as.data.frame(p_obj3[c("x", "se", "fit")])
data_df <- as.data.frame(p_obj3[c("raw", "p.resid")])


(hliF <- ggplot(sm_df, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
                alpha = 0.1) +
    geom_line(lwd= 1) +
    labs(x = " ", y = "Effect Size") +ggtitle("A. Full dataset") + ylim(-15,6) )  




#grid.arrange(elevF, imiF, hliF, elevS,imi, hli, elev1, imi1, hli1,elev3, imi3, hli3, elev2, imi2, hli2,  nrow = 5 )

grid.arrange(elevF, elevS, elev1, elev2, elevD, nrow =2)
grid.arrange(imiF, imi, imi1, imi2, imiD,nrow =2)
grid.arrange(hliF, hli, hli1, hli2, hliD, nrow =2)




pl4 <- plot(shen_gam200, select = 4, scale = 0, 
            ylab = "Temp", shade = TRUE, shade.col = "gray80", lwd = 2)
ps4 <- pl4[[4]]


ps4q <- as.data.frame(ps4[c("x",  "y")])


tempo <- expand.grid(ps4q$x, ps4q$y) 

temp2 <- cbind(tempo, ps4$fit) %>% rename(x= Var1, y = Var2, fit = "ps4$fit" ) %>%
  filter(!is.na(fit))


mycol <- c( "blueviolet","blue3" ,"turquoise3" ,"aquamarine1",
            "gray73", "yellow" , "gold1", "orange" , "orangered")

(stp <- ggplot(temp2, aes(x=x, y=y, z = fit)) + 
    geom_contour_filled(breaks = c(-1, -0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1,1.2)) + 
    xlab("3-day precipitation") + ylab("Temperature") + ggtitle("A.Full Data Set") +
    scale_fill_manual("Effect Size", values = mycol) + theme(legend.position = "none") )


grid.arrange(stp ,elevpt, ehpt,annpt,roadpt)

##### make maps -----

# predicted counts
pred_gamS_full<- predict(fds, data.frame(cov, area = rep(10000,nrow(cov))), se.fit = TRUE)

y = pred_gamS_full$fit
x = pred_gamS_full$se.fit
ZZ = "A. Full dataset"


pg2 <- pg <- data.frame(predS = y,seS = x, 
                        lat = cov$lat, lon = cov$lon,x = cov$UTMx,
                        y = cov$UTMy, elev <- cov$elev, aspect = cov$aspect, 
                        IMI = cov$IMI, std_HLI = cov$std_HLI, left = cov$left_, right = cov$right_,
                        top = cov$top_, bottom = cov$bottom)



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


(full <- ggplot() +  
           geom_tile(data = pg3, aes(x = x, y = y, fill = dummy, alpha=erDum)) + 
                             
           scale_fill_manual("P. shenandoah \n predicted counts",
                             
                             values = cols.d$colo,
                             labels = c(1000, 500, 250, 125,75,35,18,9,4,1)) +
           
           
           ylab("Latitude") + xlab("Longitude") + ggtitle(ZZ)  +
          # theme(legend.position = c(0.2, 0.8), legend.direction = "horizontal") + 
           theme(legend.position = "none")+
          scale_alpha(guide = 'none') + 
           theme(axis.text.x = element_blank(),
                 axis.text.y = element_blank())  )


grid.arrange(full, elm, elhi, annual, elmD)







