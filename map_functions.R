makemap <- function(y, x, ZZ){
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
  
  (presg<- ggplot() +  
      geom_tile(data = pg3, aes(x = x, y = y, fill = dummy, alpha=erDum)) + 
      scale_fill_manual(values = cols.d$colo) +
      
      ylab("Latitude") + xlab("Longitude") + ggtitle(ZZ)  +
      theme(legend.position = "none") + scale_alpha(guide = 'none') + 
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank()) )  }

