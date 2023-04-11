.libPaths('/data/zool-avianmacro/jesu3393/R/libraryc/')

library(packfor)
library(neutral.vp)
library(ade4)
library(dplyr)
library(colorBlindness)
library(vegan)

load('/data/zool-avianmacro/jesu3393/papers/Globallocal/nhd.RData')
################################################################
###community conversion based on observer ability

comc <- function(sampdf, a, s){
  
  all0row <- which(rowSums(sampdf[,]) == 0)###get rows with no species
  
  sampdf[all0row,] <- 1 ##set to all 1 all the subsequent calculation
  
  set.seed(2021)
  
  sampdf <- apply(sampdf, 1, function(x)sample(length(x), a, replace = T, prob = x))###observations by abundance
  
  sampdf <- lapply(1:s, function(x){##turn back to community data
    
    subsample <- data.frame(table(sampdf[,x]))
    
    subsample$Var1 <- as.character(subsample$Var1)
    
    colnames(subsample)[2] <- paste("x",x,sep = "")###different site
    
    return(subsample)
    
  })
  
  sampdf <- Reduce(function(...)merge(..., by = "Var1",all =T), sampdf)
  
  sampdf[is.na(sampdf)] <- 0
  
  rownames(sampdf) <- paste("s",sampdf$Var1,sep = "")###different species
  
  sampdf$Var1 <- NULL
  
  all0row <- paste("x",all0row,sep = "")
  
  sampdf[,all0row] <- 0
  
  sampdf[,"x"] <- NULL
  
  sampdf <- data.frame(t(sampdf))
  
  return(sampdf)
  
}

rnumb <- which(1:40000%%200 %in% c(26:175)) ##retain the central 150*150 sites

#
rnumb <- rnumb[rnumb > 5000 & rnumb < 35001]

#####################################################################
######temporal analysis

neut <- neut.nhd

y = 50##separate year calculation 

indi.dist0 <- census(neut, snap = 50)

indi.dist0 <- indi.dist0[rnumb,]

indi.dist0 <- data.frame(indi.dist0)

s4rnumber <- do.call(c, lapply(1:75, function(i)rep((75*i-74):(75*i), each = 2, times = 2)))

indi.dist0 <- rowsum(indi.dist0, as.integer(s4rnumber))

#s9rnumber <- do.call(c, lapply(1:60, function(i)rep((60*i-59):(60*i), each = 3, times = 3)))

#indi.dist0 <- rowsum(indi.dist0, as.integer(s9rnumber))

allgens <- lapply(seq(50,149,1), function(y){#get data for 100 generations
  
  cat(y,'\n')
  
  indi.dist <- census(neut, snap = y)
  
  indi.dist <- indi.dist[rnumb,]
  
  sampdf <- data.frame(indi.dist)
  
  sampdf <- rowsum(sampdf, as.integer(s4rnumber))
  
  sample.spabun <- lapply(c(50,200,500), function(a){
    
    cat(a,'\n')
    
    indi.dist0 <- comc(indi.dist0, a, s = nrow(indi.dist0))
    
    sampdf <- comc(sampdf, a, s = nrow(sampdf))
    
    betd <- lapply(1:nrow(indi.dist0), function(x){###beta diversity
      
      nc <- plyr::rbind.fill(indi.dist0[x,],sampdf[x,])###compared with the first year
      
      nc[is.na(nc)] <- 0
      
      bd1 <- vegdist(nc, method = 'jaccard', binary = T)###presence-absence data
      
      bd2 <- vegdist(nc, method = 'horn', binary = T)#
      
      bd <- cbind(bd1=bd1,bd2=bd2)
      
      return(data.frame(bd))
      
    })
    
    
    betd <- bind_rows(betd)
    ##div.eachsite <- apply(sampdf, 2, n_distinct)
    
    ##div.eachsite <- data.frame(div=div.eachsite, obserabun = a) 
    
    div.eachsite <- data.frame(sp=specnumber(sampdf, MARGIN = 1), beta = betd$bd1, beta2 = betd$bd2, site = rownames(sampdf), obserabun = a)##other indices can be included here.
    
    ##div.eachsite[all0row,"div"] <- 0 ###no species sites diversity is set to 0
    
    #return(div.eachsite)
    
  #})
  
  #div.year <- bind_rows(sample.spabun)
  
  #div.year$year <- y-49
  
  #return(div.year)
    
    div.eachsite$year <- y-49
    
    return(div.eachsite)
  
})


timeall <- bind_rows(sample.spabun)

timeall$logsp <- log10(timeall$sp)

return(timeall)

})


timeall <- bind_rows(timeall)


saveRDS(timeall,file = paste("/data/zool-avianmacro/jesu3393/papers/Globallocal/s4divbyyear3abun/n4y",".RDS",sep = ""))

rm(neut.nhd)
rm(neut)
rm(timelist)
rm(tplotdata)

save.image(paste('/data/zool-avianmacro/jesu3393/papers/Globallocal/s4divbyyear3abun/n4y','.RData', sep = ""))



