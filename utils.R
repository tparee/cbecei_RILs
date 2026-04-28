##############################################
########## TOM BASICS FUNCTUION & Co #########
##############################################

foundercolors = c(`FA.g1` = '#9CC891',
                  `FA.g2` = '#2A6218',
                  `FM.g1` = '#6DB9D2',
                  `FM.g2` = '#0E1A52',#'#415384',
                  `FB.g1` = '#DFB640',
                  `FB.g2` = "#C1443E")



######################################################
########### ggplot2 ##################################

# Personal theme

theme_perso=function(){
  require(ggplot2)
  theme_minimal()+
    theme(
      strip.text.x = element_text(color="black",face="bold"),
      strip.text.y = element_text(color="black",face="bold"),
      axis.title=element_text(face="bold"),
      axis.text=element_text(colour = "black"),
      panel.border = element_rect(colour = "black", fill=NA, size=0.4))
}


# Theme for publication from https://rdrr.io/github/HanjoStudy/quotidieR/f/

theme_Publication <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


theme_Publication2 <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour=NA,fill=NA),
            strip.text = element_text(face="bold")
    ))
  
}



theme_Publication3 <- function(base_size=7) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold",
                                      hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(fill = "white", colour = NA),
            plot.background  = element_rect(fill = "white", colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "plain"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(t=5,r=5,b=5,l=5),"pt"),
            strip.background=element_rect(colour=NA,fill=NA),
            strip.text = element_text(face="plain")
    ))
  
}



scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

# Some color to set up
cbs = c(rgb(0,0.45,0.7),
        rgb(0.9, 0.6, 0),
        rgb(0,0.6,0.5),
        rgb(0.8, 0.4, 0),
        rgb(0.35, 0.7, 0.9),
        rgb(0.8, 0.6, 0.7),
        rgb(0.95, 0.9, 0.25),
        rgb(0,0,0,),
        "gray45")


#ggplot(data.frame(x = 1:length(cbs), y =1:length(cbs),cc = cbs), 
#       aes(x=x, y =y, color=cc))+geom_point(size=3)+
#  scale_color_manual(breaks = cbs, values = cbs)+
#  theme_minimal()


# Personal color        
scale_colour_perso <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = cbs), ...)
  
}

scale_fill_perso <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = cbs), ...)
}




get.win = function(size, winsize, minsize){
  
  if(winsize < size){
    
    nblocks = size/winsize
    fullblocks = floor(nblocks)
    
    ## split column numbers into 'nblocks' groups
    SPLIT <- split(1:(fullblocks*winsize), rep(1:fullblocks, each = winsize))
    if(nblocks>fullblocks){ 
      dblock = (nblocks-fullblocks)*winsize
      SPLIT[[length(SPLIT)+1]] = 1:dblock + max(fullblocks*winsize)
    }
    
    if(length(SPLIT[[length(SPLIT)]]) < minsize & length(SPLIT) > 1){
      SPLIT[[length(SPLIT)-1]] = c(SPLIT[[length(SPLIT)]],SPLIT[[length(SPLIT)-1]])
      SPLIT = SPLIT[1:(length(SPLIT)-1)]
    }
    
    
  }else{
    
    SPLIT = list(1:size)
  }
  
  return(SPLIT)
}


se = function(x,na.rm=F){if(na.rm){x=x[!is.na(x)]};sd(x)/sqrt(length(x))}



KernelDensityCustom=function(datavector, valueout=NULL){
  # modified from https://stats.stackexchange.com/questions/108342/how-to-get-percentiles-from-empirical-density-in-r
  
  # Default: return density at each datvector value
  if(is.null(valueout)){i = datavector}else{i=valueout}
  
  hT = bw.nrd0(datavector)
  
  KDE <-  Vectorize(function(x){mean(pnorm((x-datavector)/hT))}) #CDF
  kde <- Vectorize(function(x) mean(dnorm((x-datavector)/hT)/hT)) # DF
  
  out = data.frame(value = i, density = kde(i), quantile=KDE(i))
  
  return(out)
  
  # Input = quantiles and return value
  # Interval = c(min(datavector), max(datavector))
  # quantiles = c(0.05, 0.5, 0.9)
  # cdf = KDE(datavector)
  # if(sum(quantiles < min(cdf) | quantiles >max(cdf))>0){
  #   quantiles[quantiles < min(cdf)]  = min(cdf)
  #   quantiles[quantiles > max(cdf)] = max(cdf)
  # }
  # out = sapply(quantiles, function(p){
  #   print(p)
  #   tempf <- function(t) KDE(t)-p
  #   try(uniroot(tempf,Interval)$root)
  # })
  # 
  # names(out) = quantiles
  # return(out)
  
}



localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}



#weirdFigureRatio = 1.5



scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}






relative = function(x){(x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))}
