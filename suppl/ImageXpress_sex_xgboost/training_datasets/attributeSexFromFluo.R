library(data.table)
library(writexl)

relative = function(x){(x-min(x, na.rm = T))/(max(x, na.rm=T)-min(x, na.rm=T))}

calculate_distance <- function(x1, y1, x2, y2) {
  distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  return(distance)
}

MAINPATH =  "~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/crossA/size/male_model/imgXpress_Fluo"

setwd(MAINPATH)

# Load fluorescent marker positions on imageXpress of sex-specific fluorescent worms
# Adults females with mCherry (w3) spermathecae (two dots), adult males with GFP (w2) vas deferens (one dot)
# i.e, F1 progeny of QG4602 and another QG becei strain (cross A). Not really important here
# Fluo dots are detected through "fluo_sexratio" pipeline on image taken by image Xpress
load("fluoDotsPositions.Rdata")

#  head(fluoDotsPositions)
#        X       Y Mean_Intensity Min_Value Max_Value wavelength flag.lowcontrast flag.clusters  sex                          imgID
#1 1311.539 100.645      1666.2692       630      2360         w2            FALSE         FALSE male 241107-TPsexratio-p001-m2X_A01
#2 1076.730 228.533       849.5769       401      1152         w2            FALSE         FALSE male 241107-TPsexratio-p001-m2X_A01

fluoDotsPositions$well = tstrsplit(fluoDotsPositions$imgID, "m2X_")[[2]]
fluoDotsPositions$plate = tstrsplit(fluoDotsPositions$imgID, "-")[[3]]
fluoDotsPositions$experiment = tstrsplit(fluoDotsPositions$imgID, "-")[[2]]

# Load and process the worm objects detected by cellprofiler Andersen lab's pipeline, on the same imageXpress photos
cellproOutputFile = list.files("cp_data")
cellprodata = Xpress(filedir=".",rdafile=cellproOutputFile,well.id)
cellprodata = cellprodata[["processed_data"]]
#cellprodata = subset(cellprodata , grepl('_w1.TIF', FileName_RawBF))
cellprodata = subset(cellprodata, grepl(paste("L2L3|L4"), model))
cellprodata = subset(cellprodata, !(grepl("C|E|F",Metadata_Well) & grepl("p001",Metadata_Plate)))

# Now we want to attribute the sex of worm objects in cellprodata
# by looking at the fluorescent markers within the worm objects

cellprodata = do.call(rbind, lapply(split(cellprodata, cellprodata$well.id), function(well_data){
  
  # well_data = split(cellprodata, cellprodata$well.id)[[110]]
  thisimgID = well_data$well.id[1]
  print(thisimgID)
  thisexperiment = tstrsplit(thisimgID, "_")[[1]]
  thisplate = tstrsplit(thisimgID, "_")[[2]]
  thiswell = tstrsplit(thisimgID, "_")[[3]]
 
  pos = subset(fluoDotsPositions, well == thiswell & plate == thisplate & experiment == thisexperiment)
  
  
  out = do.call(rbind, lapply(1:nrow(well_data), function(i){
    
    #print(i)
    thisobject = well_data[i,]
    
    radius = 10
    
    # Extract the coordinate of the worm object
    
    X = unlist(c(thisobject[grepl("Worm_ControlPointX_",colnames(thisobject))]))
    X=X[order( as.numeric(data.table::tstrsplit(names(X), "_")[[3]]))]
    
    Y = unlist(c(thisobject[grepl("Worm_ControlPointY_",colnames(thisobject))]))
    Y=Y[order( as.numeric(data.table::tstrsplit(names(Y), "_")[[3]]))]
    
    wormcoord = data.frame(X=X, Y=Y)
    
    fluomatch = do.call(rbind, lapply(1:nrow(wormcoord), function(j){
      d = calculate_distance(x1=wormcoord$X[j], y1=wormcoord$Y[j], x2=pos$X, y2=pos$Y)
      pos[which(d<radius),]
    }))
    
    fluomatch =  fluomatch[!duplicated(paste0( fluomatch$X,"_",fluomatch$Y)),]
    
    if(nrow(fluomatch)==0 | nrow(fluomatch) > 2){return(NULL)}
    if(length(unique(fluomatch$sex)) > 1){return(NULL)}
   
    sex = c("m", "f")[match(paste(fluomatch$sex, collapse='_'), c("male", "female_female"))]
    
    if(is.na(sex)){return(NULL)}
    thisobject$sex = sex
    thisobject
  }))
  
  return(out)
}))


writexl::write_xlsx(cellprodata, "imageXpress_sexTraining_fluoAttributed.xlsx")

#select image to be analyzed
img <- tiff::readTIFF("/Users/tomparee/Downloads/241107-TPsexratio-p002-m2X_H02_w1.TIF")

h<-dim(img)[1] # image height
w<-dim(img)[2] # image width
colnames(well_data)

well_data %>%
  ggplot2::ggplot(.) +
  ggplot2::aes(x = AreaShape_Center_X, y = AreaShape_Center_Y) +
  ggplot2::annotation_custom(grid::rasterGrob(relative(img), width=ggplot2::unit(1,"npc"), height=ggplot2::unit(1,"npc")), 0, w, 0, -h) +
  #ggplot2::geom_point(data=wormcoord, aes(x = X, y = Y), color='red')+
  ggplot2::geom_point(data=out, aes(x = Worm_ControlPointX_1, y = Worm_ControlPointY_1), color='yellow', shape=1)+
  ggplot2::geom_point(data=out, aes(x = Worm_ControlPointX_5, y = Worm_ControlPointY_5), color='yellow', shape=1)+
  ggplot2::geom_point(data=out, aes(x = Worm_ControlPointX_10, y = Worm_ControlPointY_10), color='yellow', shape=1)+
  ggplot2::geom_point(data=out, aes(x = Worm_ControlPointX_15, y = Worm_ControlPointY_15), color='yellow', shape=1)+
  ggplot2::geom_point(data=out, aes(x = Worm_ControlPointX_21, y = Worm_ControlPointY_21), color='yellow', shape=1)+
  ggplot2::geom_text(data=out, aes(x = AreaShape_Center_X, y = AreaShape_Center_Y, color=sex, label = Parent_WormObjects), size=7)+
  ggplot2::geom_rect(aes(xmin = po_AreaShape_BoundingBoxMinimum_X, xmax = po_AreaShape_BoundingBoxMaximum_X, ymin = po_AreaShape_BoundingBoxMinimum_Y, ymax = po_AreaShape_BoundingBoxMaximum_Y), color='skyblue', fill = NA, linetype='dashed')+
  ggplot2::geom_rect(aes(xmin = AreaShape_BoundingBoxMinimum_X, xmax = AreaShape_BoundingBoxMaximum_X, ymin = AreaShape_BoundingBoxMinimum_Y, ymax = AreaShape_BoundingBoxMaximum_Y), color='skyblue', fill = NA)+
  ggplot2::scale_x_continuous(expand=c(0,0),limits=c(0,w)) +
  ggplot2::scale_y_reverse(expand=c(0,0),limits=c(h,0)) +
  ggplot2::coord_equal() +
  ggplot2::theme_bw()
  

