library(ggplot2)
library(easyXpress)
library(xgboost)
source("./basics_TP.R")
library(readxl)
library(data.table)

#install.packages("devtools")


FILESIN = c("size1_Analysis-20240313.RData",
            "size2&3_Analysis-20240516.RData",
            "sizeassay4_Analysis-20241113.RData",
            "sizeassay5_6_Analysis-20241111.RData",
            "size8_Analysis-20250423.RData")

DIRPATHS = c("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size/",
             "/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size/",
             "/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size/",
             "/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size/",
             "/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size/")


wormpheno = NULL
for(i in 1:length(FILESIN)){
  filein = FILESIN[i]
  dirpath = DIRPATHS[i]
  
  data = Xpress(filedir=dirpath,rdafile=filein,well.id)
  data = subset(data[["processed_data"]], worm_length_um>300 & AreaShape_Area>1500 & model %in% c("L4_N2_HB101_100w", "L2L3_N2_HB101_100w"))
  
  if(sum(grepl("TPsizeassay3", data$FileName_RawBF))>0){
    
    wx = grepl("TPsizeassay3", data$FileName_RawBF)
    ff = data$FileName_RawBF[wx]
    data$Metadata_Date[wx]=tstrsplit(ff, "TPsizeassay3")[[1]]
    data$Metadata_Experiment[wx]="TPsizeassay3"
    
    plate=tstrsplit(ff, "-")[[2]]
    well=tstrsplit(ff, "-")[[3]]
    well=tstrsplit(well, "_")[[2]]
    well=tstrsplit(well, ".TIF")[[1]]
    
    data$Metadata_Magnification[wx]="m2X"
    data$Metadata_Group[wx]=paste0(plate,"_",well)
    data$Metadata_Plate[wx] = plate
    data$Metadata_Well[wx] = well
    data$well.id[wx]=paste0("TPsizeassay3_", plate, "_", well)
    
  }
  
  design <- read_excel("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size/meta/sizeassay_metadata.xlsx")
  design$Metadata_Well =  unlist(lapply(strsplit(design$Metadata_Well,""), function(x){
    if(length(x)<3) x = c(x[1], 0, x[2])
    x = paste(x, collapse = "")
    x
  }))
  data = merge(data, design)
  
  
  xmod = xgb.load( '/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size/sex_model/sexMale_xgb_preds')
  load( '/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size/sex_model/sexMale_xgb_preds.feature_names')
  
  
  data$probmale = predict(xmod, as.matrix(data[,feature_names]))
  data$sex[data$probmale>0.65] = "m"
  data$sex[data$probmale<0.35] = "f"
  
  
  xmod = xgb.load( '/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size/sex_model/wrongObjects_xgb_preds')
  load( '/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size/sex_model/wrongObjects_xgb_preds.feature_names')
  
  data$psuspicious = predict(xmod, as.matrix(data[,feature_names]))
  data$wrongObject = ifelse(data$psuspicious>0.5, T,F)
  
  wormpheno = rbind(wormpheno, data)
}



#write.csv(wormpheno, file = "/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size/becei_size_block1to8.csv", row.names = F)

#write.csv(wormpheno, file = "/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size/becei_size_block1to6.csv", row.names = F)

unique(wormpheno$Metadata_Experiment)

ggplot()+
  geom_histogram(data=subset(wormpheno, Metadata_Experiment=='TPsozeassay2' & Metadata_Plate == "p007" & Metadata_Well == "G09"),
                 aes(worm_length_um, color=sex, fill=sex), bins=50)+#theme_Publication2()+
  coord_cartesian(expand = 0)


ggplot()+
  geom_histogram(data=subset(wormpheno, Metadata_Experiment=='TPsozeassay2' & line == 3410 & !wrongObject),
                 aes(worm_length_um, color=sex, fill=sex), bins=50)+#theme_Publication2()+
  coord_cartesian(expand = 0)

library(ggplot2)
ggplot()+geom_histogram(data=data, aes(probmale, color=sex, fill=sex))+#theme_Publication2()+
  coord_cartesian(expand = 0)


ggplot(data=subset(wormpheno, wrongObject==F))+geom_histogram(aes(probmale, color=sex, fill=sex))+theme_Publication2()+
  coord_cartesian(expand = 0)+
  geom_vline(xintercept = 0.65, linetype='dashed', color='turquoise')+
  geom_vline(xintercept = 0.35, linetype='dashed', color="salmon")+
  xlab("Probability of being a male")


ggplot(wormpheno)+geom_histogram(aes(psuspicious, color=wrongObject, fill=wrongObject))+theme_Publication2()+
  coord_cartesian(expand = 0)+
  geom_vline(xintercept = 0.5, linetype='dashed', color='black')+
  scale_color_manual(values=c("darkgreen","red"))+
  scale_fill_manual(values=c("darkgreen","red"))+
  xlab("Probability of a wrong object")


x = subset(wormpheno, line == unique(line)[[25]] & wrongObject==F & !is.na(sex))
unique(subset(x, sex=='m' & worm_length_um > 800)$well.id)
unique(x$line)
ggplot(data=subset(x, wrongObject==F))+
  geom_histogram(aes(worm_length_um, color=sex), fill=NA)+
  coord_cartesian(expand = 0)

ggplot(data=subset(x, wrongObject==F))+
  geom_density(aes(worm_length_um, color=sex), fill=NA)+
  coord_cartesian(expand = 0)


x = subset(wormpheno, ((sex=='f' & worm_length_um < 750) | (sex=='m' & worm_length_um > 850)) & wrongObject==F)
sort(table(x$well.id), decreasing=T)[1:10]
View(subset(wormpheno, Metadata_Plate == "p001" & Metadata_Well=='B10' & Metadata_Experiment=="TPsozeassay2"))

viewWell3(df=wormpheno,
          img_dir="/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size_assay/Analysis-20240516/processed_images/",
          plate="p001", well="C09",experiment="TPsizeassayBloc", boxplot = F) 




#dt = read_excel("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size/sex_model/Book7.xlsx")
#dt$well.id = paste0(dt$Metadata_Experiment, "_", dt$Metadata_Plate, "_", dt$Metadata_Well)
#dt = dt[!duplicated(paste0(dt$well.id,'_',dt$Parent_WormObjects)),]
#dt = merge(dt, wormpheno[, -which(colnames(wormpheno) %in% c("sex", 'probmale', 'psuspicious', "wrongObject"))])

#writexl::write_xlsx(dt, "/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size/sex_model/imageXpress_sexTraining_humanAttributed2.xlsx")









