#################################################################
###############  Model Training #################################
#################################################################
library(xgboost)

library(readxl)
#source("./basics_TP.R")
dtrain_human <- read_excel("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size/sex_model/imageXpress_sexTraining_humanAttributed.xlsx")
dtrain_human2 <- read_excel("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size/sex_model/imageXpress_sexTraining_humanAttributed2.xlsx")
dtrain_human = rbind(dtrain_human,dtrain_human2[,match(colnames(dtrain_human),colnames(dtrain_human2))])

dtrain_fluo<- read_excel("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/size/sex_model/imageXpress_sexTraining_fluoAttributed.xlsx")
#design <- read_excel("rockmanlab/size_assay/sizeassay_metadata.xlsx")
#design$Metadata_Well =  unlist(lapply(strsplit(design$Metadata_Well,""), function(x){
#  if(length(x)<3) x = c(x[1], 0, x[2])
#  x = paste(x, collapse = "")
#  x
#}))

#sexflag = merge(sexflag, design)
dtrain_human$sex[!(dtrain_human$sex %in% c("m",'f'))]=NA
dtrain_human = subset(dtrain_human, !(well.id %in% c("TPSize1_p001_A07", "TPSize1_p001_A11","TPSize1_p001_C10")))
dtrain_fluo = dtrain_fluo[,match(colnames(dtrain_human), colnames(dtrain_fluo))]
dtrain_fluo$method = 'fluo'
dtrain_human$method = 'human'


#range(subset(dtrain_human, sex=='m' & CorrectDetection == "T")$worm_length_um)
#range(subset(dtrain_fluo, sex=='m')$worm_length_um)

dtrain = rbind(dtrain_human, dtrain_fluo)
colnames(dtrain)
ggplot(subset(dtrain[-which(dtrain$CorrectDetection == "F"),], !is.na(sex)),
       aes(Intensity_LowerQuartileIntensity_RawBF, color=sex))+
  geom_density()+facet_wrap(~method, ncol=1)

ggplot(subset(dtrain[-which(dtrain$CorrectDetection == "F"),], !is.na(sex)),
       aes(Intensity_UpperQuartileIntensity_RawBF,Intensity_LowerQuartileIntensity_RawBF, color=sex))+
  geom_point()+facet_wrap(~method, ncol=1)

##########################################################################################
### 1) train to recognize suspicious objects (undefined sex worms, uncorrect bounding box) 

require(xgboost)


traits = c("worm_length_um",
           "AreaShape_Area",
           "AreaShape_Perimeter",
           "AreaShape_Extent",
           "AreaShape_ConvexArea",
           "AreaShape_MinorAxisLength",
           "AreaShape_MajorAxisLength",
           "AreaShape_Compactness",
           "AreaShape_Eccentricity",
           "AreaShape_FormFactor",
           "AreaShape_MaximumRadius",
           "AreaShape_MeanRadius",
           "AreaShape_MinFeretDiameter",
           "AreaShape_MaxFeretDiameter",
           "AreaShape_EquivalentDiameter",
           "AreaShape_Compactness",
           "AreaShape_Solidity",
           "AreaShape_FormFactor",
           "AreaShape_Extent",
           "Intensity_MassDisplacement_RawBF",
           "Intensity_MeanIntensityEdge_RawBF",
           "Intensity_MADIntensity_RawBF",
           "Intensity_MedianIntensity_RawBF",
           "Intensity_MinIntensity_RawBF",
           "Intensity_IntegratedIntensity_RawBF",
           "Intensity_MeanIntensityEdge_RawBF",
           "Intensity_IntegratedIntensityEdge_RawBF",
           "Intensity_LowerQuartileIntensity_RawBF",
           "Intensity_UpperQuartileIntensity_RawBF")


dtrain_human$wrongObjects =  ifelse(!is.na(dtrain_human$sex) & as.logical(dtrain_human$CorrectDetection), 0,1)
#wtest = sample(1:nrow(dtrain), 250)
#dtest = dtrain[wtest,]
#dtrain = dtrain[-wtest,]

Xd = xgb.DMatrix(as.matrix(dtrain_human[,traits]), label = dtrain_human$wrongObjects  )


param <- list(max_depth = 4, eta = 0.3, subsample = 0.8, silent = 1, nthread = 4, objective = "binary:logistic", eval_metric = "error")
# x validation
cv <- xgb.cv(data = Xd, params = param, nfold=3, nrounds=120)

xmod = xgboost(params = param, data = Xd, nrounds = 500, early_stopping_rounds = 20)
(fimp = xgb.importance(feature_names = traits, model=xmod))

feature_names = traits
#save(feature_names, file = '/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/crossA/size/sex_model/wrongObjects_xgb_preds.feature_names')
#xgb.save(xmod, fname = '/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/crossA/size/sex_model/wrongObjects_xgb_preds')


# in sample predictions

th = 0.5
xpredt = predict(xmod, Xd)
table(dtrain_human$wrongObjects, xpredt>th)
sum(diag(table(dtrain_human$wrongObjects, xpredt>th)))/nrow(dtrain_human)


# dtest$suspicious = predict(xmod, as.matrix(dtest[,traits]))
# th = 0.5
# dtest$suspicious = ifelse(dtest$suspicious>0.5, T,F)
# dtest$wrongobject = ifelse(!is.na(dtest$sex) | as.logical(dtest$CorrectDetection), F, T)
# 
# good = sum(dtest$suspicious == dtest$wrongobject)
# falseDetection = sum(dtest$suspicious == T & dtest$wrongobject == F)
# undetected = sum(dtest$suspicious == F & dtest$wrongobject == T)
# 
#   
# test = data.frame(model = "XGB_wrongObjects", p = c(good,falseDetection,undetected)/nrow(dtest), match = factor(c("Correct","falseDetection","Undetected"), levels = c(c("Undetected","falseDetection","Correct"))))
# 
# ggplot(test, aes(model,p,fill=match))+
#   theme_Publication3()+
#   geom_bar(stat = "identity")+
#   scale_fill_manual(breaks = c("Correct","falseDetection","Undetected"), values = c("darkgreen", "salmon", "red"))+
#   scale_y_continuous(breaks = seq(0,1,0.1), expand = c(0,0))


##########################################################################################
### 2) train to recognize sexes ########################################################## 

require(xgboost)


traits = c("worm_length_um",
           "AreaShape_Area",
           "AreaShape_Perimeter",
           "AreaShape_Extent",
           "AreaShape_ConvexArea",
           "AreaShape_MinorAxisLength",
           "AreaShape_MajorAxisLength",
           "AreaShape_Compactness",
           "AreaShape_Eccentricity",
           "AreaShape_FormFactor",
           "AreaShape_MaximumRadius",
           "AreaShape_MeanRadius",
           "AreaShape_MinFeretDiameter",
           "AreaShape_MaxFeretDiameter",
           "AreaShape_EquivalentDiameter",
           "AreaShape_Compactness",
           "AreaShape_Solidity",
           "AreaShape_FormFactor",
           "AreaShape_Extent",
           "Intensity_MassDisplacement_RawBF",
           "Intensity_MeanIntensityEdge_RawBF",
           "Intensity_MADIntensity_RawBF",
           "Intensity_MedianIntensity_RawBF",
           "Intensity_MinIntensity_RawBF",
           "Intensity_IntegratedIntensity_RawBF",
           "Intensity_MeanIntensityEdge_RawBF",
           "Intensity_IntegratedIntensityEdge_RawBF",
           "Intensity_LowerQuartileIntensity_RawBF",
           "Intensity_UpperQuartileIntensity_RawBF")

dtrain = rbind(dtrain_human[,-which(colnames(dtrain_human)=="wrongObjects")], dtrain_fluo)
dtrain  = subset(dtrain , sex %in% c("m", "f") & !is.na(sex))
dtrain = dtrain[-which(dtrain$CorrectDetection == "F"),]

#wtest = sample(1:nrow(dtrain), 200)
#dtest = dtrain[wtest,]
#dtrain = dtrain[-wtest,]

Xd = xgb.DMatrix(as.matrix(dtrain[,traits]), label = as.integer(factor(dtrain$sex, labels=c(0,1)))-1)


param <- list(max_depth = 4, eta = 0.3, subsample = 0.8, silent = 1, nthread = 4, objective = "binary:logistic", eval_metric = "error")
# x validation
cv <- xgb.cv(data = Xd, params = param, nfold=3, nrounds=120)

xmod = xgboost(params = param, data = Xd, nrounds = 200, early_stopping_rounds = 20)
(fimp = xgb.importance(feature_names = traits, model=xmod))

feature_names = traits
#save(feature_names, file = '/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/crossA/size/sex_model/sexMale_xgb_preds.feature_names')
#xgb.save(xmod, fname = '/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/crossA/size/sex_model/sexMale_xgb_preds')


# in sample predictions
th = 0.5
xpredt = predict(xmod, Xd)
table(dtrain$sex, xpredt>th)
sum(diag(table(dtrain$sex, xpredt>th)))/nrow(dtrain)



#dtest$probmale = predict(xmod, as.matrix(dtest[,traits]))
#th = 0.75
#dtest$psex  = NA
#dtest$psex[dtest$probmale>th] = "m"
#dtest$psex[dtest$probmale<(1-th)] = "f"

#good = sum(dtest$sex == dtest$psex, na.rm=T)
#wrong = sum(dtest$sex != dtest$psex, na.rm=T)
#missing = nrow(dtest)-(good+wrong)


#test = data.frame(model = "XGB_sex", p = c(good,wrong,missing)/nrow(dtest), match = factor(c("Correct","Wrong","Missing"), levels = c(c("Missing","Wrong","Correct"))))

#ggplot(test, aes(model,p,fill=match))+
#  #theme_Publication3()+
#  geom_bar(stat = "identity")+
#  scale_fill_manual(breaks = c("Correct","Wrong","Missing"), values = c("darkgreen", "red", "grey"))+
#  scale_y_continuous(breaks = seq(0,1,0.1), expand = c(0,0))



# 
# 
# 
# #######################################
# ##### 1D ##############################
# #######################################
# 
# 
# minpheno  = NULL
# maxpheno = 1
# pheno = sexflag$AreaShape_Eccentricity
# 
# 
# 
# probmale1D = function(pheno, sex, minpheno=NULL, maxpheno=NULL){
#   #pheno = sexflag$worm_length_um
#   #sex = sexflag$sex
#   if(is.null(minpheno)) minpheno = min(pheno[which(sex %in% c("f", "m"))],na.rm=T)*0.8
#   if(is.null(maxpheno)) maxpheno = max(pheno[which(sex %in% c("f", "m"))],na.rm=T)*1.2
#   
#   v = seq(minpheno, maxpheno, (maxpheno-minpheno)/1000)
#   dmale = KernelDensityCustom(datavector = pheno[which(sex == "m")], valueout=v)
#   dfemale = KernelDensityCustom(datavector = pheno[which(sex %in% c("f"))], valueout=v)
#   
#   areamale = sum(((dmale[1:(nrow(dmale)-1),"density"] + d[2:nrow(dmale),"density"])/2) * diff(dmale$value), na.rm=T)
#   areafemale = sum(((dfemale[1:(nrow(dfemale)-1),"density"] + d[2:nrow(dfemale),"density"])/2) * diff(dfemale$value), na.rm=T)
#   
#   dmale$reldensity = dmale$density/areamale
#   dfemale$reldensity = dfemale$density/areafemale
#   
#   probmale =  dmale$reldensity/(dfemale$reldensity+dmale$reldensity)
#   probmale[(dmale$quantile<0.01 | dmale$quantile>0.99) & (dfemale$quantile<0.01 | dfemale$quantile>0.99)]=NA
#   
#   probmale = data.frame(value=dmale$value, probmale = probmale)
#   probmale=probmale[!is.na(probmale$probmale),]
#   colnames(probmale) = c("pheno", "probmale")
#   return(probmale)
# }
# 
# 
# 
# 
# 
# 
# 
# probmale2D = function(pheno, sexes){
#   require(MASS)
#   #pheno = as.data.frame(sexflag[,c("worm_length_um", "AreaShape_ConvexArea")])
#   #sexes = sexflag$sex
#   ranges = c(apply(pheno[which(sexes %in% c("f", "m")),], 2, range))
#   
#   # 2D Density 
#   ds = lapply(c("m","f"), function(sex){
#     phenox = pheno[which(sexes == sex),]
#     density <- MASS::kde2d(x=phenox[,1],y=phenox[,2], ,n=100, lims = ranges ) 
#     names(density)[3] = paste0("density_", sex)
#     return(density)
#   })
#   
#   ds = list(x=ds[[1]][[1]], y=ds[[1]][[2]],density_m = ds[[1]][[3]], density_f = ds[[2]][[3]])
#   ds$probmale = ds$density_m/(ds$density_m + ds$density_f)
#   
#   #Quantiles
#   qs = lapply(c("m", "f"), function(sex){
#     xbin = diff(sort(unique(ds$x)))[1]
#     ybin = diff(sort(unique(ds$y)))[1]
#     
#     x = c(ds[[paste0("density_", sex)]])
#     x=x/sum(x * xbin *ybin)
#     
#     rd = range(x)
#     steps = seq(rd[1], rd[2], (rd[2]-rd[1])/1000)
#     q = data.frame(density = steps,quantile = sapply(steps, function(s){sum(x[x<s] * xbin *ybin)}))
#     out = approx(x=q$density, y=q$quantile, xout = x)$y
#     out = array(out, dim = dim(ds[[paste0("density_", sex)]]))
#     out
#   })
#   
#   names(qs) = c("quantile_m", "quantile_f")
#   
#   filterout = qs[["quantile_m"]] < 0.001 | qs[["quantile_f"]] < 0.001
#   
#   ds[["probmale"]][filterout] = NA
#   
#   probmale = ds[c("x", "y", "probmale")]
#   
#   traitname = c(x = colnames(pheno)[1], y = colnames(pheno)[2])
#   
#   probmale[[4]]  = traitname
#   
#   names(probmale) = c("x","y","probmale","traitname")
#   
#   return(probmale)
# }
# 
# 
# ###################################
# ######### 3D ######################
# ###################################
# 
# probmale3D = function(pheno, sexes){
#   library(plotly)
#   
#   #pheno = as.data.frame(sexflag[,c("worm_length_um", "AreaShape_ConvexArea", "AreaShape_MinorAxisLength", "sex")])
#   #sexes = sexflag$sex
#   
#   ranges = c(apply(pheno[which(sexes %in% c("f", "m")),], 2, range))
#   
#   ds = lapply(c("m","f"), function(sex){
#     phenox = pheno[which(sexes == sex),]
#     density <- misc3d::kde3d(x=phenox[,1],y=phenox[,2], z = phenox[,3], ,n=75, lims = ranges ) 
#     names(density)[4] = paste0("density_", sex)
#     return(density)
#   })
#   
#   ds = list(x=ds[[1]][[1]], y=ds[[1]][[2]], z=ds[[1]][[3]],density_m = ds[[1]][[4]], density_f = ds[[2]][[4]])
#   
#   
#   ds$probmale = ds$density_m/(ds$density_m + ds$density_f)
#   
#   qs = lapply(c("m", "f"), function(sex){
#     xbin = diff(sort(unique(ds$x)))[1]
#     ybin = diff(sort(unique(ds$y)))[1]
#     zbin = diff(sort(unique(ds$z)))[1]
#     
#     x = c(ds[[paste0("density_", sex)]])
#     x=x/sum(x * xbin *ybin *zbin)
#     
#     rd = range(x)
#     steps = seq(rd[1], rd[2], (rd[2]-rd[1])/1000)
#     q = data.frame(density = steps,quantile = sapply(steps, function(s){sum(x[x<s] * xbin*ybin*zbin)}))
#     out = approx(x=q$density, y=q$quantile, xout = x)$y
#     out = array(out, dim = dim(ds[[paste0("density_", sex)]]))
#     out
#   })
#   
#   names(qs) = c("quantile_m", "quantile_f")
#   
#   filterout = qs[["quantile_m"]] < 0.001 | qs[["quantile_f"]] < 0.001
#   
#   ds[["probmale"]][filterout] = NA
#   
#   probmale = ds[c("x", "y", "z", "probmale")]
#   
#   traitname = c(x = colnames(pheno)[1], y = colnames(pheno)[2], z = colnames(pheno)[3])
#   
#   probmale[[5]]  = traitname
#   
#   names(probmale) = c("x","y","z","probmale","traitname")
#   
#   return(probmale)
# }
# 
# 
# 
# 
# 
# predictSex = function(traits ,probmale){
#   coordx = unlist(lapply(1:length(traits), function(i){
#     which.min(abs(as.numeric(probmale[[i]]) - traits[i]))
#   }))
#   
#   coordx = matrix(coordx,nrow=1)
#   probmale[["probmale"]][coordx]
# }
# 
# 
# ######################################################
# ######################################################
# ######################################################
# 
# 
# 
# 
# ######################################################
# ######################################################
# ######################################################
# 
# 
# 
# # dtest$psexXGB = predict(xmod, as.matrix(dtest[,traits]))
# # th = 0.75
# # 
# # predsex = rep(NA, nrow(dtest))
# # predsex[dtest[,"psexXGB"] > th] = "m"
# # predsex[dtest[,"psexXGB"] < 1-th] = "f"
# # good = c(sum(predsex ==dtest$sex, na.rm = T),sum(predsex != dtest$sex, na.rm = T))
# # test = data.frame(model = model, p = c(good, nrow(dtest)-sum(good))/nrow(dtest), match = factor(c("correct","wrong","missing"), levels = c(c("missing","wrong","correct"))))
# # 
# # 
# # ggplot(test, aes(model,p,fill=match))+
# #   theme_Publication3()+
# #   geom_bar(stat = "identity")+
# #   scale_fill_manual(breaks = c("correct", "wrong", "missing"), values = c("darkgreen", "#E64B35FF", "grey"))
# 
# 
# 
# 
# # dtrain = sexflag
# # custom2D = probmale2D(pheno=as.data.frame(dtrain[,c("worm_length_um", "AreaShape_ConvexArea")]),sexes = sexflag$sex)
# # custom3D = probmale3D(pheno=as.data.frame(dtrain[,c("worm_length_um", "AreaShape_ConvexArea", "AreaShape_MinorAxisLength")]),sexes = sexflag$sex)
# #  
# # #ggplot(subset(sexflag, sex %in% c("m","f")), aes(worm_length_um,AreaShape_ConvexArea,color=sex))+
# # #  geom_point()+scale_color_manual(values = c("#E64B35FF","#440154FF"))+theme_Publication3()
# # 
# # #plot_ly(x=dtrain[,"worm_length_um"], y=dtrain[,"AreaShape_ConvexArea"], z=dtrain[,"sex"], type="scatter3d", mode="markers", size=0.5, color=dtrain[,"sex"], colors = c("#E64B35FF","#440154FF"))
# # 
# # 
# # 
# # probmale = custom2D[["probmale"]]
# # 
# # dimnames(probmale)[[1]] = custom2D[["x"]]
# # dimnames(probmale)[[2]] = custom2D[["y"]]
# # 
# # probmale = reshape2::melt(probmale)
# # probmale=probmale[!is.na(probmale$value),]
# # colnames(probmale) = c("x","y","probmale")
# # 
# # ggplot(probmale, aes(x,y, color=probmale))+geom_point()+theme(legend.position = "none")+
# #   scale_colour_gradient2(high = "#440154FF",low= "#E64B35FF", mid = "purple", midpoint = 0.5)+
# #   xlim(range(custom2D$x))+ylim(range(custom2D$y))+theme_Publication3()
# # 
# # 
# # probmale = custom3D[["probmale"]]
# # 
# # dimnames(probmale)[[1]] = custom3D[["x"]]
# # dimnames(probmale)[[2]] = custom3D[["y"]]
# # dimnames(probmale)[[3]] = [["z"]]
# # 
# # probmale = reshape2::melt(probmale)
# # probmale=probmale[!is.na(probmale$value),]
# # 
# # plot_ly(x=probmale[,1], y=probmale[,2], z=probmale[,3], type="scatter3d", mode="markers", color=probmale[,4], colors = c("#E64B35FF","#440154FF"))
# # 
# # 
# # 
# # 
# # sexflag$psexXGB = predict(xmod, as.matrix(sexflag[,traits]))
# # sexflag$psex2D = apply(sexflag[,custom2D$traitname],1,function(x){predictSex(traits=x ,probmale=custom2D)})
# # sexflag$psex3D = apply(sexflag[,custom3D$traitname],1,function(x){predictSex(traits=x ,probmale=custom3D)})
# # 
# # th = 0.75
# # test = do.call(rbind, lapply(c("XGB", "2D", "3D"), function(model){
# #   wx = paste0("psex", model)
# #   xx = subset(sexflag, sex %in% c("m","f"))
# #   predsex = rep(NA, nrow(xx))
# #   predsex[xx[,wx] > th] = "m"
# #   predsex[xx[,wx] < 1-th] = "f"
# #   good = c(sum(predsex == xx$sex, na.rm = T),sum(predsex != xx$sex, na.rm = T))
# #   data.frame(model = model, n = c(good, nrow(xx)-sum(good)), match = factor(c("correct","wrong","missing"), levels = c(c("missing","wrong","correct"))))
# #   
# # }))
# # 
# # ggplot(test, aes(model,n,fill=match))+
# #   theme_Publication3()+
# #   geom_bar(stat = "identity")+
# #   scale_fill_manual(breaks = c("correct", "wrong", "missing"), values = c("darkgreen", "#E64B35FF", "grey"))
# # 
# # 
