{ #libraries and options
options("width"=200) ; options(scipen=99) ; library('parallel')#Wide screen ; no scientific output
library('kSamples') ; library('extraDistr')  ## for Anderson–Darling test ;  # for Pareto and power distribution 
library('KScorrect'); library('plyr') ; library('data.table') #for loguniform function ; for count function ; for rbindlist list ; https://cran.r-project.org/web/views/Distributions.html 
}
Quantile <- function(a) {quant_list=c(min(a),quantile(a, na.rm=TRUE,type=1,probs=c(.005,.01,.05,.1,.25,.5,.75,.9,.95,.99,.995)), max(a))} ##quantiles for output

close(file('quant.txt', open="w" )) #output catalogue cleaning

mag_limit=9
mag_limit2=8
mag_limit3=7
mag_limit4=6

wds = read.table('sample9_clean.txt', sep='\t') ; 
#wds = read.table('new_clean.txt', sep=',') ; 
wds=wds[,2:5]

wds=wds[wds[,3]<3,]
#write.table(wds, file='separation.txt', sep='\t', row.names=F)


wds=wds[wds[,1]<mag_limit,]
wds=wds[wds[,2]<mag_limit,]
mag1=ifelse(wds[,1]<wds[,2], wds[,1],wds[,2]); mag2=ifelse(wds[,1]<wds[,2], wds[,2],wds[,1]); dv_obs=wds[,3]; sep_obs=wds[,4]

model=read.table('final+.txt', sep='\t')
colnames(model)=c('m1',"m2","q","age",'metal',"V1","V2","delta",'ev1','ev2',"dist","z","ext","mag_1","mag_2","axis","sep",'th_d', 'x', 'y', 'phi', 'l', 'g12g', 'g12e')

model=model[model$delta<3,] #; model=model[model$mag_2<mag_limit,]
#model= model[sample(nrow(model), 13000), ]
write.table(model$delta, file='delta+.txt', sep='\t', row.names=F)

dv_model=model$delta ; sep_model=model$sep ; #mag1_model=model$mag_1 ; mag2_model=model$mag_2
mag1_model=ifelse(model$mag_1<=model$mag_2, model$mag_1,model$mag_2)
mag2_model=ifelse(model$mag_1<=model$mag_2, model$mag_2,model$mag_1)

ks=ks.test(dv_obs,dv_model)$statistic[[1]] ;ks_test=ks.test(dv_obs, dv_model)$p.value ; 
ad_test_1=ad.test(dv_obs, dv_model)$ad[5] ; ad_test_2=ad.test(dv_obs, dv_model)$ad[6] 
ks_sep=ks.test(sep_obs, sep_model)$p.value; ad1_sep=ad.test(sep_obs, sep_model)$ad[5] ; ad2_sep=ad.test(sep_obs, sep_model)$ad[6] 
ks_mag1=ks.test(mag1, mag1_model)$p.value ; ad1_mag1=ad.test(mag1, mag1_model)$ad[5] ; ad2_mag1=ad.test(mag1, mag1_model)$ad[6]
ks_mag2=ks.test(mag2, mag2_model)$p.value ; ad1_mag2=ad.test(mag2, mag2_model)$ad[5] ; ad2_mag2=ad.test(mag2, mag2_model)$ad[6]

a=round(length(subset(model, model$mag_1<mag_limit4 & model$mag_2<mag_limit4)$delta)/58,2) ##60
b=round(length(subset(model, model$mag_1<mag_limit3 & model$mag_2<mag_limit3)$delta)/168,2) ##173
c=round(length(subset(model, model$mag_1<mag_limit2 & model$mag_2<mag_limit2)$delta)/460,2) ##478
d=round(length(dv_model)/length(mag1),2)

cat('N9/N8',length(dv_model)/length(subset(model, model$mag_1<mag_limit2 & model$mag_2<mag_limit2)$delta),'\n')
cat('mag_limit', mag_limit,'\n')
cat('Length, data,model',length(mag1),d,'\n') # Length of model sample and ks.test in terminal

cat('mag 6,7,8: ',a,b,c,'\n') 

cat('delta magnitude',ks_test, ad_test_1, ad_test_2,'\n')
cat('Primary mag',ks_mag1, ad1_mag1, ad2_mag1,'\n')
cat('Secondary mag',ks_mag2, ad1_mag2, ad2_mag2,'\n')
cat('separation', ks_sep, ad1_sep, ad2_sep,'\n')
#percent= ecdf(model$dist); cat('% at dist <100, 50 and 20 pc, median distance', percent(100),percent(50), percent(20), median(model$dist),'\n')
#length(subset(model, model$dist<400)$dist)/length(dv_model)
e=median(model$dist)
f=median(subset(model, model$mag_1<mag_limit3 & model$mag_2<mag_limit3)$dist)
i=median(subset(model, model$mag_1<mag_limit2 & model$mag_2<mag_limit2)$dist)
cat('median distance ', e,f,i,'\n')

g=quantile(model$dist, 0.25)
h=quantile(model$dist, 0.75)
summary(model$mag_1)
wilcox=wilcox.test(dv_obs,dv_model)$p.value
wil1=wilcox.test(mag1,mag1_model)$p.value
wil2=wilcox.test(mag2,mag2_model)$p.value
wils=wilcox.test(sep_obs,sep_model)$p.value

quant_list=signif(round(apply(model, 2, Quantile),3),3) #quantiles output
colnames(quant_list)=c('\tm1',"m2","q","age",'metal',"V1","V2","delta",'ev1','ev2',"dist","z","ext","mag_1","mag_2","axis","sep",'th_d', 'x', 'y', 'phi', 'l', 'g12g', 'g12e')
write.table(quant_list, file='quant.txt', append=T, sep='\t')

close=model[model$sep<median(sep_model),]
wide=model[model$sep>median(sep_model),]

quant_list=signif(round(apply(close, 2, Quantile),3),3) #quantiles output
colnames(quant_list)=c('\tm1',"m2","q","age",'metal',"V1","V2","delta",'ev1','ev2',"dist","z","ext","mag_1","mag_2","axis","sep",'th_d', 'x', 'y', 'phi', 'l', 'g12g', 'g12e')
write.table(quant_list, file='quant.txt', append=T, sep='\t')


quant_list=signif(round(apply(wide, 2, Quantile),3),3) #quantiles output
colnames(quant_list)=c('\tm1',"m2","q","age",'metal',"V1","V2","delta",'ev1','ev2',"dist","z","ext","mag_1","mag_2","axis","sep",'th_d', 'x', 'y', 'phi', 'l', 'g12g', 'g12e')
write.table(quant_list, file='quant.txt', append=T, sep='\t')

stats=data.table(length(dv_model), a,b,c,d,g,e,h,f,i, round(d/c,2), round(d/b,2), round(d/a,2))
colnames(stats)=c('length','ratio6','ratio7','ratio8','ratio9', '.25','median d9', '.75','median d7','median d8','9/8','9/7','9/6')
write.table(stats, file='quant.txt', append=T, sep='\t', row.names=F)

test_dv=round(data.table(ks_test, ad_test_1, ad_test_2, wilcox),3)
colnames(test_dv)=c('ks_dv','ad1','ad2','wil')
write.table(test_dv, file='quant.txt', append=T, sep='\t', row.names=F)

test_m1=round(data.table(ks_mag1, ad1_mag1, ad2_mag1, wil1),3)
colnames(test_m1)=c('ks_m1','ad1','ad2','wil')
write.table(test_m1, file='quant.txt', append=T, sep='\t', row.names=F)

test_m2=round(data.table(ks_mag2, ad1_mag2, ad2_mag2, wil2),3)
colnames(test_m2)=c('ks_m2','ad1','ad2','wil')
write.table(test_m2, file='quant.txt', append=T, sep='\t', row.names=F)

test_s=round(data.table(ks_sep, ad1_sep, ad2_sep,wils),3)
colnames(test_m2)=c('ks_sep','ad1','ad2','wil')
write.table(test_s, file='quant.txt', append=T, sep='\t', row.names=F)


postscript(file='deltaMag.ps'); #Delta mag plot of model and observation
#title=paste(scenario, m_split_1, imf_1, power_1, power_2, power_3, length(dv_model), round(ks_test,3)) #title for plot 
p=hist(dv_obs, plot=FALSE, breaks=seq(0,15,.5)) ; plot(p$mids, p$density, xlim=c(0,max(dv_obs)), ylim=c(0,1),pch=23, col='blue', cex=2) 
par(new=TRUE); p=hist(dv_model, plot=FALSE, breaks=seq(0,15,.5)) ; plot(p$mids, p$density, xlim=c(0,max(dv_obs)), ylim=c(0,1), main='title') ; dev.off() 

postscript(file='Mag1.ps'); #Delta mag plot of model and observation
#title=paste(scenario, m_split_1, imf_1, power_1, power_2, power_3, length(dv_model), round(ks_test,3)) #title for plot 
p=hist(mag1, plot=FALSE, breaks=seq(-6,9,.5)); plot(p$mids, p$density, xlim=c(-6,9), ylim=c(0,0.6),pch=23, col='blue', cex=2)
par(new=TRUE); p=hist(mag1_model, plot=FALSE, breaks=seq(-6,9,.5)); plot(p$mids, p$density, xlim=c(-6,9), ylim=c(0,0.6),main='title') ; dev.off() 

postscript(file='Mag2.ps'); #Delta mag plot of model and observation
#title=paste(scenario, m_split_1, imf_1, power_1, power_2, power_3, length(dv_model), round(ks_test,3)) #title for plot 
p=hist(mag2, plot=FALSE, breaks=seq(-6,9,.5)); plot(p$mids, p$density, xlim=c(-6,9), ylim=c(0,0.6),pch=23, col='blue', cex=2)
par(new=TRUE); p=hist(mag2_model, plot=FALSE, breaks=seq(-6,9,.5)); plot(p$mids, p$density, xlim=c(-6,9), ylim=c(0,0.6),main='title') ; dev.off() 

postscript(file='sep.ps'); #Delta mag plot of model and observation
#title=paste(scenario, m_split_1, imf_1, power_1, power_2, power_3, length(dv_model), round(mass_dens_2_total[[2]],5), round(ks_sep,3))
p=hist(sep_obs, plot=FALSE, breaks=10**(seq(-1,2,.1))); plot(p$mids, p$counts/length(sep_obs), xlim=c(.8,15), ylim=c(0,.2),pch=23, col='blue', cex=2, log='x')
par(new=TRUE); p=hist(sep_model, plot=FALSE, breaks=10**(seq(-1,2,.1))); plot(p$mids, p$counts/length(sep_model), xlim=c(.8,15), ylim=c(0,.2),main='title', log='x') ; dev.off()

stellar_density=read.table('exit_table+.txt', sep='\t')
density= data.table(sum(stellar_density[1]), sum(stellar_density[2]), sum(stellar_density[3]), sum(stellar_density[4]),sum(stellar_density[5]),sum(stellar_density[6]),sum(stellar_density[7]),sum(stellar_density[8]))
colnames(density)=c('0.72_p','0.72_s','1_p','1_s','1.6_p','1.6_s','2_p','2_s')
write.table(density, file='quant.txt', append=T, sep='\t', row.names=F)

#file.rename('quant.txt','quant+.txt')


