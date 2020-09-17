{ #libraries and options
options("width"=200) ; options(scipen=99) ; library('parallel')#Wide screen ; no scientific output
library('kSamples') ; library('extraDistr')  ## for Anderson–Darling test ;  # for Pareto and power distribution 
library('KScorrect'); library('plyr') ; library('data.table') #for loguniform function ; for count function ; for rbindlist list ; https://cran.r-project.org/web/views/Distributions.html 
}
Quantile <- function(a) {quant_list=c(min(a),quantile(a, na.rm=TRUE,type=1,probs=c(.005,.01,.05,.1,.25,.5,.75,.9,.95,.99,.995)), max(a))} ##quantiles for output

close(file('quant.txt', open="w" )) #output catalogue cleaning

tests='NO'
mag_limit=9
mag_limit8=8
mag_limit7=7
mag_limit6=6

model=read.table('final+.txt', sep='\t')
blue=read.table('tycho3.txt', sep=',')[,1]
color=read.table('tycho3.txt', sep=',')[,3]
mag1=read.table('tycho3.txt', sep=',')[,2]
#model=model[1:120,]
colnames(model)=c('m1',"age", 'metal',"V1", 'B1','ev1',"dist","z","ext","mag_1",'blue','color','th_d', 'x', 'y', 'phi', 'l', 'g12g', 'g12e')
count(model$th_d)

model=model[model$mag_1<mag_limit,] 

mag1_model=model$mag_1; blue_model=model$blue; color_model=model$color

length9=length(model$m1)
length8=length(subset(model, model$mag_1<mag_limit8)$m1)
length7=length(subset(model, model$mag_1<mag_limit7)$m1)
length6=length(subset(model, model$mag_1<mag_limit6)$m1)
length_pole=length(subset(model,abs(model$phi)>1.0472)$dist)
length_equator=length(subset(model,abs(model$phi)<0.08727)$dist)

cat('length, 6mag', length6,'\n')
cat('length, 9mag', length9,'\n')
cat('length, pole', length_pole,'\n')
cat('9/8, 8/7, 7/6:',length9/length8, length8/length7, length7/length6,'\n')

#ad.test(mag1, mag1_model)
#cat('mag_limit', mag_limit,'\n')
#cat('mag 6,7,8,9: ',length(subset(model, model$mag_1<mag_limit6)$m1),length(subset(model, model$mag_1<mag_limit7)$m1),length(subset(model, model$mag_1<mag_limit8)$m1),length(model$m1),'\n') 


dist9=median(model$dist)
dist7=median(subset(model, model$mag_1<mag_limit7)$dist)
dist79=median(subset(model, model$mag_1>mag_limit7)$dist)
q79=quantile(subset(model, model$mag_1>mag_limit7)$dist, c(0.25,0.5,0.75))
equ79=quantile(subset(model, model$mag_1>mag_limit7 & abs(model$phi)<0.08727)$dist, c(0.25,0.5,0.75))

pole79=quantile(subset(model, model$mag_1>mag_limit7 & abs(model$phi)>1.0472)$dist, c(0.25,0.5,0.75))

equator_color=quantile(subset(model,abs(model$phi)<0.08727)$color, c(0.25,0.5,0.75))
pole_color=quantile(subset(model,abs(model$phi)>1.0472)$color, c(0.25,0.5,0.75))

#cat('median distance 9', dist9,'\n')
#cat('median distance 7', dist7,'\n')
summary(subset(model, model$mag_1>mag_limit7))
cat('median distance 7-9', q79,'\n')
cat('pole median distance 7-9', pole79,'\n')




stellar_density=read.table('exit_table.txt', sep='\t')
cat(sum(stellar_density[1]), sum(stellar_density[2]), sum(stellar_density[3]), sum(stellar_density[4]),'\n')
density= data.table(sum(stellar_density[1]), sum(stellar_density[2]), sum(stellar_density[3]), sum(stellar_density[4]),sum(stellar_density[5]),sum(stellar_density[6]),sum(stellar_density[7]))
colnames(density)=c('0.72','1','1.6','2', 'full', 'local','thick')


if (tests=='YES')
{
ks_mag1=ks.test(mag1, mag1_model)$p.value ; ad1_mag1=ad.test(mag1, mag1_model)$ad[5] ; ad2_mag1=ad.test(mag1, mag1_model)$ad[6]
#ks_blue=ks.test(blue, blue_model)$p.value ; ad1_blue=ad.test(blue, blue_model)$ad[5] ; ad2_blue=ad.test(blue, blue_model)$ad[6]
#ks_color=ks.test(color, color_model)$p.value ; ad1_color=ad.test(color, color_model)$ad[5] ; ad2_color=ad.test(color, color_model)$ad[6]
ks_blue=2; ad1_blue=2; ad2_blue=2; ks_color=2; ad1_color=2; ad2_color=2
}

if (tests=='NO')
{ks_mag1=2 ; ad1_mag1=2 ; ad2_mag1=2; ks_blue=2; ad1_blue=2; ad2_blue=2; ks_color=2; ad1_color=2; ad2_color=2}
quant_list=signif(round(apply(model, 2, Quantile),3),3) #quantiles output
colnames(quant_list)=c('\tm1',"age",'metal',"V1",'B1','ev1',"dist","z","ext","mag_1",'blue','B-V','th_d', 'x', 'y', 'phi', 'l', 'g12g', 'g12e')
write.table(quant_list, file='quant.txt', append=T, sep='\t')


quant_list=signif(round(apply(subset(model,abs(model$phi)<0.08727), 2, Quantile),3),3) #quantiles output
colnames(quant_list)=c('\tm1',"age",'metal',"V1",'B1','ev1',"dist","z","ext","mag_1",'blue','B-V','th_d', 'x', 'y', 'phi', 'l', 'g12g', 'g12e')
write.table(quant_list, file='quant.txt', append=T, sep='\t')

quant_list=signif(round(apply(subset(model,abs(model$phi)>1.0472), 2, Quantile),3),3) #quantiles output
colnames(quant_list)=c('\tm1',"age",'metal',"V1",'B1','ev1',"dist","z","ext","mag_1",'blue','B-V','th_d', 'x', 'y', 'phi', 'l', 'g12g', 'g12e')
write.table(quant_list, file='quant.txt', append=T, sep='\t')

quant_list=signif(round(apply(subset(model,model$th_d==2), 2, Quantile),3),3) #quantiles output
colnames(quant_list)=c('\tm1',"age",'metal',"V1",'B1','ev1',"dist","z","ext","mag_1",'blue','B-V','th_d', 'x', 'y', 'phi', 'l', 'g12g', 'g12e')
write.table(quant_list, file='quant.txt', append=T, sep='\t')

stats=data.table(length6, length9, round(length9/length8,3), round(length8/length7,3), round(length7/length6,3), round(q79), round(equ79), round(pole79), round(length_pole/length9,3), round(length_equator/length9,3), round(equator_color,3), round(pole_color,3))
colnames(stats)=c('length6','length9','9/8 mag','8/7 mag','7/6 mag','dist 7-9', 'pole dist 7-9','equator distance 7-9', 'pole fr.','equator fr.', 'color eq.', 'color pole')
write.table(stats, file='quant.txt', append=T, sep='\t', row.names=F)

write.table(density, file='quant.txt', append=T, sep='\t', row.names=F)

tests=data.table(ks_mag1, ad1_mag1, ad2_mag1, ks_blue, ad1_blue, ad2_blue, ks_color, ad1_color, ad2_color)
tests=round(tests,3)
colnames(tests)= c('ks','ad1','ad2','blue','','','color','','')
write.table(tests, file='quant.txt', append=T, sep='\t', row.names=F)

summary(color)


postscript(file='Mag1.ps'); #Delta mag plot of model and observation
#title=paste(scenario, m_split_1, imf_1, power_1, power_2, power_3, length(dv_model), round(ks_test,3)) #title for plot 
p=hist(mag1, plot=FALSE, breaks=seq(-6,9,.5)); plot(p$mids, p$density, xlim=c(-6,9), ylim=c(0,0.6),pch=23, col='blue', cex=2)
par(new=TRUE); p=hist(mag1_model, plot=FALSE, breaks=seq(-6,9,.5)); plot(p$mids, p$density, xlim=c(-6,9), ylim=c(0,0.6),main='title') ; dev.off() 

postscript(file='blue.ps'); #Delta mag plot of model and observation
p=hist(blue, plot=FALSE, breaks=seq(-6,15,.5)); plot(p$mids, p$density, xlim=c(-6,15), ylim=c(0,0.6),pch=23, col='blue', cex=2)
par(new=TRUE); p=hist(blue_model, plot=FALSE, breaks=seq(-6,15,.5)); plot(p$mids, p$density, xlim=c(-6,15), ylim=c(0,0.6),main='title') ; dev.off() 

postscript(file='color.ps'); #Delta mag plot of model and observation
p=hist(color, plot=FALSE, breaks=seq(-1,7,.5)); plot(p$mids, p$density, xlim=c(-1,7), ylim=c(0,1),pch=23, col='blue', cex=2)
par(new=TRUE); p=hist(color_model, plot=FALSE, breaks=seq(-1,7,.5)); plot(p$mids, p$density, xlim=c(-1,7), ylim=c(0,1),main='title') ; dev.off() 


SELECT
list.source_id , gaiadr2_complements.geometric_distance.r_est
from
list
join
gaiadr2_complements.geometric_distance
on
list.source_id = gaiadr2_complements.geometric_distance.source_id

