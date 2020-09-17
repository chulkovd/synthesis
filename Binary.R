Sys.sleep(5)

clock=Sys.time()
{ #libraries and options
options("width"=200) ; options(scipen=99) ; library('parallel')#Wide screen ; no scientific output
library('kSamples') ; library('extraDistr')  ## for Anderson–Darling test ;  # for Pareto and power distribution 
library('KScorrect'); library('plyr') ; library('data.table') #for loguniform function ; for count function ; for rbindlist list ; https://cran.r-project.org/web/views/Distributions.html 
}

Quantile <- function(a) {quant_list=c(min(a),quantile(a, na.rm=TRUE,type=1,probs=c(.005,.01,.05,.1,.25,.5,.75,.9,.95,.99,.995)), max(a))} ##quantiles for output
scenario='RP' ## RP, PCP, SCP, single
if (scenario != 'single') {isochrones <- fread('track2.dat', head=TRUE, sep="|")}
if (scenario == 'single') {isochrones <- fread('2colors.dat', head=TRUE, sep="|")}
summary(isochrones)
summary(subset(isochrones, isochrones$ev==1))
isochrones=subset(isochrones,isochrones$ev<8.5)
close(file('final.txt', open="w" )) #output catalogue cleaning
close(file('exit_table.txt', open="w" )) #output catalogue cleaning
close(file('ages.txt', open="w" )) #output catalogue cleaning
close(file('metals.txt', open="w" )) #output catalogue cleaning

mass_dens_1_total=mass_dens_2_total=mass_dens_1=mass_dens_2=thick_length=neigh_length=nb_=th_d=generation=NULL
k1_q_pdf=k2_q_pdf=q1_length=q2_length=q3_length=0 ; num=1 ##initialization

m_min=.3; m_split_1=0.5; m_split_2=5; m_max=68; m_step=10**-4 #minimum and maximum mass, splits; step for mass pdf
imf_1='power' ; imf_mean=0.07; sigma_imf=0.7 ; mean_imf=log10(imf_mean); ## power or lognorm, # parameters for lognormal mass distribution (only small masses supported)  

power_1=-1.3; power_2=power_3=-2
RPtwins=T;twin_fr= 0.025; twin_sigma=0.01 ##parameters for power distribution: f(m) ~ m^power
unresolved=F ; mult=.8 # doesn't work for single scenario

cores=10 ; image=seq(cores) #number of runs in parallel
n=1000000; num_limit=7000  #  n=500000; num_limit=2000 single stars 
minimal_age=0.004 ; thick_age=10; maximum_age=13 ; sfr=0.1 ##4*10**6 < age < 1.35*10**10
thick_age=log10(thick_age)+9.005

max_dist=2000; Sun=8500 ; thick_disk_height=900 ; thin_disk_length=2500 ; thick_disk_length=3500 ##Galactic parameters 

#z0=50; t0=5*10**7 ; scale=1; z_distr='sech'; height='not fixed'; thin_disk_height=190
#z0=50; t0=5.55*10**8 ; scale=5/3; z_distr='sech'; ; 
thin_disk_height=200
height='not fixed'; z_distr='sech'
k_ext=1; k_height=1


mag_limit_1=mag_limit_2=9; delta_mag_limit=15 # Limiting apparent magnitude and magnitude difference
sep_min=0.8; sep_max=15 ; axis_min=5; axis_split=5000; axis_max=sep_max*max_dist; ax=-1.2 ; ap=-1.6 #Limits on axis (separation in AU) and separation in arcsec
local=50; mass_dens_list=c(.72,1,1.6,2) #local neighbourhood distance ; cuts for mass density - see https://arxiv.org/abs/1704.05063 , page 12

q_min=.1; q_split_1=.3 ; q_split_2=.95 ; q_step=.001 #q split positions; minimal q; step in q
small_q='power'; medium_q='power'; large_q='power' #choice of q distribution, power or gauss
beta1=beta2=-0.6; beta3=0; jump=1.4 #parameters for power distribution: f(q) ~ q ^ beta ;  twin star excessean1=.08; sigma1=0.7; mean2=.6; sigma2=0.6; mean3=1.6; sigma3=0.6 # parameters for Gaussian q distribution
generation_plot=F; plot=F #Plots info
xlim=c(0.1,10) ; ylim=c(-3.5,5.5) #  limits for approximation plot, masses vs abs mag
q_plot_number=100 ; q_plot_split=100 # plot every _th graphic of q ; split q in columns
m_plot_number=100 ; m_plot_split=10**seq(-1.1,1.9,.1); #log m plot
plot_number=1 #plot every _th graphic approximation

a0=0.0012; a1=0.0003; a2=55*pi/180 ; z_a=70; gamma=19*pi/180; a0g=0.0012; a1g=0.0011; a2g=130*pi/180; zeta_a=50 #extinction #alpha=0.0016 ; beta=100

if (scenario=='SCP') {m_max=2*m_max; m_min=2*m_min; m_split_1=2*m_split_1}
if (scenario =='RP' & RPtwins=='T') {m_max=2*m_max}
#check=ifelse(scenario=='single',T,F)
check=T

Testing_function = function(image) { ##Function starts

if (num_limit< 5 & image == 1 & scenario !='single') {cat('\n','axis distr',ax,ap)}

if (scenario=='PCP' | scenario=='SCP') { # Creating q pdf for PCP and SCP
c1_q=seq(0, q_split_1, q_step); c2_q=seq(q_split_1, q_split_2, q_step); c3_q=seq(q_split_2, 1, q_step) #small, medium, large q

if (small_q=='power') ###small q pdf
{if (beta1>(-1)) {pdf_q_1=dpower(c1_q,q_split_1, beta1+1) ; if (num_limit< 5 & image == 1) {cat('\nsmall q - power', beta1,' ; ')}}
if (beta1==(-1)) {pdf_q_1=dlunif(c1_q,q_min,q_split_1); if (num_limit< 5 & image == 1) {cat('\nsmall q -1 ;')}}
if (beta1<(-1)) {pdf_q_1=dpareto(c1_q,-beta1-1, q_min) ; if (num_limit< 5 & image == 1) {cat('\nsmall q - pareto', beta1,' ; ')}}}
if (small_q=='gauss') {pdf_q_1=dnorm(c1_q,mean1, sigma1) ; if (num_limit< 5 & image == 1) {cat('\nsmall q - gauss', mean1, sigma1,' ; ')}} #normal distribution

if (medium_q=='power') ###medium q
{if(beta2>(-1)) {pdf_q_2=dpower(c2_q,q_split_2, beta2+1) ; if (num_limit< 5 & image == 1) {cat('medium q - power', beta2,' ; ')}}
if (beta2==(-1)) {pdf_q_2=dlunif(c2_q,q_split_1,q_split_2); if (num_limit< 5 & image == 1) {cat('medium q -1 ;')}} 
if (beta2<(-1)) {pdf_q_2=dpareto(c2_q,-beta2-1, q_split_1) ; if (num_limit< 5 & image == 1) {cat('medium q - pareto', beta2,' ; ')}}}
if (medium_q=='gauss') {pdf_q_2=dnorm(c2_q, mean2, sigma2) ; if (num_limit< 5 & image == 1) {cat('medium q - gauss', mean2, sigma2, ' ; ')}} #normal distribution

if (large_q=='power') ###large q
{if (beta3>-1) {pdf_q_3=dpower(c3_q,1, beta3+1) ; if (num_limit< 5 & image == 1) {cat('large q - power', beta3,' ; ')}}
if (beta3==(-1)) {pdf_q_3=dlunif(c3_q,q_split_2,1); if (num_limit< 5 & image == 1) {cat('large q -1 ;')}}
if (beta3<(-1)) {pdf_q_3=dpareto(c3_q,-beta3-1, q_split_2) ; if (num_limit< 5 & image == 1) {cat('large q - pareto', beta3,' ; ')}}}
if (large_q=='gauss') {pdf_q_3=dnorm(c3_q, mean3, sigma3) ; if (num_limit< 5 & image == 1) {cat('large q - gauss', mean3, sigma3)}} #normal distribution
if (num_limit< 5 & image == 1) {cat('step', jump)}
k1_q_pdf=pdf_q_1[length(pdf_q_1)-1]/pdf_q_2[1] ; k2_q_pdf=pdf_q_2[length(pdf_q_2)-1]/pdf_q_3[1]  ###Adjustment of piecewise function ; bracket for non-RP if
}

{# Creating mass probability distribution function
m1=seq(m_min,m_split_1,m_step) ; m2=seq(m_split_1,m_split_2,m_step) ; m3=seq(m_split_2,m_max,m_step) # Small, medium, large masses 

if (imf_1=='power') #small m pdf
{if (power_1<(-1)) {d1=dpareto(m1, -power_1-1,m_min) ; if (num_limit< 5 & image == 1) {cat('\nsmall m - pareto, split', power_1, m_split_1,' ; ')}}
if (power_1==(-1)) {d1=dlunif(m1, m_min, m_split_1) ; if (num_limit< 5 & image == 1) {cat('\nsmall m -1 ;')}}
if (power_1>(-1)) {d1=dpower(m1, m_split_1,power_1+1) ; if (num_limit< 5 & image == 1) {cat('\nsmall m, power', power_1,' ; ')}}}
if (imf_1=='lognorm') {d1=dlnorm(m1,mean_imf, sigma_imf) ; if (num_limit< 5 & image == 1) {cat('\nsmall m - gauss', mean_imf, sigma_imf,' ; ')}} #lognormal distribution

if (power_2<(-1)) {d2=dpareto(m2, -power_2-1,m_split_1) ; if (num_limit< 5 & image == 1) {cat('medium m, pareto', power_2,' ; ')}} #medium m pdf
if (power_2==(-1)) {d2=dlunif(m2, m_split_1, m_split_2) ; if (num_limit< 5 & image == 1) {cat('medium m -1 ;')}}
if (power_2>(-1)) {d2=dpower(m2, m_split_2,power_2+1) ; if (num_limit< 5 & image == 1) {cat('medium m, power', power_2, ' ; ')}}

if (power_3<(-1)) {d3=dpareto(m3, -power_3-1,m_split_2) ; if (num_limit< 5 & image == 1) {cat('large m, pareto' , power_3)}} #large m pdf
if (power_3==(-1)) {d3=dlunif(m3, m_split_2, m_max) ; if (num_limit< 5 & image == 1) {cat('large m -1 ;')}} 
if (power_3>(-1)) {d3=dpower(m3, m_max,power_3+1) ; if (num_limit< 5 & image == 1) {cat('large m, power', power_3)}}

k1=d1[length(d1)-1]/d2[1] ; k2=d2[length(d2)-1]/d3[1] ; #Adjustment of piecewise function for imf
}

if (num_limit< 5 & image == 1) {cat('\nk for q and m pdf:', k1_q_pdf, k2_q_pdf, k1,k2, '\n')} ##terminal output, m and q pdf

if (sfr==0) 
{ages=round(log10(runif(num_limit, minimal_age*10**9, maximum_age*10**9)),2) #Uniform age distribution

}
#if (image==1) {cat('ages',summary(ages),'\n')}}

if (sfr>0)
{
flat=runif(num_limit, 0,1)
ages=10**9*(1/sfr)*log(flat*(exp(sfr*maximum_age)-exp(sfr*minimal_age))+exp(sfr*minimal_age))

#ages=rpower(num_limit, maximum_age-minimal_age, sfr+1) #decresing sfr

ages=round(log10(minimal_age+ages),2)
if (image==1) {cat('ages',summary(ages),'\n')}}

count_ages=count(ages); write.table(count_ages, file='ages.txt', append=T, col.names=FALSE, row.names=FALSE) ##Age statistics, for output

metals1=round(rnorm(num_limit, 0, 0.1),1) ; metals1=ifelse(metals1<(-0.45), -0.4,metals1) ; metals1=ifelse(metals1> 0.45, 0.4,metals1)
metals2=round(rnorm(num_limit, -0.1, 0.1),1) ; metals2=ifelse(metals2<(-.55), -0.5,metals2) ; metals2=ifelse(metals2>0.35, 0.3,metals2)
metals3=round(rnorm(num_limit, -0.1, 0.15),1) ; metals3=ifelse(metals3<(-.75), -0.7,metals3) ; metals3=ifelse(metals3>0.55, 0.5,metals3)
metals4=round(rnorm(num_limit, -0.2, 0.2),1) ;  metals4=ifelse(metals4<(-1.05), -1.0,metals4); metals4=ifelse(metals4>0.65, 0.6,metals4)
metals5=round(rnorm(num_limit, -0.5, 0.3),1) ; metals5=ifelse(metals2<(-1.65), -1.6,metals4); metals5=ifelse(metals4>0.65, 0.6,metals4)
metals=ifelse(ages<9.305, metals1, -10) #< 2 billion
metals=ifelse(ages>9.305 & ages<9.605, metals2,metals)# 2-4 billion
metals=ifelse(ages>9.605 & ages<9.925, metals3,metals)#4-8.5 billion
metals=ifelse(ages>9.925 & ages<10.005, metals4,metals) #8.5-10 billion 
metals=ifelse(ages>10.005, metals5, metals) #>10 billion

#count_metal=count(metals); write.table(count_metal, file='metal.txt', append=T, col.names=FALSE, row.names=FALSE) ##Age statistics, for output

##random generation
repeat { #outer cycle start (num is the variable, num_limit)

time=ages[num]
feh=metals[num]

system_list=rho=z=NULL; na_count=0  

if (image == 1) {cat('Outer cycle iteration №', num,'\n')}


if (scenario=='PCP' | scenario=='SCP') { # q generation start for PCP and SCP

if (small_q=='power') #small q generation
{if (beta1>(-1)) {q1=rpower(n, q_split_1, beta1+1)}
if(beta1==(-1)) {q1=rlunif(n,q_min,q_split_1)}
if (beta1<(-1)) {q1=rpareto(n,-beta1-1, q_min)}}
if (small_q=='gauss') {q1=rnorm(n, mean1, sigma1)}
q1=subset(q1, q1>=q_min & q1<=q_split_1)

if (medium_q=='power') #medium q generation
{if (beta2>(-1)) {q2=rpower(n*k1_q_pdf, q_split_2, beta2+1)}
if(beta2==(-1)) {q2=rlunif(n*k1_q_pdf,q_split_1,q_split_2)}
if (beta2<(-1)) {q2=rpareto(n*k1_q_pdf,-beta2-1, q_split_1)}}
if (medium_q=='gauss') {q2=rnorm(n*k1_q_pdf, mean2, sigma2)}

q2=subset(q2, q2>q_split_1 & q2<=q_split_2)

if (large_q=='power') # large q generation
{if (beta3>(-1)) {q3=rpower(n*k1_q_pdf*k2_q_pdf*jump, 1, beta3+1)}
if(beta3==(-1)) {q3=rlunif(n*k1_q_pdf*k2_q_pdf*jump,q_split_2,1)}
if (beta3<(-1)) {q3=rpareto(n*k1_q_pdf*k2_q_pdf*jump,-beta3-1, q_split_2)}}
if (large_q=='gauss') {q3=rnorm(n*k1_q_pdf*k2_q_pdf*jump, mean3, sigma3)}
q3=subset(q3, q3>q_split_2 & q3<=1)
q=c(q1,q2,q3)} ##q merge



if (imf_1=='power') ###mass generation , all scenarios 
{if (power_1<(-1)) {mass1=rpareto(n, -power_1-1, m_min)} #low mass generation
if (power_1==(-1)) {mass1=rlunif(n,m_min,m_split_1)}
if (power_1>(-1)) {mass1=rpower(n, m_split_1, power_1+1)}}
if (imf_1=='lognorm') {mass1=rlnorm(n, mean_imf, sigma_imf)}
mass1=subset(mass1, mass1 >= m_min & mass1 <= m_split_1)

if (power_2<(-1)) {mass2=rpareto(n*k1, -power_2-1, m_split_1)} #medium mass generation
if (power_2==(-1)) {mass2=rlunif(n*k1,m_split_1,m_split_2)}
if (power_2>(-1)) {mass2=rpower(n*k1, m_split_2, power_2+1)}
mass2=subset(mass2, mass2 > m_split_1 & mass2 <= m_split_2)

if (power_3<(-1)) {mass3=rpareto(n*k1*k2, -power_3-1, m_split_2)} #high mass generation
if (power_3==(-1)) {mass3=rlunif(n*k1*k2,m_split_2,m_max)}
if (power_3>(-1)) {mass3=rpower(n*k1*k2, m_max, power_3+1)}
mass3=subset(mass3, mass3 > m_split_2 & mass3 <= m_max)
mass=c(mass1,mass2,mass3)


 #pairing functions
###Omitting of excess elements (number of m should equal number of q), also randomization
if (scenario=='PCP' | scenario=='SCP') {number=min(length(mass),length(q)) ; q=sample(q, number) ; mass=sample(mass, number)}
if (scenario=='RP') {mass=sample(mass, 6*length(mass)%/%6); number=length(mass)/2} #length of masses array should be even for further division in two

#Pairing mechanisms
if (scenario=='single') {mass_1=mass}
if (scenario=='PCP') {mass_1=mass ; mass_2=mass*q} #PCP: m1, q
if (scenario=='SCP') {mass_1=mass/(1+q) ; mass_2=mass-mass_1} #SCP: m1+m2, q
if (scenario=='RP') 
{
if (RPtwins==F) {half=split(mass, sample(rep(1:2, length(mass)/2))) ; mass_1=half[[1]] ; mass_2=half[[2]]}
if (RPtwins==T) {half=split(mass, sample(rep(1:3, length(mass)/3))) ; mass_1=half[[1]] ; mass_2=half[[2]]; mass_twin=half[[3]]}
} # RP: m1, m2; Divide masses in two parts


if (scenario=='RP') 
{if (RPtwins==T) {twn= floor(twin_fr*length(mass_twin)); q=0.5+abs(rnorm(twn, 0,twin_sigma))} #RP+twins

if (RPtwins==T) {mass_twin=sample(mass_twin, twn); mass_twin1=mass_twin*q ; mass_twin2=mass_twin-mass_twin1 ; mass_1=c(mass_1, mass_twin1); mass_2=c(mass_2, mass_twin2)}
}
if (scenario=='SCP') {mass_2[2*mass_2<m_min]=NA} # Secondary mass limit for every pairing mechanism
if (scenario=='RP') {mass_2[mass_2<m_min]=NA} # Secondary mass limit for every pairing mechanism
if (scenario=='PCP') {mass_2[mass_2<m_min]=NA} # Secondary mass limit for every pairing mechanism

if (scenario!='single')
{generation=data.table(mass_1,mass_2)
generation=generation[!is.na(generation$mass_1)] 
generation=generation[!is.na(generation$mass_2)] 
mass_1=generation$mass_1
mass_2=generation$mass_2}

one_age=isochrones[age==time & metal==feh]
m_max_iso = max(one_age$mass) ; m_min_iso = min(one_age$mass) ; 
mass_1[mass_1>m_max_iso]=NA ;  if (scenario!='single') {mass_2[mass_2>m_max_iso]=NA} #Omit masses larger than maximum mass on the isochrone
mass_2[mass_2<m_min_iso]=NA ;  if (scenario!='single') {mass_2[mass_2>m_max_iso]=NA} #Omit masses larger than maximum mass on the isochrone

#{write(c(image, num, time, feh),append=T,'error2.txt')}
f=approxfun(one_age[,mass], one_age[,V], method='linear')
fev=approxfun(one_age[,mass],one_age[,ev], method='linear')
if (scenario=='single') {fblue=approxfun(one_age[,mass], one_age[,B], method='linear')}
#{write(c(image, num, time, feh, f(0.5), f(0.6)),append=T,'error.txt')}

if (scenario!='single') {star1=ifelse(mass_1>=mass_2, mass_1, mass_2) ;  star2=ifelse(mass_1>=mass_2, mass_2, mass_1)} #mass1 > mass2 condition for binary
if (scenario=='single') {star1=mass_1}

if (unresolved==T & scenario!='single')
{multipleA=runif(length(star1)); multipleB=runif(length(star1)) 
star1=ifelse(multipleA>mult,star1/2,star1); 
mag_abs_1=ifelse(multipleA>mult,f(star1)-0.75,f(star1))
star2=ifelse(multipleB>mult,star2/2,star2); 
mag_abs_2=ifelse(multipleB>mult,f(star2)-0.75,f(star2))}

if (unresolved==T & scenario=='single')
{multipleA=runif(length(star1)); star1=ifelse(multipleA>mult,star1,star1); mag_abs_1=ifelse(multipleA>mult,f(star1)-0.75,f(star1))}

if (unresolved!=T) {mag_abs_1=f(star1) ; if (scenario!='single') {mag_abs_2=f(star2)} } 

ev1=round(fev(star1),2) ; if (scenario!='single'){ev2=round(fev(star2),2)}  #absolute mag and luminosity class for a given mass
if (scenario=='single') {blue_abs=round(fblue(star1),2)}

if (scenario!='single') {system=data.table(star1, star2, time, feh, mag_abs_1, mag_abs_2, ev1, ev2, num)} # List of systems for a given age
if (scenario=='single') {system=data.table(star1, time, feh, mag_abs_1, ev1, num, blue_abs)} # List of systems for a given age
na_count=na_count+sum(!complete.cases(system)) ; system=system[complete.cases(system),] ##count and omit NA masses

if (plot==T) {if (num %% plot_number == 0) ### plot approximation of magnitude vs mass
{png(filename=paste('approx_',num, '.png')) #mass generation plot
plot(star1, mag_abs_1, xlim=xlim, ylim=ylim, xlab='mass', ylab='Absolute magnitude', col=2+ev1) ### Check approximation plot for a given age 
par(new=TRUE) ; curve(f, main=time, xlim=xlim, ylim=ylim, xlab=NULL, ylab=NULL); dev.off()}} 

system_list[[num]] = system ; #Create list of systems for various ages, proceed to the next age , End of time cycle

if (scenario == 'PCP' | scenario == 'SCP') {q1_length=length(q1) ; q2_length=length(q2); q3_length=length(q3)} ## for RP number preset to 0; output control only
if (num_limit< 5 & image == 1 & scenario!='single') {cat('Length of last generated pieces of q and m',q1_length, q2_length, q3_length, length(mass1), length(mass2), length(mass3) , number,'\n')}

system_list=rbindlist(system_list) #creating system_list, make use of 'data.table' library
if (scenario!='single') {mass_1=system_list[[1]] ; mass_2=system_list[[2]] ; age=system_list[[3]] ; metal=system_list[[4]] ; mag_abs_1=system_list[[5]] ; mag_abs_2=system_list[[6]] ; ev1=system_list[[7]] ; ev2=system_list[[8]]} #Systems
if (scenario=='single') {mass_1=system_list[[1]] ; age=system_list[[2]] ; metal=system_list[[3]] ; mag_abs_1=system_list[[4]] ; ev1=system_list[[5]] ; blue_abs=system_list[[7]] } #Systems

if (num_limit< 5 & image == 1) {cat('system_list size, Mb',object.size(system_list)/1024**2,'\n')}
nb=length(mass_1); if (num_limit< 5 & image == 1) {cat('Number of binaries with generated mass, number of NA',nb, na_count,'\n')}  #number of binaries for further calculations

z2=rexp(1.02*nb/pexp(max_dist, 1/thick_disk_height), 1/thick_disk_height) ; 
z2=sample(z2[z2<max_dist],nb);  #Thick disk z generation

if (height=='not fixed')
{
#z_h=z0+((10**age)/t0)**scale
#z_h=z0*(1+(10**age)/t0)**scale

#z_h=(1+(10**age)/t0)**(2/3) * z0
#if (num_limit< 5 & image==1) {cat('summary z_h:', summary(z_h),'\n')}


z_h=ifelse(age<=8.7, (35**2+(175*10**(age-9))**2)**0.5, 0) 
z_h=ifelse(age>8.7, 78*(1+10**(age-9))**0.5,z_h)

#z_h=ifelse(age<=8.7, (45**2+(417*10**(age-9))**2)**0.5, 0) #original form
#z_h=ifelse(age>8.7, 177*(1+10**(age-9))**0.5,z_h)
#z_h=z_h*k_height


#if (num_limit< 5 & image==1) {cat('summary z_h:', summary(z_h),'\n')}
if (z_distr=='sech')
{uniform=runif(1.05*nb) ; z_age=abs(log(uniform) - log(1-uniform)) *z_h}

if (z_distr=='exp')
{z_age=-log(runif(1.05*nb))*z_h}

z1=sample(z_age[z_age<max_dist],nb)}

if (height=='fixed' & z_distr=='exp') {z1=rexp(1.02*nb/pexp(max_dist, 1/thin_disk_height), 1/thin_disk_height)
z1=sample(z1[z1<max_dist],nb)} 

if (height=='fixed' & z_distr=='sech') 
{uniform=runif(1.05*nb) ; z1=abs(log(uniform) - log(1-uniform)) * thin_disk_height
z1=sample(z1[z1<max_dist],nb)}

if (num_limit< 5 & image==1) {cat('summary z thin disk:', summary(z1),'\n')}

r=seq(Sun-max_dist,Sun+max_dist,0.1) #radial exponential distribution, value 
u1=1-(1+r/thin_disk_length)*exp(-r/thin_disk_length) ; u2=1-(1+r/thick_disk_length)*exp(-r/thick_disk_length) #argument
radial1<-approxfun(u1,r) ; radial2<-approxfun(u2,r) ## approximation function 
test1=runif(nb, min(u1), max(u1)) ; test2=runif(nb, min(u2), max(u2)) ##random generation of argument
rho1=radial1(test1) ; rho2=radial2(test2); #distance from galactic centre

teta=runif(nb,-asin(max_dist/Sun),asin(max_dist/Sun)) #galactocentric angle uniform distribution
th_d=ifelse(age<thick_age,1,2)
rho=ifelse(th_d==1,rho1,rho2); z=ifelse(th_d==1, z1,z2); sign=sample(c(-1,1),nb,replace=T)

dist=(rho**2+Sun**2-2*rho*Sun*cos(teta)+z**2)**.5 ; dist[dist>max_dist]=NA ; 
if (num_limit< 5 & image == 1) {cat('distance', quantile(dist, na.rm=TRUE),'\n')} #distance from the Sun
module=5*log10(dist)-5 ; phi=asin(z/dist)*sign ; #ext=(1-exp(-dist*sin(phi)/beta))*alpha*beta/sin(phi) #module of distance, Extinction barometric law
x=rho*sin(teta); y=rho*cos(teta)-Sun; l=atan2(x,-y); l=ifelse(l<0,l+2*pi,l)

### extinction
gould_b=asin( (sin(phi)) * (cos(gamma)) - (sin(gamma)) * (cos(phi)) * (sin(l+pi/2)) ) #identical to eq => az, but l=l+90
gould_l1=asin((cos(phi)*sin(l+pi/2)*cos(gamma)+sin(phi)*sin(gamma))/cos(gould_b))
gould_l2 =acos(cos(phi)*cos(l+pi/2)/cos(gould_b))

gould_l=ifelse(gould_l1<0, -gould_l2-pi/2, NA)
gould_l=ifelse(gould_l1>0, gould_l2-pi/2, gould_l)
gould_l=ifelse(gould_l<0,gould_l+2*pi, gould_l)

zeta=dist*sin(gould_b); 
g12e=(1-exp(-z/z_a))*(z_a/z)*(a0+a1*sin(l+a2))*dist; 
distg=ifelse(dist>600,600,dist)
g12g=(1-exp(-abs(zeta)/zeta_a))*(zeta_a/abs(zeta))*(a0g+a1g*sin(2*gould_l+a2g))*distg
ext=k_ext*(g12e+g12g) #ext=(1-exp(-abs(z)/beta))*alpha*beta*dist/z


if (num_limit< 5 & image==1) {cat('ext',Quantile(ext),'\n')}

if (scenario!='single') { ### axis and separation
#axis=runif(nb, log10(axis_min), log10(axis_max)) ; axis=10**axis ; cat('axis, AU', Quantile(axis),'\n') # Uniform separation distribution in log a

c1=seq(axis_min,axis_split) ; c2=seq(axis_split, axis_max)
if (ax==(-1)) {pdf_a1=dlunif(c1,axis_min,axis_split)} 
if (ax<(-1)) {pdf_a1=dpareto(c1, -ax-1, axis_min)} 
pdf_a2=dpareto(c2, -ap-1, axis_split) 

k=pdf_a1[length(pdf_a1)]/pdf_a2[1]
if (ax==(-1)) {axis1=rlunif(nb,axis_min,axis_split)}
if (ax<(-1)) {axis1=rpareto(2.1*nb, -ax-1,axis_min) }
axis1=axis1[axis1<axis_split] ; axis1=sample(axis1,nb)

axis2=rpareto(nb*k,-ap-1, axis_split); 
axis2=axis2[axis2<axis_max] ; axis=sample(c(axis1,axis2),nb)


if (num_limit< 5 & image == 1) {cat('axis, AU', Quantile(axis),'\n')} # Separation distribution
sep=axis/dist ; sep[sep>sep_max] = NA ; sep[sep<sep_min] = NA ; if (num_limit< 5 & image == 1) {cat('separation, arcsec', quantile(sep, na.rm=TRUE),'\n')}#Calculation and limits on separation
} # not single bracket

app_mag_1=mag_abs_1+module+ext; if (scenario!='single') {app_mag_2=mag_abs_2+module+ext}; # Apparent mag from distance, absolute mag and extinction
if (scenario=='single') {blue_app=blue_abs+module+ext+ext/3.1}

#app_mag_1[app_mag_1>mag_limit_1]=NA ; app_mag_2[app_mag_2>mag_limit_2]=NA 
if (scenario!='single') {delta_mag=abs(app_mag_1-app_mag_2)} ; #delta_mag[delta_mag>delta_mag_limit]=NA #Imposing limits on delta-m in the model

if (scenario!='single') {position=data.table(mass_1, mass_2, mass_2/mass_1, age, metal, mag_abs_1, mag_abs_2, delta_mag, ev1, ev2, dist, z, ext, app_mag_1, app_mag_2, axis, sep, th_d, x, y, phi, l, g12g, g12e)} #Catalogue
if (scenario=='single') {position=data.table(mass_1, age, metal, mag_abs_1, blue_abs, ev1, dist, z, ext, app_mag_1, blue_app,blue_app-app_mag_1, th_d, x, y, phi, l, g12g, g12e)}
#fwrite(position, file='position0.txt', sep='\t', na='***') #catalogue

if (image == 1) {cat('position_list size, Mb',object.size(position)/1024**2,'\n')}
position=position[!is.na(position$dist)] #distances contained NA

#some plots
if (num==1 & scenario!='single' & plot==T)
{
png('galaxy1.png' ); smoothScatter(position$x,position$y);dev.off()
png('galaxy2.png' ); smoothScatter(position$x,position$z);dev.off()
png('galaxy3.png' ); smoothScatter(position$y,position$z);dev.off()
png('galaxy4.png' ); hist(position$x, breaks=50); dev.off()
png('galaxy5.png' ); hist(position$y, breaks=50); dev.off()
png('galaxy6.png' ); hist(position$z, breaks=50); dev.off()
png('galaxy7.png' ); hist(position$dist, breaks=50); dev.off()
png('galaxy8.png' ); hist(position$phi*180/pi, breaks=50); dev.off()
png('galaxy9.png' ); hist(position$l*180/pi, breaks=50); dev.off()
png('galaxy10.png' ); smoothScatter(position$phi*180/pi,position$g12e);dev.off()
png('galaxy11.png' ); smoothScatter(position$l*180/pi,position$g12e);dev.off()
png('galaxy12.png' ); smoothScatter(position$dist,position$g12e);dev.off()
png('galaxy13.png' ); smoothScatter(position$z,position$g12e);dev.off()
png('galaxy14.png' ); smoothScatter(position$phi*180/pi,position$g12g);dev.off()
png('galaxy15.png' ); smoothScatter(position$l*180/pi,position$g12g);dev.off()
png('galaxy16.png' ); smoothScatter(position$dist,position$g12g);dev.off()
png('galaxy17.png' ); smoothScatter(position$z,position$g12g);dev.off()
png('galaxy18.png' ); smoothScatter(position$phi*180/pi,position$ext,nbin=200);dev.off()
png('galaxy19.png' ); smoothScatter(position$l*180/pi,position$ext);dev.off()
png('galaxy20.png' ); smoothScatter(position$dist,position$ext);dev.off()
png('galaxy21.png' ); smoothScatter(position$z,position$ext);dev.off()
png('galaxy22.png' ); hist(log10(position$mass_1), breaks=50); dev.off()
png('galaxy23.png' ); hist(log10(position$mass_2), breaks=50); dev.off()
png('galaxy24.png' ); hist(position$V3, breaks=50); dev.off()
png('galaxy24+.png' ); hist(position$V3[position$mass_1>1], breaks=50); dev.off()
png('galaxy25.png' ); hist(position$age, breaks=50); dev.off()
png('galaxy26.png' ); hist(position$mag_abs_1, breaks=50); dev.off()
png('galaxy27.png' ); hist(position$mag_abs_2, breaks=50); dev.off()
png('galaxy28.png' ); hist(position$delta_mag, breaks=50); dev.off()
png('galaxy29.png' ); hist(position$app_mag_1, breaks=50); dev.off()
png('galaxy30.png' ); hist(position$app_mag_2, breaks=50); dev.off()
}
#write.table(position, file='position.txt', append=FALSE, col.names=T, row.names=FALSE, na="\"\"",sep='\t')
position$app_mag_1[position$app_mag_1>mag_limit_1]=NA ; 

if (scenario!='single') 
{position$app_mag_2[position$app_mag_2>mag_limit_2]=NA; 
position$delta_mag[position$delta_mag>delta_mag_limit]=NA}#Limiting apparent magnitude



###plots, black- generating, green- mass1 and mass2, blue- primary mass, red- secondary mass
###black - generating, green - all, blue-mass>1
if (num==1 & generation_plot==T & scenario!='single' & unresolved!=T) #mass and q plots  
{positions=subset(position, position$mass_1>1); positions5=subset(position, position$mass_1>5); mass_comb=c(position$mass_1, position$mass_2) 

postscript(file=paste('mass_', scenario, imf_1, power_2,'_',beta1,'.ps')) #mass generation plot
title=paste(scenario, imf_1, power_2); xlab='Mass'; ylab='Counts' #title for plot 
p=hist(mass, plot=FALSE, breaks=10**seq(-1,2,.1)) ; plot(p$mids, p$density, pch='.', col='black', cex=10, log='xy',ylim=c(10**-6,10), main=title, xlab=xlab, ylab=ylab); par(new=TRUE);
p=hist(mass_comb, plot=FALSE, breaks=10**seq(-1,2,.1)) ; plot(p$mids, p$density, pch=23, col='green', cex=2, log='xy', ylim=c(10**-6,10), main=title, xlab=xlab, ylab=ylab) ; par(new=TRUE);
p=hist(position$mass_1, plot=FALSE, breaks=10**seq(-1,2,.1)) ; plot(p$mids, p$density, pch='*', col='blue', cex=2, log='xy', ylim=c(10**-6,10), main=title, xlab=xlab, ylab=ylab) ; par(new=TRUE);
p=hist(position$mass_2, plot=FALSE, breaks=10**seq(-1,2,.1)) ; plot(p$mids, p$density, pch='*', col='red', cex=2, log='xy', ylim=c(10**-6,10), main=title, xlab=xlab, ylab=ylab); dev.off()

postscript(file=paste('q_',scenario, imf_1, power_2,'_',beta1,'.ps')); #q generation plot
title=paste(scenario, imf_1, power_2); xlab='q = mass_2/mass_1'; ylab='Counts' #title for plot 
if (scenario=='PCP' | scenario=='SCP') {p=hist(q, plot=FALSE, breaks=seq(0,1,.05)) ; plot(p$mids, p$density, pch='.', col='black', cex=10, log='y', main=title, xlab=xlab, ylab=ylab,ylim=c(10**-3,10)) ; par(new=TRUE)}
p=hist(position$V3, plot=FALSE, breaks=seq(0,1,.05)) ; plot(p$mids, p$density, pch='*', col='green', cex=2, log='y', main=title, xlab=xlab, ylab=ylab,ylim=c(10**-3,10)) ; par(new=TRUE)
p=hist(positions$V3, plot=FALSE, breaks=seq(0,1,.05)) ; plot(p$mids, p$density, pch='*', col='blue', cex=2, log='y', main=title, xlab=xlab, ylab=ylab,ylim=c(10**-3,10)) ; par(new=TRUE);
#p=hist(positions5$V3, plot=FALSE, breaks=seq(0,1,.05)) ; plot(p$mids, p$density, pch='*', col='violet', cex=2, log='y', main=title, xlab=xlab, ylab=ylab,ylim=c(10**-3,10)) 
dev.off()} ##plots



if (check==T)
{

neigh=subset(position, position$dist<local) ; thick=subset(neigh, neigh$th_d==2) ### local neighbourhood sample
#fwrite(neigh, file='neigh.txt', sep='\t', na='***', append=T) #catalogue
neigh_length[num]=length(neigh$dist); thick_length[num]= length(thick$th_d) ###length of local sample
if (num_limit< 5 & image == 1) {cat('local and full star number',length(neigh$dist), length(position$dist),'\n')}; nb_[num]=length(position$dist)}
if (scenario!='single' & check==T)
{### Local density calculation for comparison with Bovy
j=1; for (mass_dens_cut in mass_dens_list) ###Local mass density calculation, i is the variable
{mass_sample_1=subset(neigh, neigh[[1]]>mass_dens_cut) ; mass_sample_2=subset(mass_sample_1, mass_sample_1[[2]]>mass_dens_cut); #select systems with primary mass larger than given
mass_dens_1[j]=sum(mass_sample_1[[1]]) ; mass_dens_2[j]=sum(mass_sample_1[[1]])+sum(mass_sample_2[[2]]); j=j+1} #Proceed to next mass cut value
mass_dens_1_total[[num]]=mass_dens_1 ; mass_dens_2_total[[num]]=mass_dens_2} ### local neighbourhood sample finish

if (scenario=='single' & check==T)
{j=1; for (mass_dens_cut in mass_dens_list) ###Local mass density calculation, i is the variable
{neigh=subset(neigh, neigh$ev1<1.5)
mass_sample_1=subset(neigh, neigh[[1]]>mass_dens_cut) #select systems with primary mass larger than given
mass_dens_1[j]=sum(mass_sample_1[[1]]); j=j+1} #Proceed to next mass cut value
mass_dens_1_total[[num]]=mass_dens_1
} ### local neighbourhood sample finish

position=position[complete.cases(position),] #Exclude NaNs
fwrite(position,file=paste('final.txt'), sep='\t', na='***', append=T, col.names=F) ###Final catalogue, parallel


num=num+1; #creating final list
if (num>num_limit) {break} } ## end of outer cycle

if (check==T)
{
if (scenario!='single')
{mass_dens_1_total=data.frame(matrix(unlist(mass_dens_1_total), nrow=num_limit, byrow=T)); mass_dens_2_total=data.frame(matrix(unlist(mass_dens_2_total), nrow=num_limit, byrow=T)); 
mass_dens_1_total=10**4*mass_dens_1_total/(4/3*pi*local**3) ; mass_dens_2_total=10**4*mass_dens_2_total/(4/3*pi*local**3)}

if (scenario=='single')
{mass_dens_1_total=data.frame(matrix(unlist(mass_dens_1_total), nrow=num_limit, byrow=T))
mass_dens_1_total=10**4*mass_dens_1_total/(4/3*pi*local**3)
}}

if (check==T)
{if (scenario!='single')
{stellar_density=data.frame(mass_dens_1_total[[1]], mass_dens_2_total[[1]], mass_dens_1_total[[2]], mass_dens_2_total[[2]],mass_dens_1_total[[3]],mass_dens_2_total[[3]],mass_dens_1_total[[4]],mass_dens_2_total[[4]])}
if (scenario=='single')
{stellar_density=data.frame(mass_dens_1_total[[1]], mass_dens_1_total[[2]], mass_dens_1_total[[3]], mass_dens_1_total[[4]], nb_, neigh_length, thick_length)}

stellar_density=round(stellar_density,5)
fwrite(stellar_density, col.names=FALSE, row.names=FALSE, file='exit_table.txt', sep='\t', append=T, quote=F)} #Writing output to the file
} # Function end


#arg=expand.grid(image, scenario, m_split_1, m_split_2, imf_1); arg #arguments for function
#mapply(Testing_function, arg[[1]],arg[[2]], arg[[3]], arg[[4]], arg[[5]], arg[[6]], arg[[7]], arg[[8]])
#mapply(Testing_function, image)
mcmapply(Testing_function, mc.cores=cores, image) #parallel run

file.rename('final.txt','final+.txt')
file.rename('exit_table.txt','exit_table+.txt')

Sys.time()-clock

