

#################################################################################
# Model information
# Model name: Bayesian vine copula (BVC) model
# Model developer: Zhiyong Liu  
# Contact:   271683042@qq.com/ liuzhiy25@mail.sysu.edu.cn
# Version:   v1.0
# Language: R 
# Date: 2021-03-17
#article: Environmental Modelling & Software: A Hybrid Bayesian Vine Model for Water Level Prediction , 2021
#################################################################################


library(xts)
library(hydroTSM)
library(hydroGOF)
library(VineCopula)
library(fitdistrplus)
library(copula)
library(ensembleBMA)
library(CDVineCopulaConditional)

#develop a function  for fitting the best marginal distributions
best_fit_dist<-function (x,x1=10,x2=0.5) { 
  set.seed(1)
  library(fitdistrplus)
  # we know these data are normally distributed... 
  
  # let's compute some fits...
  require(MASS)
  dat=x
  normal<-tryCatch({  fitdist(dat,"norm",method = "mle")}, 
                   error = function(err) {return(NA)})
  
  gamma<-tryCatch({ fitdist(dat,"gamma",method = "mle")}, 
                  error = function(err) {return(NA)})
  
  lognormal<-tryCatch({  fitdist(dat,"lnorm",method = "mle")}, 
                      error = function(err) {return(NA)})
  weibull<-tryCatch({ fitdist(dat,"weibull",method = "mle")}, 
                    error = function(err) {return(NA)})
  
  
  
  fits=list( normal, gamma, lognormal,weibull )
  sim={}
  p_value={}
  for (i in 1:4 ) { 
    if   (!is.na (fits[[ i]][1])  ) {  
      #sim=cbind(sim, gofstat(fits[[ i]]  )$chisq )
      #sim=cbind(sim, gofstat(fits[[ i]]  )$chisq )
      sim =tryCatch({cbind(sim, gofstat(fits[[ i]]  )$chisq) },
                    error = function(err) {return(cbind(sim,-i))} )  
      
      
      #p_value =cbind(p_value, gofstat(fits[[ i]]  )$chisqpvalue )
      # p_value =tryCatch({cbind(p_value, gofstat(fits[[ i]]  )$chisqpvalue) },
      #                   error = function(err) {return(cbind(p_value,i))} )   
      
      
    } else {
      sim=cbind(sim, NA )
      #p_value=cbind(p_value,NA)
      
    }
  }
  
  chi=rbind(sim)
  
  
  colnames(chi)<- c('normal',"gamma", "lognormal",'Weibull' )
  # get the logliks for each model...
  #max_log=sapply(fits, function(i) i$loglik)
  #max_index<- which.max( max_log )
  max_index<- which.min( sim)
  
  #max_index=5
  paras=fits[[max_index]]
  if (max_index==1) {
    
    v_cdf=pnorm(x,mean=paras$estimate[1],sd=paras$estimate[2] )
    v1=pnorm(x1,mean=paras$estimate[1],sd=paras$estimate[2] )
    v2=qnorm(x2,mean=paras$estimate[1],sd=paras$estimate[2] )
    
  } else if (max_index==4) {
    v_cdf=pweibull(x,shape=paras$estimate[1],scale=paras$estimate[2] )
    v1=pweibull(x1,shape=paras$estimate[1],scale=paras$estimate[2] )
    v2=qweibull(x2,shape=paras$estimate[1],scale=paras$estimate[2] )
    
    
    
  }      else if (max_index==2) {
    
    v_cdf=pgamma(x,shape=paras$estimate[1],rate=paras$estimate[2] )
    v1=pgamma(x1,shape=paras$estimate[1],rate=paras$estimate[2] )
    v2=qgamma(x2,shape=paras$estimate[1],rate=paras$estimate[2] )
  }      else if (max_index==3) {
    
    v_cdf= plnorm(x,meanlog=paras$estimate[1],sdlog=paras$estimate[2] )
    v1= plnorm(x1,meanlog=paras$estimate[1],sdlog=paras$estimate[2] )
    v2= qlnorm(x2,meanlog=paras$estimate[1],sdlog=paras$estimate[2] )
  }    
  
  results<-list (v_cdf=v_cdf,v1=v1,v2=v2 , chi=chi)
  return(results)
  
  
}


#station2=read.table("D:\\SYSU_science\\SYSU_珠江口专项\\珠江口\\data\\input2020\\heyuan_water_level_daily.dat",
#                   sep=" ",header=F)
station2=read.table("D:\\SYSU_science\\SYSU_珠江口专项\\珠江口\\data\\input2020\\longchuan_monthly.txt",
                    sep=" ",header=F) 

#for daily 
month_seris= seq(as.Date("1989/1/1"), as.Date("2011/12/31"), "month")
#month_seris=format(month_seris, "%Y%m%d")

station2<-xts:::xts(as.vector(station2), order.by=month_seris)
#station3=station2['1989/2005']
#station3=station2
station3=station2



# start
station2=read.table("D:\\SYSU_science\\SYSU_珠江口专项\\珠江口\\data\\input2020\\heyuan_water_level_daily.dat",
                    sep=" ",header=F) 

#for daily 
month_seris= seq(as.Date("1989/1/1"), as.Date("2011/12/31"), "days")
#month_seris=format(month_seris, "%Y%m%d")

station2<-xts:::xts(as.vector(station2), order.by=month_seris)
#station3=station2['1989/2005']
#station3=station2
station3=daily2monthly(station2, FUN=mean, na.rm=TRUE)

index1=4 #for the dimension 
for (jj1 in 0:2) { #for lag 1 to lag 3 month ahead 
  lag=jj1  
  length1=length(station3)
  b1=as.numeric(station3[(index1+lag):(length1)])
  for (i in  (index1-1):1) { # the first colomn is the predicted one 
    
    # b0=as.numeric(station3[(i+lag):(length1+i-4-lag)])
    b0=as.numeric(station3[i:(length1+i-4-lag)])
    b1=cbind(b1, b0)
  }
  head(b1)
  
  
  
  b_cdf={}
  for (jj in 1:index1)
  { 
    b_cdf=cbind( b_cdf, best_fit_dist(b1[,jj])$v_cdf)
  }
  
  #b_cdf[,]
  
  #data=cbind(b_cdf[,4], b_cdf[,1], b_cdf[,2], b_cdf[,3])
  #data=cbind(b_cdf[,4], b_cdf[,3], b_cdf[,1], b_cdf[,2])
  
  # daily time series 
  #month_seris= seq(as.Date("1989/1/4"), as.Date("2011/12/31"), "days")
  
  # monthly time series 
  #month_seris= seq(as.Date("1989/4/1"), as.Date("2011/12/31"), "months")
  month_seris= seq(as.Date(paste("1989/",4+jj1,"/1/",  sep="")), as.Date("2011/12/31"), "months")
  
  #month_seris=format(month_seris, "%Y%m%d")
  
  aa1<-xts:::xts((b_cdf[,]), order.by=month_seris)
  #station3=station2['1989/2011']
  data=aa1
  data1=aa1['1989/2006']
  colnames(data1) <- c("Y1","Y2","X3","X4")
  colnames(data) <- c("Y1","Y2","X3","X4")
  
  
  # for 6 vine model for 4-d case  
  ListVines <- CDVineCondListMatrices(data,Nx=3,"CVine")
  
  final.ouput={}
  
  for (i in  1:length(ListVines$CVine)) {
    Matrix=ListVines$CVine[[i]]
    Matrix
    ## Not run:
    # Fit copula families for the defined vine:
    RVM <- CDVineCondFit(data1,Nx=3,Matrix=Matrix)
    summary(RVM)
    RVM$Matrix
    RVM$family
    # check
    identical(RVM$Matrix,Matrix)
    #  plot(RVM)
    
    
    # Set the values of the conditioning variables as those used for the calibration.
    # Order them with respect to RVM$Matrix, considering that is a C-Vine
    
    d=dim(RVM$Matrix)[1]
    cond1 <- data[,RVM$Matrix[(d+1)-1,(d+1)-1]]
    cond2 <- data[,RVM$Matrix[(d+1)-2,(d+1)-2]]
    cond3 <- data[,RVM$Matrix[(d+1)-3,(d+1)-3]]
    condition <- cbind(cond1,cond2, cond3)
    # Simulate the variables
    Sim1={}
    for (m1 in 1:50) { 
      Sim0 <- CDVineCondSim(RVM,condition)
      head(Sim0)
      Sim1=cbind(Sim1,Sim0[,1])
    }
    
    Sim2=rowMeans(Sim1)
    # Plot the simulated variables over the observed
    Sim2 <- data.frame(Sim2)
    #data<-as.matrix(data)
    # overplot(Sim,as.data.frame(b_cdf))
    
    
    final=best_fit_dist(b1[,1], x2=Sim2[,1])$v2
    final.ouput=cbind(final.ouput,final)
    print(i)
  }
  
  
  newdata_v0<-xts:::xts((b1[,]), order.by=month_seris)
  
  final.ouput1<-xts:::xts((final.ouput[,]), order.by=month_seris)
  #newdata_v=newdata_v0['2006/2011'][, 1]
  #final.ouput2=final.ouput1['2006/2011'][, ]
  newdata_v=newdata_v0 [, 1]
  final.ouput2=final.ouput1 [, ]
  
  cor(final.ouput2[,5],newdata_v[,1])
  plot(final.ouput2[,3],newdata_v[,1])
  
  
  
  
  inital1=00
  trainmonths=7 #could be changed 
  month_seris1=row.names(data.frame(final.ouput2[,1]))
  month_seris1=format(as.Date(month_seris1), "%Y%m%d")
  TestData <- ensembleData( forecasts = data.frame(final.ouput2) ,
                            dates = month_seris1   ,
                            observations = as.numeric(newdata_v),
                            forecastHour = inital1,
                            initializationTime = "00")
  
  TestFit <- ensembleBMAgamma( TestData, trainingDays = trainmonths,control = controlBMAgamma(startupSpeed = 1))
  
  TestForc <- quantileForecast( TestFit, TestData,quantiles = c( 0.025,0.05,0.166,0.5,0.833,0.95,0.975))
  
  # cor(tempTestForc[,4],newdata_v[traindays:length(newdata_v),1])
  #plot(tempTestForc[,4], newdata_v[traindays:length(newdata_v),1])
  
  TestForc1=xts:::xts(TestForc[,], order.by=as.Date(month_seris1))
  # write all output 
  # all_forecast=cbind(obs=newdata_v['2007/2011'], 
  #                    final.ouput2['2007/2011'],
  #                    tempTestForc1['2007/2011'],
  #                    anfis=anfis_predict 
  #                    )
  all_forecast=cbind(obs=newdata_v , 
                     final.ouput2 ,
                     TestForc1  )
  write.csv(all_forecast , file = file.path(paste("D:\\SYSU_science\\SYSU_珠江口专项\\珠江口\\data\\output2020\\server\\0423\\monthly\\monthly_training\\",
                                                  "Nouse_monthly_longchuan_Forc_ANFIS_lag",jj1,".csv", sep="")),
            sep=",", col.names = T, row.names = F ,append=F)
  
  
  
}



#------------ 
#calculate the measure metrics 

orginal_path="D:\\SYSU_science\\SYSU_珠江口专项\\珠江口\\data\\output2020\\server\\0423\\monthly\\" 

setwd(orginal_path)
files_outlier<-list.files(pattern='*csv')

library(hydroTSM)
library(hydroGOF)
rmse2=  {}
cor2=  {}
nse2= {}
for (i1 in 1:3) {  
  data1=read.csv( files_outlier[i1])[,]
  length(data1[,3])
  month_seris= seq(as.Date("2007/1/1"), as.Date("2011/12/31"), "months")
  # month_seris=format(month_seris, "%Y%m%d")
  
  data1<-xts:::xts(as.matrix(cbind(data1 )), order.by=as.Date(month_seris))
  
  
  copula_mean=xts:::xts(rowMeans(data1[, 2:7]), order.by=as.Date(month_seris))
  data1_1=cbind( (data1[,c(1,2:7,11,15) ]),copula_mean=copula_mean)
  head(data1_1)
  #data1_2=as.matrix(data1_1)
  data1_2=data1_1[,c(1,8,10, 9)]
  head(data1_2)
  measures_all1={}
  observations=data1_2[,1]
  measures_all={}
  for ( j in 2:4 ) { 
    rmse1= rmse(data1_2[,j], observations   ) 
    cor1= cor(data1_2[,j], observations  ) *cor(data1_1[,j], observations  ) 
    nse1= NSE(data1_2[,j], observations ) 
    
    measures=cbind(rmse1,cor=cor1,nse1  )
    measures_all=rbind(measures_all,measures  )
  }
  rmse2=rbind(rmse2, measures_all[,1])
  cor2=rbind(cor2, measures_all[,2])
  nse2=rbind(nse2, measures_all[,3])
} 
measures_all2=cbind(nse2,cor2, rmse2 )
write.table(measures_all2, "clipboard-128", sep='\t')
#--------------
#plot the time series  monthly  
save_dir="D:\\SYSU_science\\SYSU_珠江口专项\\珠江口\\data\\output2020\\server\\0423\\figures\\"
mypath_ave <- file.path( save_dir, 
                         paste('monthly_plot',  ".JPG", sep=""))
jpeg(file=mypath_ave,width = 4500 , height =  1500 , res = 200 )
par(mar=c(5,6,4,2)+0.1)
par(mfrow=c(1,2))

orginal_path="D:\\SYSU_science\\SYSU_珠江口专项\\珠江口\\data\\output2020\\server\\0423\\monthly\\" 

setwd(orginal_path)
files_outlier<-list.files()

library(hydroTSM)
library(hydroGOF)

data1=read.csv( files_outlier[1])[,]
length(data1[,3])
month_seris= seq(as.Date("2007/1/1"), as.Date("2011/12/31"), "months")
# month_seris=format(month_seris, "%Y%m%d")

data1<-xts:::xts(as.matrix(cbind(data1 )), order.by=as.Date(month_seris))


copula_mean=xts:::xts(rowMeans(data1[, 2:7]), order.by=as.Date(month_seris))
data1_1=cbind( (data1[,c(1,2:7,11,15) ]),copula_mean=copula_mean)
head(data1_1)
#data1_2=as.matrix(data1_1)
data1_2=data1_1[,c(1,8,10, 9)]
head(data1_2)


period='2007/2011'
data2=data1[period]
head(data2)
lower=as.numeric(data2$X0.05 )
upper=as.numeric(data2$X0.95)
seris2=data2$X0.5
x1=seq(as.Date("2007/1/1"), as.Date("2011/12/31"), "months")

copula_mean=rowMeans(data2[, 2:7])


plot.zoo(seris2, main = "(a) ", xaxt="n",  ylim=c(min(lower), max(upper)),  cex=2,cex.lab=2,  cex.axis=1.5,cex.main=2,
         xlab='Time',ylab='Water level (m)',   font.main=2, font.lab=2, font.sub=2) 
polygon(c(x1, rev(x1)), c(lower, rev(upper)),  col='gray', border=NA)
axis.Date(1,at=pretty(index(seris2),12)  ,
          labels=format(pretty(index(seris2),12),format="  %Y-%m"),
          las=1.5, cex.axis=1.5)
lines(x1,seris2,col=' red' , lwd=2)
lines(x1,copula_mean,col='darkgreen' , lwd=2, lty=1)
lines(x1,data2$anfis,col='blueviolet', lwd=2, lty=1)
points(x1,data2$b1,col='blue', pch=20)
# legend("topright", legend=c("90% uncertainty intervals"), fill=c("gray"), cex=1.5, 
#        density=c(NA), bty="n")

legend("topright", legend=c("90% uncertainty intervals"), fill=c("gray"), cex=1.5,  
       density=c(NA), bty="n")

legend('topleft', legend =c(' ','BVC','Arithmetic mean','ANFIS') , col = c('white','red','darkgreen','blueviolet'),  cex=1.5,lwd=3,
       lty = c(NA,1,1,1), bty="n")   #, trace = TRUE)

legend('topleft', legend =c('     Observations' ) , col = c('blue' ), cex=1.5,  
       pch=20, bty="n")   #, trace = TRUE)

#  longchuan

data1=read.csv( files_outlier[4])[,]
length(data1[,3])
month_seris= seq(as.Date("2007/1/1"), as.Date("2011/12/31"), "months")
# month_seris=format(month_seris, "%Y%m%d")

data1<-xts:::xts(as.matrix(cbind(data1 )), order.by=as.Date(month_seris))


copula_mean=xts:::xts(rowMeans(data1[, 2:7]), order.by=as.Date(month_seris))
data1_1=cbind( (data1[,c(1,2:7,11,15) ]),copula_mean=copula_mean)
head(data1_1)
#data1_2=as.matrix(data1_1)
data1_2=data1_1[,c(1,8,10, 9)]


period='2007/2011'
data2=data1[period]
head(data2)
lower=as.numeric(data2$X0.05 )
upper=as.numeric(data2$X0.95)
seris2=data2$X0.5
x1=seq(as.Date("2007/1/1"), as.Date("2011/12/31"), "months")

copula_mean=rowMeans(data2[, 2:7])


plot.zoo(seris2, main = "(b) ", xaxt="n",  ylim=c(min(lower), max(upper)),  cex=2,cex.lab=2,  cex.axis=1.5,cex.main=2,
         xlab='Time',ylab='Water level (m)',   font.main=2, font.lab=2, font.sub=2) 
polygon(c(x1, rev(x1)), c(lower, rev(upper)),  col='gray', border=NA)
axis.Date(1,at=pretty(index(seris2),12)  ,
          labels=format(pretty(index(seris2),12),format="  %Y-%m"),
          las=1.5, cex.axis=1.5)
lines(x1,seris2,col=' red' , lwd=2)
lines(x1,copula_mean,col='darkgreen' , lwd=2, lty=1)
lines(x1,data2$anfis,col='blueviolet', lwd=2, lty=1)
points(x1,data2$b1,col='blue', pch=20)
# legend("topright", legend=c("90% uncertainty intervals"), fill=c("gray"), cex=1.5, 
#        density=c(NA), bty="n")

legend("topright", legend=c("90% uncertainty intervals"), fill=c("gray"), cex=1.5,  
       density=c(NA), bty="n")

legend('topleft', legend =c(' ','BVC','Arithmetic mean','ANFIS') , col = c('white','red','darkgreen','blueviolet'),  cex=1.5, lwd=3,
       lty = c(NA,1,1,1), bty="n")   #, trace = TRUE)

legend('topleft', legend =c('     Observations' ) , col = c('blue' ),  cex=1.5,  
       pch=20, bty="n")   #, trace = TRUE)

dev.off()

