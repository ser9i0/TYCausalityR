#!/usr/bin/env Rscript

library(tseries)
library(urca)
library(aod)
library(fUnitRoots)
library(zoo)
library(vars)
source("stationarity.R",chdir=T)


# TEST DE CAUSALIDAD

# Test de causalidad de Granger cuando se trata de series temporales no-estacionarias 
# (ya sean cointegradas o no) basado en el procedimiento descrito por Toda-Yamamoto (T-Y).

# El matiz mas importante es que si se usa el test de Wald para comprobar restricciones lineales
# en un modelo VAR (vector autoregression) y una de las series es no-estacionaria, entonces la
# estadística del test de Wald no sigue la usual distribución chi-cuadrado asintotica bajo el valor nulo

# La causalidad se analizará entre dos series temporales.
# Por un lado la serie de estudio, y por otro la serie objetivo.

args<-commandArgs(trailingOnly=TRUE)

if(length(args)==0){
    stop("Invalid arguments: datafile_path init_col (2) end_col (-1) target_col (1).")
}

dataset<-read.csv( args[1] )
init_col<-as.numeric(args[2])
end_col<-as.numeric(args[3])
target_col<-as.numeric(args[4])

order_integration <- vector(mode="list", length=length(names(dataset)))
names(order_integration) <- names(dataset)
target_name<-names(dataset)[target_col]

if( end_col == -1 )
{
  end_col<-ncol(dataset)
}

#message("init_col=",init_col,", end_col=",end_col,", target_col=",target_col)

if(init_col>end_col)
{
  stop("Wrong INIT (",init_col,") and/or END (",end_col,") columns value.")
}

#if((target_col >= init_col) && (target_col <= end_col))
#{
#  stop("Wrong TARGET (",target_col,"), it must be out of the range [INIT,END] ([",init_col,",",end_col,"]).")
#}

# Plot all rows
#plot(seq(1,nrow(dataset),1), dataset$"target1", type="l")

# Plot first 10 values of "feature1"
#plot(seq(1,10,1), cof[1:10,]$"feature1", type="l")

message("# STATIONARITY TESTS ############################")

# First, check stationarity of target time series
message("..TARGET: ",target_name," stationarity:")
stationarity<-stationary_test(dataset[,target_col])
target_order_integration<-stationarity$max_order_integration
order_integration[[{names(dataset)[target_col]}]]<-target_order_integration
message("-----------------------------------------------")

for( i in init_col:end_col )
{
  if( i == target_col )
  {
    next
  }

  col_name<-names(dataset)[i]
  message("..",col_name," stationarity:")

  stationarity<-stationary_test(dataset[,i])

  # Johansen Cointegration Test
  cjtest<-ca.jo(dataset[,c(i,target_col)])
  r1test<-cjtest@teststat[[1]] # r <= 1
  r1_perc10<-cjtest@cval[[1]]
  r0test<-cjtest@teststat[[2]] # r = 0
  r0_perc10<-cjtest@cval[[2]]
  #message("r<=1: ",r1test," r<=1 perc10: ",r1_perc10," - r=0: ",r0test," r=0 perc10: ",r0_perc10)
  if( (r1test < r1_perc10) && (r0test < r0_perc10) )
  {
    message("....Cointegration: time-series ",target_name," and ",col_name," are not cointegrated.")
  }
  else if( (r1test < r1_perc10) && (r0test >= r0_perc10) )
  {
    message("....Cointegration: time-series ",target_name," and ",col_name," are not cointegrated although r=0 hypothesis has been rejected.")
  }
  else
  {
    message("....Cointegration: time-series ",target_name," and ",col_name," are cointegrated with significance level >=10%.")
  }
  

  if( stationarity$max_order_integration >= target_order_integration )
  {
    order_integration[[{names(dataset)[i]}]]<-stationarity$max_order_integration
  }
  else
  {
    order_integration[[{names(dataset)[i]}]]<-target_order_integration
    message("INFO: Order of integration of TARGET time-series is higher (",target_order_integration,"), this one will be used.")
  }

  message("-----------------------------------------------")
}

# After testing stationarity, then set up a VAR model
# Regardless of the order of integration, the data will not be differenced

pdf(paste(target_name,"causality.pdf",sep="_"))
plot( seq(1,nrow(dataset),1), 
      dataset[,target_col], 
      #ylim=c(0.5,2), 
      type="l",
      lty = 1,
      lwd = 2,
      col="black",
      main=paste(target_name,"causality analysis"),
      col.lab="red",
      cex.lab=0.75,
      )
par(new = TRUE)

colors<-colors()
graph_colors<-c("black")
graph_labels<-c(target_name)

message("# SET UP VAR-model ##############################")
for( i in init_col:end_col )
{
  if( i == target_col )
  {
    next
  }
  # VARselect returns infomation criteria and final prediction error for 
  # sequential increasing the lag order up to a VAR(p)-proccess.
  var_select<-VARselect(dataset[,c(i,target_col)],lag=20,type="both")

if(FALSE) {
  message("-- v1 ---------------------------------------------")
  v1<-VAR(dataset[,c(i,target_col)],p=1,type="both")
  print(serial.test(v1))
  print(1/roots(v1)[[1]])
  print(1/roots(v1)[[2]])
  message("-- v2 ---------------------------------------------")
  v2<-VAR(dataset[,c(i,target_col)],p=2,type="both")
  print(serial.test(v2))
  print(1/roots(v2)[[1]])
  print(1/roots(v2)[[2]])
  message("-- v3 ---------------------------------------------")
  v3<-VAR(dataset[,c(i,target_col)],p=3,type="both")
  print(serial.test(v3))
  print(1/roots(v3)[[1]])
  print(1/roots(v3)[[2]])
  message("-- v4 ---------------------------------------------")
  v4<-VAR(dataset[,c(i,target_col)],p=4,type="both")
  print(serial.test(v4))
  print(1/roots(v4)[[1]])
  print(1/roots(v4)[[2]])
  message("-- v5 ---------------------------------------------")
  v5<-VAR(dataset[,c(i,target_col)],p=5,type="both")
  print(serial.test(v5))
  print(1/roots(v5)[[1]])
  print(1/roots(v5)[[2]])
  message("-- v6 ---------------------------------------------")
  v6<-VAR(dataset[,c(i,target_col)],p=6,type="both")
  print(serial.test(v6))
  print(1/roots(v6)[[1]])
  print(1/roots(v6)[[2]])
  message("-- v7 ---------------------------------------------")
  v7<-VAR(dataset[,c(i,target_col)],p=7,type="both")
  print(serial.test(v7))
  print(1/roots(v7)[[1]])
  print(1/roots(v7)[[2]])
  message("-- v8 ---------------------------------------------")
  v8<-VAR(dataset[,c(i,target_col)],p=8,type="both")
  print(serial.test(v8))
  print(1/roots(v8)[[1]])
  print(1/roots(v8)[[2]])
  message("-- v9 ---------------------------------------------")
  v9<-VAR(dataset[,c(i,target_col)],p=9,type="both")
  print(serial.test(v9))
  print(1/roots(v9)[[1]])
  print(1/roots(v9)[[2]])
  message("-- v10 ---------------------------------------------")
  v10<-VAR(dataset[,c(i,target_col)],p=10,type="both")
  print(serial.test(v10))
  print(1/roots(v10)[[1]])
  print(1/roots(v10)[[2]])
  message("-- v11 ---------------------------------------------")
  v11<-VAR(dataset[,c(i,target_col)],p=11,type="both")
  print(serial.test(v11))
  print(1/roots(v11)[[1]])
  print(1/roots(v11)[[2]])
  message("-- v12 ---------------------------------------------")
  v12<-VAR(dataset[,c(i,target_col)],p=12,type="both")
  print(serial.test(v12))
  print(1/roots(v12)[[1]])
  print(1/roots(v12)[[2]])
  message("-- v13 ---------------------------------------------")
  v13<-VAR(dataset[,c(i,target_col)],p=13,type="both")
  print(serial.test(v13))
  print(1/roots(v13)[[1]])
  print(1/roots(v13)[[2]])
  message("-- v14 ---------------------------------------------")
  v14<-VAR(dataset[,c(i,target_col)],p=14,type="both")
  print(serial.test(v14))
  print(1/roots(v14)[[1]])
  print(1/roots(v14)[[2]])
  message("-- v15 ---------------------------------------------")
  v15<-VAR(dataset[,c(i,target_col)],p=15,type="both")
  print(serial.test(v15))
  print(1/roots(v15)[[1]])
  print(1/roots(v15)[[2]])
  message("-- v16 ---------------------------------------------")
  v16<-VAR(dataset[,c(i,target_col)],p=16,type="both")
  print(serial.test(v16))
  print(1/roots(v16)[[1]])
  print(1/roots(v16)[[2]])
  message("-- v17 ---------------------------------------------")
  v17<-VAR(dataset[,c(i,target_col)],p=17,type="both")
  print(serial.test(v17))
  print(1/roots(v17)[[1]])
  print(1/roots(v17)[[2]])
  message("-- v18 ---------------------------------------------")
  v18<-VAR(dataset[,c(i,target_col)],p=18,type="both")
  print(serial.test(v18))
  print(1/roots(v18)[[1]])
  print(1/roots(v18)[[2]])
  message("-- v19 ---------------------------------------------")
  v19<-VAR(dataset[,c(i,target_col)],p=19,type="both")
  print(serial.test(v19))
  print(1/roots(v19)[[1]])
  print(1/roots(v19)[[2]])
  message("-- v20 ---------------------------------------------")
  v20<-VAR(dataset[,c(i,target_col)],p=20,type="both")
  print(serial.test(v20))
  print(1/roots(v20)[[1]])
  print(1/roots(v20)[[2]])
}

  # Among the infomation criteria, we find:
  # - AIC (Akaike Information criterion)
  # - SC (Schwarz information criterion)
  # - HQ (Hannan-Quin information criterion)
  # - FPE (Final Prediction Error criterion)
  # The function returns the lag length that minimizes the values of the criteria for the VAR model.
  # AIC and SC are mostly used, although they may give conflicting results

  # The lag length ought to be set such that the VAR residuals are free of autocorrelation, 
  # even if this implies longer lags than suggested by the information criteria
  # - serial.test(var_model) (the higher P value, the better)
  # - 1/roots(var_model)[[1]] ( must be > 1)
  # - 1/roots(var_model)[[2]] ( must be > 1)

  ord_int<-order_integration[[{names(dataset)[i]}]]
  AIC_selection<-var_select$selection[1]
  p_value<-AIC_selection + ord_int
  message("..",names(dataset)[i]," var selection: ",var_select$selection," p_value (k+dmax): ",p_value)

  var_model<-VAR(dataset[,c(i,target_col)],p=p_value,type="both")

  # Serial correlation: check that the VAR residuals are free of autocorrelation
  if( serial.test(var_model)[["serial"]][["p.value"]] <= 0.1 )
  {
    message("WARNING!: The VAR residuals are not free of autocorrelation. Please check the VAR model chosen with lag = AIC(",AIC_selection,")+ m(",ord_int,").")
  }
  # Stability analysis
  if( (1/roots(var_model)[[1]] <=1) || (1/roots(var_model)[[2]] <=1) )
  {
    warning("WARNING!: Stability analysis failed. Please check the VAR model chosen with lag = AIC(",AIC_selection,") + m(",ord_int,").")
  }

  wt1<-wald.test(b=coef(var_model$varresult[[1]]), Sigma=vcov(var_model$varresult[[1]]), Terms=c(seq(2, AIC_selection*2, 2)))
  wt2<-wald.test(b=coef(var_model$varresult[[2]]), Sigma=vcov(var_model$varresult[[2]]), Terms=c(seq(1, AIC_selection*2, 2)))

  wt1_p<-wt1[["result"]][["chi2"]][["P"]]
  wt2_p<-wt2[["result"]][["chi2"]][["P"]]

  if( wt1_p <= 0.1 )
  {
    message( ">>> ",names(dataset)[target_col]," Granger-causes ",names(dataset)[i]," at the 10% significance level (p = ",wt1_p,")." )
  }

  if( wt2_p <= 0.1 )
  {
    message( ">>> ",names(dataset)[i]," Granger-causes ",names(dataset)[target_col]," at the 10% significance level (p = ",wt2_p,")." )
    plot( seq(1,nrow(dataset),1), 
          dataset[,i],
          type="l",
          col = colors[i*3], 
          lty = 1, 
          lwd = 1,
          xlab = "",
          ylab = "")
    par(new = TRUE)
    graph_colors<-c(graph_colors,colors[i*3])
    graph_labels<-c(graph_labels,paste(names(dataset)[i],"(p ",round(wt2_p,3),")"))
  }

  message("-----------------------------------------------")
}

#plot(seq(1,nrow(dataset),1), dataset[,4], ylim=c(0.5,2), type="l", col="black", lwd=2)
#lines(seq(1,nrow(dataset),1), dataset[,7], col="blue", lty=2, lwd=1)
legend("topleft",graph_labels,col=graph_colors,lty=1,bty="n")


