library(maptools)
library(sf)
library(raster)
library(readxl)
library(writexl)
library(dplyr)
library(rugarch)
library(tseries)
library(fUnitRoots)
library(FinTS)
library(car)
library(Kendall)
library(RobustLinearReg)
library(mosum)
library(lomb)
# Required formulas throughout the script
median_by_row <- function(row) {
  return(median(row, na.rm = TRUE))
}
mean_by_row <- function(row) {
  return(mean(row, na.rm = TRUE))
}
min_by_row <- function(row) {
  return(min(row, na.rm = TRUE))
}
max_by_row <- function(row) {
  return(max(row, na.rm = TRUE))
}
minimum_values<-function(vec,num){
  sorted_vec <- sort(vec)
  min_values <- sorted_vec[1:num]
  return(min_values)
}
detect_local_maxima <- function(x) {
  maxima_indices <- which(diff(sign(diff(x))) == -2) + 1
  return(maxima_indices)
}
calculate_negative_semivariance <- function(returns) {
  negative_returns <- returns[returns < mean(returns)]
  if (length(negative_returns) > 0) {
    return(sum(negative_returns^2) / length(negative_returns))
  } else {
    return(0)
  }
}

################################ FORMULAS ########################

#Number of Seasons using the Lomb-Scargle periodogram
Seasonality<-function(Number_years,NDVI_TS){#Number_years=Number of complete years // NDVI_TS=Time series of the NDVI
  lsp_result<-lsp(NDVI_TS,plot=F)
  peak_index <- which.max(lsp_result$power)
  peak_frequency <- lsp_result$scanned[peak_index]
  peak_period1 <- 1 / peak_frequency
  peak_period1<-round(peak_period1) #Length of the season
  length_year<-length(NDVI_TS)/Number_years #Length of the year
  seasons_per_year<-peak_period1/length_years
  return(seasons_per_year)
}

TS<-readRDS("C:/PHD_REGIONES/ts_ejemplo.rds")
NDVI_TS<-TS[36:(length(TS)-2)]
NSeasons<-1
Number_years<-21
pct_amplitude<-20

#Extract SOS, EOS, DOS and peak for each year

pheno_parameters<-function(NSeasons,NDVI_TS,Number_years,pct_amplitude){#NSeasons= 1 or 2 // NDVI_TS=Time series of the NDVI // Number_years=Number of complete years  // pct_amplitude= percentage of the amplitude
  pct_amplitude<-pct_amplitude/100
  length_year<-length(NDVI_TS)/Number_years #How many period does a year have
  dates<-rep(1:length_year,times=Number_years)
  year<-rep(c(1:Number_years),each=length_year)
  df<-as.data.frame(cbind(year,dates,NDVI_TS))
  colnames(df)<-c("Year","Dates","TS")
  
  if(NSeasons==1){
    Season<-data.frame(matrix(ncol=5,nrow=0))
    for (j in 1:(Number_years-1)){
      NDVI_year<-subset(df,df$Year==j)
      max_index <- which.max(NDVI_year$TS)
      peak_number<-max_index
      row_max<-as.numeric(rownames(subset(df,df$Year==j & df$Dates==max_index)))
      period_before<-max(round(length_year/2,0),max_index) #This row and the following one are used to take at least the calendar year.
      period_after<-max(round(length_year/2,0),length_year-max_index)
      NDVI_cicle<-df[max(1,(row_max-period_before+1)):(row_max+period_after-1),]
      
      max_value<-max(NDVI_cicle$TS)
      max_index<-which.max(NDVI_cicle$TS)
      min_value_before<-min(NDVI_cicle$TS[1:max_index])
      min_index_before<-which(NDVI_cicle$TS==min_value_before)
      min_value_after<-min(NDVI_cicle$TS[max_index:length(NDVI_cicle$TS)])
      min_index_after<-which(NDVI_cicle$TS==min_value_after)
      
      NDVI_cicle<-NDVI_cicle[min_index_before:min_index_after,]
      max_index<-which.max(NDVI_cicle$TS)
      Threshold_before<-min_value_before+pct_amplitude*(max_value-min_value_before)
      Threshold_after<-min_value_after+pct_amplitude*(max_value-min_value_after)
      SOS_index <- which.min(abs(NDVI_cicle$TS[1:max_index] - Threshold_before))
      EOS_index <- which.min(abs(NDVI_cicle$TS[(max_index + 1):length(NDVI_cicle$TS)] - Threshold_after)) + max_index
      SOS_date<-NDVI_cicle[SOS_index,2]
      EOS_date<-NDVI_cicle[EOS_index,2]
      EOS_date<-ifelse(EOS_date>length_year,length_year+EOS_date,EOS_date)
      DOS<-EOS_date-SOS_date
      Season1<-c(j,SOS_date,EOS_date,DOS,df[peak_number,2])
      Season<-rbind(Season,Season1)
      colnames(Season)<-c("Year","SOS","EOS","DOS","PEAK")
    }
  }else if (NSeasons==2){
    Season<-data.frame(matrix(ncol=9,nrow=0))
    for (j in 1:(Number_years-1)){
      NDVI_year<-subset(df,df$Year==j)
      N_Max_local<-length(detect_local_maxima(NDVI_year$TS))
      if (N_Max_local==2){
        TS<-as.numeric(NDVI_year$TS)
      }else{
        TS<-as.numeric(NDVI_year$TS)
        while (N_Max_local>2){
          TS1<-append(TS[1],TS[1])
          TS1<-append(TS1,TS)
          TS1<-append(TS1,TS[length(TS)])
          TS1<-append(TS1,TS[length(TS)])
          TS <- zoo::rollmean(TS1, k = 5)
          N_Max_local<-length(detect_local_maxima(TS))
        }
      }
      Max_local<-detect_local_maxima(TS)
      
      if (length(Max_local)!=2){
        return(print(paste0("In year ",j, " there are not 2 seasons")))
      }else{
        Pos_Max_season1<-Max_local[1]
        peak1_number<-Pos_Max_season1
        Pos_Max_season2<-Max_local[2]
        peak2_number<-Pos_Max_season2
      }
      Max_season1<-max(NDVI_year$TS[1:(Pos_Max_season1+2)]) #Margin of 2 for the moving average window
      Max_season2<-max(NDVI_year$TS[(Pos_Max_season2-2):(length(NDVI_year$TS)-2)])
      
      #Season1
      period_before<-max(round(length_year/2,0),Pos_Max_season1) #This row and the following one are used to take at least the calendar year.
      period_after<-max(round(length_year/4,0),Pos_Max_season2-Pos_Max_season1)
      row_max<-as.numeric(rownames(subset(df,df$Year==j & df$Dates==Pos_Max_season1))) 
      NDVI_cicle<-df[max(1,(row_max-period_before+1)):(row_max+period_after-1),]
      
      max_index<-which(NDVI_cicle$TS==Max_season1)
      min_value_before<-min(NDVI_cicle$TS[1:max_index])
      min_index_before<-which(NDVI_cicle$TS==min_value_before)
      min_value_after<-min(NDVI_cicle$TS[max_index:length(NDVI_cicle$TS)])
      min_index_after<-which(NDVI_cicle$TS==min_value_after)
      NDVI_cicle<-NDVI_cicle[min_index_before:min_index_after,]
      max_index<-which(NDVI_cicle$TS==Max_season1)
      Threshold_before<-min_value_before+pct_amplitude*(Max_season1-min_value_before)
      Threshold_after<-min_value_after+pct_amplitude*(Max_season1-min_value_after)
      SOS1_index <- which.min(abs(NDVI_cicle$TS[1:max_index] - Threshold_before))
      EOS1_index <- which.min(abs(NDVI_cicle$TS[(max_index + 1):length(NDVI_cicle$TS)] - Threshold_after)) + max_index
      SOS1_date<-NDVI_cicle[SOS1_index,2]
      EOS1_date<-NDVI_cicle[EOS1_index,2]
      DOS1<-EOS1_date-SOS1_date
      
      #Season2
      period_before<-max(round(length_year/4,0),Pos_Max_season2-Pos_Max_season1) 
      period_after<-max(round(length_year/2,0),length_year-Pos_Max_season2)
      row_max<-as.numeric(rownames(subset(df,df$Year==j & df$Dates==Pos_Max_season2))) 
      NDVI_cicle<-df[max(1,(row_max-period_before+1)):(row_max+period_after-1),]
      
      max_index<-which(NDVI_cicle$TS==Max_season2)
      min_value_before<-min(NDVI_cicle$TS[1:max_index])
      min_index_before<-which(NDVI_cicle$TS==min_value_before)
      min_value_after<-min(NDVI_cicle$TS[max_index:length(NDVI_cicle$TS)])
      min_index_after<-which(NDVI_cicle$TS==min_value_after)
      NDVI_cicle<-NDVI_cicle[min_index_before:min_index_after,]
      max_index<-which(NDVI_cicle$TS==Max_season2)
      Threshold_before<-min_value_before+pct_amplitude*(Max_season2-min_value_before)
      Threshold_after<-min_value_after+pct_amplitude*(Max_season2-min_value_after)
      SOS2_index <- which.min(abs(NDVI_cicle$TS[1:max_index] - Threshold_before))
      EOS2_index <- which.min(abs(NDVI_cicle$TS[(max_index + 1):length(NDVI_cicle$TS)] - Threshold_after)) + max_index
      SOS2_date<-NDVI_cicle[SOS2_index,2]
      EOS2_date<-NDVI_cicle[EOS2_index,2]
      EOS2_date<-ifelse(EOS2_date>length_year,length_year+EOS2_date,EOS2_date)
      DOS2<-EOS2_date-SOS2_date
      
      Season<-rbind(Season,c(j,SOS1_date,EOS1_date,DOS1,peak1_number,SOS2_date,EOS2_date,DOS2,peak2_number))
      colnames(Season)<-c("Year","SOS1","EOS1","DOS1","PEAK1","SOS2","EOS2","DOS2","PEAK2")
    }
  }else{
      return(print(paste0("NSeasons has to be 1 or 2")))
  }
  return(Season)
}

#Is there a trend in SOS, EOS, DOS and PEAK? Mean and medians
trend_pheno_parameters<-function(Season){ #The season argument has to be a dataframe such as the one obtained using the "pheno_parameters" formula
  if (ncol(Season)==5){ #There is only 1 season
    trend<-data.frame(matrix(nrow=0,ncol=16))
    
    theil_sen<-theil_sen_regression(SOS~Year,data=Season)
    coef_SOS<-theil_sen$coefficients[2]
    pvalue_SOS<- as.numeric(MannKendall(Season$SOS)$sl)  #Mann-Kendall test
    
    theil_sen<-theil_sen_regression(EOS~Year,data=Season)
    coef_EOS<-theil_sen$coefficients[2]
    pvalue_EOS<- as.numeric(MannKendall(Season$EOS)$sl)  #Mann-Kendall test
    
    theil_sen<-theil_sen_regression(DOS~Year,data=Season)
    coef_DOS<-theil_sen$coefficients[2]
    pvalue_DOS<- as.numeric(MannKendall(Season$DOS)$sl)  #Mann-Kendall test
    
    theil_sen<-theil_sen_regression(PEAK~Year,data=Season)
    coef_PEAK<-theil_sen$coefficients[2]
    pvalue_PEAK<- as.numeric(MannKendall(Season$PEAK)$sl)  #Mann-Kendall test
    
    SOS_median<-round(median(Season$SOS),0)
    EOS_median<-round(median(Season$EOS),0)
    DOS_median<-round(median(Season$DOS),0)
    PEAK_median<-round(median(Season$PEAK),0)
    
    SOS_mean<-round(mean(Season$SOS),0)
    EOS_mean<-round(mean(Season$EOS),0)
    DOS_mean<-round(mean(Season$DOS),0)
    PEAK_mean<-round(mean(Season$PEAK),0)
    
    trend<-rbind(trend,c(SOS_mean,SOS_median,EOS_mean,EOS_median,DOS_mean,DOS_median,PEAK_mean,PEAK_median,coef_SOS,pvalue_SOS,coef_EOS,pvalue_EOS,coef_DOS,pvalue_DOS,coef_PEAK,pvalue_PEAK))
    colnames(trend)<-c("SOS_mean","SOS_median","EOS_mean","EOS_median","DOS_mean","DOS_median","PEAK_mean","PEAK_median","coef_SOS","pvalue_SOS","coef_EOS","pvalue_EOS","coef_DOS","pvalue_DOS","coef_PEAK","pvalue_PEAK")
    return(trend)
    
  } else if (ncol(Season)==9){#There are 2 seasons
    trend<-data.frame(matrix(nrow=0,ncol=44))
    
    Season$SOSC<-Season$SOS1  #2 seasons together
    Season$EOSC<-Season$EOS2  #2 seasons together
    Season$DOSC<- Season$EOS2- Season$SOS1  #2 seasons together
    
    #Season 1
    theil_sen<-theil_sen_regression(SOS1~Year,data=Season)
    coef_SOS1<-theil_sen$coefficients[2]
    pvalue_SOS1<- as.numeric(MannKendall(Season$SOS1)$sl)  #Mann-Kendall test
    
    theil_sen<-theil_sen_regression(EOS1~Year,data=Season)
    coef_EOS1<-theil_sen$coefficients[2]
    pvalue_EOS1<- as.numeric(MannKendall(Season$EOS1)$sl)  #Mann-Kendall test
    
    theil_sen<-theil_sen_regression(DOS1~Year,data=Season)
    coef_DOS1<-theil_sen$coefficients[2]
    pvalue_DOS1<- as.numeric(MannKendall(Season$DOS1)$sl)  #Mann-Kendall test
    
    theil_sen<-theil_sen_regression(PEAK1~Year,data=Season)
    coef_PEAK1<-theil_sen$coefficients[2]
    pvalue_PEAK1<- as.numeric(MannKendall(Season$PEAK1)$sl) 
    
    SOS1_median<-round(median(Season$SOS1),0)
    EOS1_median<-round(median(Season$EOS1),0)
    DOS1_median<-round(median(Season$DOS1),0)
    PEAK1_median<-round(median(Season$PEAK1),0)
    
    SOS1_mean<-round(mean(Season$SOS1),0)
    EOS1_mean<-round(mean(Season$EOS1),0)
    DOS1_mean<-round(mean(Season$DOS1),0)
    PEAK1_mean<-round(mean(Season$PEAK1),0)
    
    #Season 2
    theil_sen<-theil_sen_regression(SOS2~Year,data=Season)
    coef_SOS2<-theil_sen$coefficients[2]
    pvalue_SOS2<- as.numeric(MannKendall(Season$SOS2)$sl)  #Mann-Kendall test
    
    theil_sen<-theil_sen_regression(EOS2~Year,data=Season)
    coef_EOS2<-theil_sen$coefficients[2]
    pvalue_EOS2<- as.numeric(MannKendall(Season$EOS2)$sl)  #Mann-Kendall test
    
    theil_sen<-theil_sen_regression(DOS2~Year,data=Season)
    coef_DOS2<-theil_sen$coefficients[2]
    pvalue_DOS2<- as.numeric(MannKendall(Season$DOS2)$sl)  #Mann-Kendall test
    
    theil_sen<-theil_sen_regression(PEAK2~Year,data=Season)
    coef_PEAK2<-theil_sen$coefficients[2]
    pvalue_PEAK2<- as.numeric(MannKendall(Season$PEAK2)$sl) 
    
    SOS2_median<-round(median(Season$SOS2),0)
    EOS2_median<-round(median(Season$EOS2),0)
    DOS2_median<-round(median(Season$DOS2),0)
    PEAK2_median<-round(median(Season$PEAK2),0)
    
    SOS2_mean<-round(mean(Season$SOS2),0)
    EOS2_mean<-round(mean(Season$EOS2),0)
    DOS2_mean<-round(mean(Season$DOS2),0)
    PEAK2_mean<-round(mean(Season$PEAK2),0)
    
    #Complete
    theil_sen<-theil_sen_regression(SOSC~Year,data=Season)
    coef_SOSC<-theil_sen$coefficients[2]
    pvalue_SOSC<- as.numeric(MannKendall(Season$SOSC)$sl)  #Mann-Kendall test
    
    theil_sen<-theil_sen_regression(EOSC~Year,data=Season)
    coef_EOSC<-theil_sen$coefficients[2]
    pvalue_EOSC<- as.numeric(MannKendall(Season$EOSC)$sl)  #Mann-Kendall test
    
    theil_sen<-theil_sen_regression(DOSC~Year,data=Season)
    coef_DOSC<-theil_sen$coefficients[2]
    pvalue_DOSC<- as.numeric(MannKendall(Season$DOSC)$sl)  #Mann-Kendall test
    
    SOSC_median<-round(median(Season$SOSC),0)
    EOSC_median<-round(median(Season$EOSC),0)
    DOSC_median<-round(median(Season$DOSC),0)
    
    SOSC_mean<-round(mean(Season$SOSC),0)
    EOSC_mean<-round(mean(Season$EOSC),0)
    DOSC_mean<-round(mean(Season$DOSC),0)
    
    trend<-rbind(trend,c(SOS1_mean,SOS1_median,EOS1_mean,EOS1_median,DOS1_mean,DOS1_median,PEAK1_mean,PEAK1_median,coef_SOS1,pvalue_SOS1,coef_EOS1,pvalue_EOS1,coef_DOS1,pvalue_DOS1,coef_PEAK1,pvalue_PEAK1,SOS2_mean,SOS2_median,EOS2_mean,EOS2_median,DOS2_mean,DOS2_median,PEAK2_mean,PEAK2_median,coef_SOS2,pvalue_SOS2,coef_EOS2,pvalue_EOS2,coef_DOS2,pvalue_DOS2,coef_PEAK2,pvalue_PEAK2,SOSC_mean,SOSC_median,EOSC_mean,EOSC_median,DOSC_mean,DOSC_median,coef_SOSC,pvalue_SOSC,coef_EOSC,pvalue_EOSC,coef_DOSC,pvalue_DOSC))
    colnames(trend)<-c("SOS1_mean","SOS1_median","EOS1_mean","EOS1_median","DOS1_mean","DOS1_median","PEAK1_mean","PEAK1_median","coef_SOS1","pvalue_SOS1","coef_EOS1","pvalue_EOS1","coef_DOS1","pvalue_DOS1","coef_PEAK1","pvalue_PEAK1","SOS2_mean","SOS2_median","EOS2_mean","EOS2_median","DOS2_mean","DOS2_median","PEAK2_mean","PEAK2_median","coef_SOS2","pvalue_SOS2","coef_EOS2","pvalue_EOS2","coef_DOS2","pvalue_DOS2","coef_PEAK2","pvalue_PEAK2","SOSC_mean","SOSC_median","EOSC_mean","EOSC_median","DOSC_mean","DOSC_median","coef_SOSC","pvalue_SOSC","coef_EOSC","pvalue_EOSC","coef_DOSC","pvalue_DOSC")
    return(trend)
    
  } else{
    return(print("There has to be 1 or 2 seasons per year"))
  }
    
  
}

####### ANALYZE WET SEASON ###########

parameters_wet_season<-function(SOS,EOS,NDVI_TS,Number_years){
  length_year<-length(NDVI_TS)/Number_years #How many periods does a year have
  dates<-rep(1:length_year,times=Number_years)
  year<-rep(c(1:Number_years),each=length_year)
  df<-as.data.frame(cbind(year,dates,NDVI_TS))
  df_comp<-df
  colnames(df)<-c("Year","Dates","TS")
  window<-EOS-SOS+1
  if (EOS>length_year){ #In case the EOS is greater than the length of the year
    EOS<- EOS - length_year
    df<-subset(df,dates<=EOS | dates>=SOS)
    df<-df[(EOS+1):nrow(df),]
  }else{
    EOS<-EOS
    df<-subset(df,dates<=EOS & dates>=SOS) 
  }
  #STL decomposition
  start_date<-as.Date(df[1,2])
  TS<-ts(df$TS,start=start_date,frequency = window)
  stl_TS<-stl(TS,s.window=7,t.window=window)
  df_stl<-as.data.frame(stl_TS$time.series)
  
  ###ANALYZE TREND COMPONENT
  #Trend
  Trend<-df_stl["trend"]
  Trend1<-as.numeric(Trend$trend)
  df_trend<-as.data.frame(cbind(Trend1,c(1:length(Trend1))))
  mann_test<- as.numeric(MannKendall(Trend1)$sl)
  theil_sen<-theil_sen_regression(Trend1~V2,data=df_trend)
  trend_coef<-theil_sen$coefficients[2]
  #Mosum
  lm_model <- lm(Trend1~V2,data=df_trend)
  residuals <- residuals(lm_model)
  mosum_test <- mosum(residuals, G=0.2)
  pvalue<-mosum_test$cpts.info$p.value
  mosum_result<- sum(pvalue < 0.05)
  
  ###ANALYZE SEASON COMPONENT
  min1 <- as.data.frame(tapply(df$TS, df$Year, min))
  min1$Year<-c(1:nrow(min1))
  colnames(min1)[1]<-"TS"
  
  pvalue_min<- as.numeric(MannKendall(min1$TS)$sl)  #Mann-Kendall test
  theil_sen<-theil_sen_regression(TS~Year,data=min1)
  coef_min<-theil_sen$coefficients[2]
  
  max1 <- as.data.frame(tapply(df$TS, df$Year, max))
  max1$Year<-c(1:nrow(max1))
  colnames(max1)[1]<-"TS"
  
  pvalue_max<- as.numeric(MannKendall(max1$TS)$sl)  #Mann-Kendall test
  theil_sen<-theil_sen_regression(TS~Year,data=max1)
  coef_max<-theil_sen$coefficients[2]
  
  amplitude<-as.data.frame(cbind(max1,min1$TS))
  colnames(amplitude)[c(1,3)]<-c("max","min")
  amplitude$ampl<-amplitude$max-amplitude$min
  
  pvalue_amp<- as.numeric(MannKendall(amplitude$ampl)$sl)  #Mann-Kendall test
  theil_sen<-theil_sen_regression(ampl~Year,data=amplitude)
  coef_amp<-theil_sen$coefficients[2]
  
  #REMAINDER (negative_semivar and standard deviation)
  remains<-as.data.frame(cbind(df$Year,(df_stl$remainder+df_stl$trend)))
  colnames(remains)<-c("year","remainder")
  remains$remainder<-as.numeric(remains$remainder)
  remains$year<-as.numeric(remains$year)
  remains_sd<-aggregate(remainder~year,data=remains,FUN=sd)
  neg_semivar<- remains %>% 
    group_by(year) %>%
    summarize(Neg_semivar=calculate_negative_semivariance(remainder))
  pvalue_sd<- as.numeric(MannKendall(remains_sd$remainder)$sl)  #Mann-Kendall test
  theil_sen<-theil_sen_regression(remainder~year,data=remains_sd)
  coef_sd<-theil_sen$coefficients[2]
  pvalue_semivar<- as.numeric(MannKendall(neg_semivar$Neg_semivar)$sl)  #Mann-Kendall test
  theil_sen<-theil_sen_regression(Neg_semivar~year,data=neg_semivar)
  coef_semivar<-theil_sen$coefficients[2]
  
  #Rolling standard deviation and negative semivariance (COMPLETE YEAR)
  start_date<-as.Date(df_comp[1,2])
  TS<-ts(df_comp$NDVI_TS,start=start_date,frequency = length_year)
  stl_TS_comp<-stl(TS,s.window=7,t.window=window)
  df_stl_comp<-as.data.frame(stl_TS_comp$time.series)
  
  remainder_comp<-as.numeric(df_stl_comp$remainder)+as.numeric(df_stl_comp$trend)
  rolling_sd<-zoo::rollapply(remainder_comp,width=length_year,FUN=sd)
  rolling_semivar<-zoo::rollapply(remainder_comp,width=length_year,FUN=calculate_negative_semivariance)
  
  rolling<-as.data.frame(cbind(c(1:length(rolling_sd)),rolling_sd,rolling_semivar))
  
  pvalue_sd_rolling<- as.numeric(MannKendall(rolling$rolling_sd)$sl)  #Mann-Kendall test
  theil_sen<-theil_sen_regression(rolling_sd~V1,data=rolling)
  coef_sd_rolling<-theil_sen$coefficients[2]
  
  pvalue_semivar_rolling<- as.numeric(MannKendall(rolling$rolling_semivar)$sl)  #Mann-Kendall test
  theil_sen<-theil_sen_regression(rolling_semivar~V1,data=rolling)
  coef_semivar_rolling<-theil_sen$coefficients[2]
  
  results<-as.data.frame(t(c(trend_coef,mann_test,mosum_result,coef_amp,pvalue_amp,coef_sd,pvalue_sd,coef_semivar,pvalue_semivar,coef_sd_rolling,pvalue_sd_rolling,coef_semivar_rolling,pvalue_semivar_rolling)))
  colnames(results)<-c("TREND_theil_sen","TREND_mann_test","MOSUM", "SEASON_ampl_theil_sen", "SEASON_ampl_mann_test","REMAINDER_sd_theil_sen","REMAINDER_sd_mann_test","REMAINDER_negsemivar_theil_sen","REMAINDER_negsemivar_mann_test","SD_rolling_theil_sen","SD_rolling_mann_test","NEGSEMIVAR_rolling_theil_sen","NEGSEMIVAR_rolling_mann_test")
  return(results)
}


