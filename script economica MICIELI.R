###############################################################
###                                                         ###
###   Studente: MICIELI OTTAVIO                             ###
###   Matricola: 214209                                     ###
###                                                         ###
###   Script tesina analisi serie storica di tipo economico ###
###   Esame MODELLI E TECNICHE DI PREVISIONE                ###
###   Professore PIERFRANCESCO PERRI                        ###
###                                                         ###
###############################################################

####                 PACCHETTI UTILIZZATI                  ####
library(fBasics)
library(forecast)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(lmtest)
library(Metrics)
library(moments)
library(seastests)
library(tidyverse)
library(TSA)
library(tseries)
library(TSstudio)

####           FUNZIONI UTILIZZATE CREATE AD HOC           ####

# Funzione per l'esecuzione della valutazione del modello circa le previsioni con approccio
# sliding window

sliding_window<-function(serie, arima_par=c(0,0,0), ntest){
   mat<-matrix(ncol=3)
   for(i in c(seq(from=round(length(serie)/2),to=length(serie)-ntest,by=ntest))){
      training<-window(serie,start=1+i-1000,end=i)
      test<-window(serie,start=i+1,end=i+ntest)
      model<-Arima(training,order=arima_par)
      forec<-forecast(model,h=ntest)
      mat<-rbind(mat,c(rmse(test,forec$mean),mae(test,forec$mean),mape(test,forec$mean)))
   }
   return(data.frame(RMSE=mean(mat[-1,1]),MAE=mean(mat[-1,2]),MAPE=mean(mat[-1,3])))
}

####              FUNZIONI GRAFICHE UTILIZZATE             ####

# layout grafico per serie storica

plot_ts<-function(serie,subtitle=NULL){
   data<-data.frame(t=1:length(serie),value=serie)
   ggplot(data,aes(x=t,y=value))+
      geom_line(color="#F4B183",cex=1)+
      labs(title="Serie storica di tipo economico", subtitle=subtitle, y=bquote(y[t]))+
      theme_minimal()+
      theme(
         plot.title=element_text(size=20,hjust=0.5),
         plot.subtitle=element_text(size=12,hjust=0.5,color="grey30"),
         axis.title=element_text(size=15,hjust=0.5)
      )
}

plot_norm<-function(serie,titolo1,titolo2,row){
   data<-data.frame(t=1:length(serie),value=serie)
   media=round(mean(data$value),4)
   hist<-ggplot(data=data,aes(x=value))+
      geom_histogram(aes(y=..density..),bins=50,size=0.1,color="white",fill="#F4B183")+
      geom_vline(aes(xintercept=mean(value)),cex=1.5)+
      stat_function(fun=dnorm,args=list(mean=mean(data$value),sd=sd(data$value)),cex=1,color="gray25")+
      labs(title=titolo1,x=bquote(y[t]))+
      annotate(geom="text",x=mean(data$value)+15,y=0.0015,label=bquote("Media: "~.(media)),angle=270,size=3)+
      theme_bw()+
      theme(
         plot.title=element_text(size=15,hjust=0.5),
         axis.title=element_text(size=10),
         panel.grid=element_line(linetype="dashed",color="grey70")
      )
   
   qq<-ggplot(data=data,aes(sample=value))+geom_qq_line(cex=1)+geom_qq(color="#F4B183")+
      theme_bw()+
      labs(title=titolo2,x="Teoriche",y="Osservate")+
      theme(
         plot.title=element_text(size=15,hjust=0.5),
         axis.title=element_text(size=10),
         panel.grid=element_line(linetype="dashed",color="grey80")
      )
   
   grid.arrange(hist,qq,nrow=row)
}

plot_corr<-function(serie,lag.max,ci,horiz,subtitle=NULL){
   data<-data.frame(acf=acf(serie,lag=lag.max,plot=F)$acf,pacf=pacf(serie,lag=lag.max,plot=F)$acf,lag=1:lag.max)
   acf_plot<-ggplot(data,aes(x=lag,y=acf))+
               geom_col(fill="#F4B183",width=0.7)+
               geom_hline(yintercept=qnorm((1-0.95)/2)/sqrt(length(serie)),color="grey50",linetype="dashed",cex=1)+
               geom_hline(yintercept=-qnorm((1-0.95)/2)/sqrt(length(serie)),color="grey50",linetype="dashed",cex=1)+
               geom_hline(yintercept=0,color="black",cex=0.8)+
               scale_x_continuous(breaks=seq(0,max(data$lag),5))+
               scale_y_continuous(limits=c(min(data$acf,data$pacf),1))+
               theme_bw()+
               labs(title="Estimated Autocorrelation function",subtitle=subtitle,x="Lag",y=bquote(hat(rho[k])))+
               theme(
                  plot.title=element_text(size=15,hjust=0.5),
                  plot.subtitle=element_text(size=10,hjust=0.5,color="grey50"),
                  axis.title=element_text(size=10),
                  axis.text=element_text(size=7),
                  panel.grid=element_line(linetype="dashed",color="grey90")
               )
   pacf_plot<-ggplot(data,aes(x=lag,y=pacf))+
               geom_col(fill="#F4B183",width=0.7)+
               geom_hline(yintercept=qnorm((1-0.95)/2)/sqrt(length(serie)),color="grey50",linetype="dashed",cex=1)+
               geom_hline(yintercept=-qnorm((1-0.95)/2)/sqrt(length(serie)),color="grey50",linetype="dashed",cex=1)+
               geom_hline(yintercept=0,color="black",cex=0.8)+
               scale_x_continuous(breaks=seq(0,max(data$lag),5))+
               scale_y_continuous(limits=c(min(data$acf,data$pacf),1))+
               theme_bw()+
               labs(title="Estimated Partial Autocorrelation function",subtitle=subtitle,x="Lag",y=bquote(hat(pi[k])))+
               theme(
                  plot.title=element_text(size=15,hjust=0.5),
                  plot.subtitle=element_text(size=10,hjust=0.5,color="grey50"),
                  axis.title=element_text(size=10),
                  axis.text=element_text(size=7),
                  panel.grid=element_line(linetype="dashed",color="grey90")
               )
   if(horiz==T){grid.arrange(acf_plot,pacf_plot,nrow=1)}
   else{grid.arrange(acf_plot,pacf_plot,nrow=2)}
}

plot_lag<-function(serie,lag.max){
   gglagplot(serie,lags=lag.max,do.lines=F,colour=F,diag.col="grey50")+
      theme_bw()+
      labs(title="Lag plot")+
      theme(
         plot.title=element_text(size=25,hjust=0.5),
         axis.title=element_blank(),
         axis.text=element_blank(),
         panel.grid=element_blank(),
         axis.ticks=element_blank(),
         strip.background=element_rect(fill="#F4B183"),
         strip.text=element_text(color="white")
      )
}

plot_boxSeas<-function(serie,freq,subtitle=NULL){
   data<-data.frame(value=as.numeric(serie),lag=factor(rep(seq(from=1,to=freq,by=1),length.out=length(serie))))
   ggplot(data,aes(x=lag,y=value))+
      geom_boxplot(outlier.colour="#F4B183",outlier.size=2)+
      scale_x_discrete(breaks=seq(from=1,to=freq,by=2))+
      theme_bw()+
      labs(title="Boxplot stagionale",subtitle=subtitle,x="Period",y=NULL)+
      theme(
         plot.title=element_text(size=12,hjust=0.5),
         plot.subtitle=element_text(size=10,hjust=0.5,color="grey50"),
         axis.title=element_text(size=8),
         axis.text=element_text(size=5),
         panel.grid=element_line(linetype="dashed",color="grey90")
      )
}

plot_tsfit<-function(serie,fit,subtitle=NULL){
   data<-data.frame(t=1:length(serie),value=serie,fitted=fit)
   ggplot(data)+
      geom_line(aes(x=t,y=fitted,colour="Serie teorica"),cex=1)+
      geom_line(aes(x=t,y=value,colour="Serie originale"),cex=1,alpha=0.7)+
      labs(title="Valori osservati e valori teorici", subtitle=subtitle, y=bquote(y[t]~" / "~hat(y[t])))+
      scale_colour_manual("",breaks=c("Serie originale","Serie teorica"),
                          values=c("Serie originale"="grey50","Serie teorica"="#F4B183"))+
      theme_minimal()+
      theme(
         plot.title=element_text(size=20,hjust=0.5),
         plot.subtitle=element_text(size=12,hjust=0.5,color="grey30"),
         axis.title=element_text(size=15,hjust=0.5),
         legend.position=c(.95,.95),
         legend.justification=c("right","top"),
         legend.box.just="right",
         legend.margin=margin(6,6,6,6)
      )
}

plot_res<-function(res,subtitle=NULL){
   data<-data.frame(value=res,lag=1:length(res))
   ggplot(data,aes(x=lag,y=res))+
      geom_point()+
      theme_bw()+
      labs(title="Scatterplot dei residui",subtitle=subtitle,x="t",y=bquote(hat(epsilon[t])))+
      theme(
         plot.title=element_text(size=15,hjust=0.5),
         plot.subtitle=element_text(size=12,hjust=0.5,color="grey50"),
         axis.title=element_text(size=15),
         axis.text=element_text(size=10),
         panel.grid=element_line(linetype="dashed",color="grey90")
      )
}

plot_rescorr<-function(res,lag.max,subtitle=NULL){
   data<-data.frame(acf=acf(res,lag=lag.max,plot=F)$acf,lag=1:lag.max)
   ggplot(data,aes(x=lag,y=acf))+
      geom_col(fill="#F4B183",width=0.7)+
      geom_hline(yintercept=qnorm((1-0.95)/2)/sqrt(length(res)),color="grey50",linetype="dashed",cex=1)+
      geom_hline(yintercept=-qnorm((1-0.95)/2)/sqrt(length(res)),color="grey50",linetype="dashed",cex=1)+
      geom_hline(yintercept=0,color="black",cex=0.8)+
      scale_x_continuous(breaks=seq(0,max(data$lag),5))+
      scale_y_continuous(limits=c(min(data$acf,-0.3),max(data$acf,0.3)))+
      theme_bw()+
      labs(title="Estimated Autocorrelation function",subtitle=subtitle,x="Lag",y=bquote(hat(rho[k])))+
      theme(
         plot.title=element_text(size=20,hjust=0.5),
         plot.subtitle=element_text(size=15,hjust=0.5,color="grey50"),
         axis.title=element_text(size=15),
         axis.text=element_text(size=10),
         panel.grid=element_line(linetype="dashed",color="grey90")
      )
}

plot_resnorm<-function(res,subtitle=NULL,row=T){
   data<-data.frame(t=1:length(res),value=res)
   media=round(mean(data$value),4)
   hist<-ggplot(data=data,aes(x=value))+
      geom_histogram(aes(y=..density..),bins=50,size=0.1,color="white",fill="#F4B183")+
      geom_vline(aes(xintercept=mean(value)),cex=1.5)+
      stat_function(fun=dnorm,args=list(mean=mean(data$value),sd=sd(data$value)),cex=1,color="gray25")+
      labs(title="Istogramma dei residui",subtitle=subtitle,x=bquote(hat(epsilon[t])))+
      annotate(geom="text",x=mean(data$value)+5,y=0.35,label=bquote("Media: "~.(media)),size=4)+
      theme_bw()+
      theme(
         plot.title=element_text(size=20,hjust=0.5),
         plot.subtitle=element_text(size=10,hjust=0.5,color="grey50"),
         axis.title=element_text(size=10),
         panel.grid=element_line(linetype="dashed",color="grey70")
      )
   
   qq<-ggplot(data=data,aes(sample=value))+geom_qq_line(cex=1)+geom_qq(color="#F4B183")+
      theme_bw()+
      labs(title="QQ plot",subtitle=subtitle,x="Teoriche",y="Osservate")+
      theme(
         plot.title=element_text(size=15,hjust=0.5),
         plot.subtitle=element_text(size=10,hjust=0.5,color="grey50"),
         axis.title=element_text(size=10),
         panel.grid=element_line(linetype="dashed",color="grey80")
      )
   
   if(row==T){grid.arrange(hist,qq,nrow=1)} else {grid.arrange(hist,qq,nrow=2)}
}

####                LETTURA DEI DATI DA FILE               ####

setwd("C:/Users/ottav/Desktop/Progetto Perri/Serie economica")
load("dati_febbraio_2021.csv")
ts<-ts(Y); rm(Y)

####          ANALISI DESCRITTIVA SERIE STORICA            ####

format(basicStats(ts),scientific=F)
plot_ts(ts,"Valori osservati") 
plot_norm(ts,"Istogramma delle frequenze della serie","QQ plot serie storica",1)

ksnormTest(ts)
shapiroTest(ts)
jarqueberaTest(ts)
dagoTest(ts)

agostino.test(ts)
anscombe.test(ts)

####    ANALISI STAZIONARIETà E STAGIONALITà DELLA SERIE   ####

# Verifica stazionarietà serie storica originale
adf.test(ts,alternative = "stationary")
adf.test(ts,alternative = "explosive")
kpss.test(ts)
pp.test(ts, alternative = "stationary")
pp.test(ts, alternative = "explosive")

plot_corr(ts,lag.max=50,ci=0.95,horiz=T)

#Individuazione dell'ordine della differenziazione
ndiffs(ts)
ts_d1<-diff(ts)

adf.test(ts_d1,alternative = "stationary")
adf.test(ts_d1,alternative = "explosive")
kpss.test(ts_d1)
pp.test(ts_d1, alternative = "stationary")
pp.test(ts_d1, alternative = "explosive")

#Individuazione della componente stagionale all'interno della serie storica

plot_ts(ts_d1,bquote("Serie storica differenziata"~nabla~Y[t]))
plot_corr(ts_d1,lag.max=15,ci=0.95,horiz=F,bquote("Serie storica differenziata"~nabla~Y[t]))
plot_lag(ts_d1,48)

isSeasonal(ts_d1, freq=52)
isSeasonal(ts_d1, freq=12)
isSeasonal(ts_d1, freq=7)
isSeasonal(ts_d1, freq=4)
isSeasonal(ts_d1, freq=3)

isSeasonal(ts_d1, freq=52, test="kw")
isSeasonal(ts_d1, freq=12, test="kw")
isSeasonal(ts_d1, freq=7, test="kw")
isSeasonal(ts_d1, freq=4, test="kw")
isSeasonal(ts_d1, freq=3, test="kw")

grid.arrange(plot_boxSeas(ts_d1,52,"Frequenza 52 (serie storica settimanale)"),
             plot_boxSeas(ts_d1,12,"Frequenza 12 (serie storica mensile)"),
             nrow=2)

grid.arrange(plot_boxSeas(ts_d1,7,"Frequenza 7 (serie storica giornaliera, periodo di riferimento:settimana)"),
             plot_boxSeas(ts_d1,4,"Frequenza 4 (serie storica trimestrale)"),#messo solo per mantenere la forma
             nrow=2)

grid.arrange(plot_boxSeas(ts_d1,4,"Frequenza 4 (serie storica trimestrale)"),
             plot_boxSeas(ts_d1,3,"Frequenza 3 (serie storica quadrimestrale)"), 
             nrow=1)


####          INDIVIDUAZIONE DEL MIGLIOR MODELLO           ####

auto.arima(ts)
auto.arima(ts_d1)

# Ricerchiamo manualmente il processo generatore dei dati
p<-0:5
q<-0:5
par_mod<-expand.grid(p,q)
names(par_mod)<-c("p", "q")
par_mod$d<-1 

arima_search<-lapply(1:nrow(par_mod), function(i){
   md<-NULL
   md<-arima(ts, order=c(par_mod$p[i], 1, par_mod$q[i]))
   results<-data.frame(p=par_mod$p[i], d=1, q=par_mod$q[i], AIC=md$aic, BIC=BIC(md))
}) %>% bind_rows() %>% arrange(AIC)
head(arima_search)

# stessa funzione ma viene applicata la funzione Arima in forecast
Arima_search<-lapply(1:nrow(par_mod), function(i){
   md<-NULL
   md<-Arima(ts, order=c(par_mod$p[i], 1, par_mod$q[i]))
   results<-data.frame(p=par_mod$p[i], d=1, q=par_mod$q[i], AIC=md$aic, BIC=BIC(md))
}) %>% bind_rows() %>% arrange(AIC)
head(Arima_search)

ARIMA_2.1.0<-Arima(ts,order=c(2,1,0)); summary(ARIMA_2.1.0); coeftest(ARIMA_2.1.0); coefci(ARIMA_2.1.0)
ARIMA_1.1.2<-Arima(ts,order=c(1,1,2)); summary(ARIMA_1.1.2); coeftest(ARIMA_1.1.2); coefci(ARIMA_1.1.2)

plot_tsfit(ts,ARIMA_2.1.0$fitted,"Modello ARIMA(2,1,0)")
plot_tsfit(ts,ARIMA_1.1.2$fitted,"Modello ARIMA(1,1,2)")

####                  ANALISI DEI RESIDUI                  ####

# Analisi residui Modello ARIMA(2,1,0)

res_ARIMA_2.1.0<-residuals(ARIMA_2.1.0)
plot_res(res_ARIMA_2.1.0, "Modello ARIMA(2,1,0)")
plot_rescorr(res_ARIMA_2.1.0,50,"Modello ARIMA(2,1,0)")
plot_resnorm(res_ARIMA_2.1.0,"Modello ARIMA(2,1,0)",row=F)

Box.test(residuals(ARIMA_2.1.0),lag=50,type="Box-Pierce")
Box.test(residuals(ARIMA_2.1.0),lag=50,type="Ljung-Box")

ksnormTest(res_ARIMA_2.1.0)
shapiroTest(res_ARIMA_2.1.0)
jarqueberaTest(res_ARIMA_2.1.0)
dagoTest(res_ARIMA_2.1.0)

McLeod.Li.test(object=ARIMA_2.1.0) #no eteroschedasticità nei residui

tsdiag(ARIMA_2.1.0,gof.lag=50)

# Analisi residui Modello ARIMA(1,1,2)

res_ARIMA_1.1.2<-residuals(ARIMA_1.1.2)
plot_res(res_ARIMA_1.1.2, "Modello ARIMA(1,1,2)")
plot_rescorr(res_ARIMA_1.1.2,50,"Modello ARIMA(1,1,2)")
plot_resnorm(res_ARIMA_1.1.2,"Modello ARIMA(1,1,2)",row=F)

Box.test(residuals(ARIMA_1.1.2),lag=50,type="Box-Pierce")
Box.test(residuals(ARIMA_1.1.2),lag=50,type="Ljung-Box")

ksnormTest(res_ARIMA_1.1.2)
shapiroTest(res_ARIMA_1.1.2)
jarqueberaTest(res_ARIMA_1.1.2)
dagoTest(res_ARIMA_1.1.2)

McLeod.Li.test(object=ARIMA_1.1.2) #no eteroschedasticità nei residui

tsdiag(ARIMA_1.1.2,gof.lag=50)

####         ANALISI SULL'ACCURATEZZA DI PREVISIONE        ####

#la serie ricopre il ruolo sia di training test che di test set
forecast::accuracy(ARIMA_2.1.0)
forecast::accuracy(ARIMA_1.1.2)

#vengono escluse le ultime 5 osservazioni
train<-window(ts,start=1,end=1996)
test_set<-window(ts,start=1997,end=length(ts))
model_2.1.0<-Arima(train,order=c(2,1,0))
fore_2.1.0<-forecast(model_2.1.0,h=length(test_set))
forecast::accuracy(fore_2.1.0,test_set)

model_1.1.2<-Arima(train,order=c(1,1,2))
fore_1.1.2<-forecast(model_1.1.2,h=length(test_set))
forecast::accuracy(fore_1.1.2,test_set)

#utilizziamo la tecnica dello sliding window
sliding_window(ts,arima_par=c(2,1,0),5)
sliding_window(ts,arima_par=c(1,1,2),5)

#Analizziamo per quale ampiezza campionaria il modello ARIMA(2,1,0) performa meglio dell'ARIMA(2,1,1)
matrice<-matrix(ncol=6)
for(j in 1:50){ #impiega un po' di tempo a calcolare
   A2.1.0<-sliding_window(ts,arima_par=c(2,1,0),j)
   A1.1.2<-sliding_window(ts,arima_par=c(1,1,2),j)
   matrice<-rbind(matrice,c(A2.1.0$RMSE,A2.1.0$MAE,A2.1.0$MAPE,A1.1.2$RMSE,A1.1.2$MAE,A1.1.2$MAPE))
}
matrice<-as.data.frame(matrice[-1,])
result<-data.frame(RMSE=mean(matrice[,1]-matrice[,4]),MAE=mean(matrice[,2]-matrice[,5]),
                   MAPE=mean(matrice[,3]-matrice[,6])); result
RMSE<-0
MAE<-0
MAPE<-0
for(h in 1:nrow(matrice)){
   if(matrice[h,1]<matrice[h,4]){RMSE<-RMSE+1}
   if(matrice[h,2]<matrice[h,5]){MAE<-MAE+1}
   if(matrice[h,3]<matrice[h,6]){MAPE<-MAPE+1}
}

####         ANALISI SPERIMENTALI PER SCELTA MODELLO       ####

theo_acf<-as.numeric(ARMAacf(ar=c(ARIMA_2.1.0$coef[1],ARIMA_2.1.0$coef[2]),lag.max=50))
theo_acf1<-as.numeric(ARMAacf(ar=c(ARIMA_1.1.2$coef[1]),ma=c(ARIMA_1.1.2$coef[2],ARIMA_1.1.2$coef[3]),lag.max=50))
acf<-Acf(ts_d1,lag.max=50,plot=F)
data.frame(Mod_ARIMA_2.1.0=mean(abs(acf$acf-theo_acf)),Mod_ARIMA_1.1.2=mean(abs(acf$acf-theo_acf1)))

theo_pacf<-ARMAacf(ar=c(ARIMA_2.1.0$coef[1],ARIMA_2.1.0$coef[2]),lag.max=50,pacf=T)
theo_pacf1<-ARMAacf(ar=c(ARIMA_1.1.2$coef[1]),ma=c(ARIMA_1.1.2$coef[2],ARIMA_1.1.2$coef[3]),lag.max=50,pacf=T)
pacf<-Acf(ts_d1,lag.max=50,plot=F,type="partial")
data.frame(Mod_ARIMA_2.1.0=mean(abs(pacf$acf-theo_pacf)),Mod_ARIMA_1.1.2=mean(abs(pacf$acf-theo_pacf1)))

prev_ARIMA_1.1.2<-forecast(ARIMA_1.1.2,h=50)

autoplot(window(ts,start=1950),cex=1)+autolayer(prev_ARIMA_1.1.2,cex=1,colour="#F4B183")+
   labs(title="Previsioni su modello ARIMA(1,1,2)", y=bquote(y[t]))+
   theme_minimal()+
   theme(
      plot.title=element_text(size=20,hjust=0.5),
      plot.subtitle=element_text(size=12,hjust=0.5,color="grey30"),
      axis.title=element_text(size=15,hjust=0.5),
      panel.grid=element_line(linetype="dashed",color="grey90")
   )