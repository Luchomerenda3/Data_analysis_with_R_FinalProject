setwd("/media/luciano/LUCHO/UNCuyo-FCEN/Análisis de datos con R/Informe_final/")

source("functions.R")
library(lubridate)
library(tidyverse)
library(tsbox)
library(gganimate)
library(gifski)
library(glue)
library(broom)
#library(forecast)
#library(ggTimeSeries)
#library(TSstudio)
#library(sugrrants)

#USEFUL FACTS
#The POSIXct class stores date/time values as the number of seconds since January 1, 1970
#When using ts_lag() we have to use strings as delta_t! in non-regular data,
#it works with numeric integer but only on regular data.

######## 90x90 pixel Boxes Data ######### We wont use it for the moment.

#search & read 36 events 90x90_pix boxes computed data
#dat_90x90<-list.files(pattern = "^swp45_(.*)csv$", full.names = TRUE)
#evs_90x90<-map(dat_90x90,read_csv)

#converting original date format in every tibble
#for (i in seq_along(evs_90x90))
#{
#  evs_90x90[[i]]$Date_obs<-as_datetime(evs_90x90[[i]]$Date_obs)
#}


######## 30x30 pixel Boxes Data #########

#search & read 36 events 30x30_pix boxes computed data
dat_30x30<-list.files(pattern = "^swp_(.*)csv$", full.names = TRUE)
evs_30x30<-map(dat_30x30,read_csv)

#converting original date format in every tibble 
#and normalize it (in days units) to the event time
for (i in seq_along(evs_30x30))
  {
   evs_30x30[[i]]$Date_obs<-as_datetime(evs_30x30[[i]]$Date_obs)
   evs_30x30[[i]]<-rename(evs_30x30[[i]],time=Date_obs)
  }

#We need to convert each parameter to a time series object
#first lets create 36 length list for each parameter that will contain
#a time series of the evolution of that parameter for each event!
usfz<-list(36)
mcdz<-list(36)
mchz<-list(36)
tchz<-list(36)
absnchz<-list(36)

#Then we asign them values 
for (i in c(1:36))
   {
    #First USFz (unsigned total flux)
    usfz[[i]]<- evs_30x30[[i]] %>%
                select(time,USFz) %>%
                rename(value=USFz)
    #Mean current density
    mcdz[[i]]<- evs_30x30[[i]] %>%
                select(time,MCDz) %>%
                rename(value=MCDz) 
    #Mean Current Helicity
    mchz[[i]]<- evs_30x30[[i]] %>%
                select(time,MCHz) %>%
                rename(value=MCHz)
    #Total Unsigned Current Helicity
    tchz[[i]]<- evs_30x30[[i]] %>%
                select(time,TCHz) %>%
                rename(value=TCHz)
    #Absolute Value of the Net Current Helicity
    absnchz[[i]]<- evs_30x30[[i]] %>%
                   select(time,ABSNCHz) %>%
                   rename(value=ABSNCHz) 
   }

#Leemos data sobre CMES identificadas para la región!
#Corregimos la fecha y la transformamos como objeto datetime
cmes_fer<-read_csv("cmeinfo_fer.csv") %>%
  mutate(Date=as.Date(Date,'%d/%m/%Y')) %>%
  mutate(Date=str_c(as.character(Date),' ',Time),Time=NULL) %>%
  mutate(Date=as_datetime(Date)) %>%
  rename(time=Date)
cmes_fsel<-read_csv("5rotcme_locationsFV.csv") %>%
  mutate(Date=str_c(as.character(Date),' ',Hour),Hour=NULL) %>%
  mutate(Date=as_datetime(Date)) %>%
  rename(time=Date)

#De todas las cmes identificadas de la región seleccionamos las que podemos 
cme_fer_fv<-inner_join(cmes_fsel,cmes_fer) %>%
            rename(Angular_Width=`Angular Width`,
                   Central_Position_Angle=`Central Pos. Angle`,
                   Linear_speed=`Linear Speed`,
                   Kinetic_Energy=`Kinetic Energy`,
                   Potencial_Energy=`Potential Energy`)

#Now we create Morphologic & Dynamics Time-series
angwid<-cme_fer_fv %>%
        select(time,Angular_Width)
cpa<-cme_fer_fv %>%
     select(time,Central_Position_Angle)
lisp<-cme_fer_fv %>%
      select(time,Linear_speed)
mas<-cme_fer_fv %>%
     select(time,Mass) 
poe<-cme_fer_fv %>%
     select(time,Potencial_Energy)
kie<-cme_fer_fv %>%
     select(time,Kinetic_Energy)
erup_times<-cme_fer_fv$time

#Ya tenemos los dos time series para analizar 
#Los time-series morfologicos y/o dinamicos,
#y las list de time-series de parametros fotosfericos
#de cada evento.

#Lets make a delta_t vector with del_t values
#Determiné tomar un máximo de 3 días antes del evento, en intervalos de 1 hora
#para disminuir lo más posible la cantidad de datos sin perder cadencia de las
#observaciones que es de 12 min. es decir redujimos la cadencia en 5 veces.

delta_t<-snum2str(seq(0,3600*24*3,by=3600),"sec")
delta_tn<-seq(0,3600*24*3,by=3600)
#testing example
#b<-sw_md(angwid,usfz,delta_t,erup_times)
#traceback()
##ITS WORKING!!!

#Ok so we have n=5 photosferic quantities & m=6 morpho/dynamic parameters.
#This mean that mean can obtain n*m =30 ways of combinations for this parameters
#each combination gives us 73 linear regressions models, 73 correlation factors 
#per each delta_t and 73 tibbles with the final photospheric parameter vs 
#morpho/dynamic for each event!. Thats a lot of data to analize! even more!
#we have the same photospheric quantities computed for 90x90 pixels boxes i.e.
#we can make 30 more combination in the same way as before Thats HUGE amount of data.

#To simplify we'll choose 2 photospheric parameters and 2 morpho/dynamic quantities 
#The choice is based on experimental evidence of correlation between magnetic flux and 
#Helicity with the eruption periods.

#We chose usfz & tchz as space_weather properties. and ang_wid & kie as morpho/dynamic parameters.
#4 posibles ways 

usfz_vs_angwid<-sw_md(angwid,usfz,delta_t,erup_times,delta_tn)
usfz_vs_kie<-sw_md(kie,usfz,delta_t,erup_times,delta_tn)

tchz_vs_angwid<-sw_md(angwid,tchz,delta_t,erup_times,delta_tn)
tchz_vs_kie<-sw_md(kie,tchz,delta_t,erup_times,delta_tn)


###BELOW THIS LINE JUST FOR PLOTTING PURPOSES!

##Correlation factor evolution plot for all combinations.

#usfz_vs_angwid
usfzangwid_corfac<-ggplot(usfz_vs_angwid$Corr_factor)+
                   xlab("Delta t [h]")+
                   ylab("Correlation factor \"r\"")+
                   theme(panel.background = element_rect(fill = "white", colour = "black"),
                   panel.grid.major = element_line(colour = "gray90"),
                   panel.grid.minor = element_line(colour = "gray90"))+
                   geom_point(aes(delta_tnum, value ),size=0.7)+
                   scale_x_continuous(breaks =seq(0,3600*24*3,by=4*3600),
                                      labels=snum2str(seq(0,3600*24*3,by=4*3600),"hour"))+
                   geom_line(aes(delta_tnum, value))
usfzangwid_corfac
ggsave("usfzangwid_corfac.png",plot=usfzangwid_corfac ,dpi="retina")



#usfz_vs_kie
usfzkie_corfac<-ggplot(usfz_vs_kie$Corr_factor)+
                xlab("Delta t [h]")+
                ylab("Correlation factor \"r\"")+
                theme(panel.background = element_rect(fill = "white", colour = "black"),
                panel.grid.major = element_line(colour = "gray90"),
                panel.grid.minor = element_line(colour = "gray90"))+
                geom_point(aes(delta_tnum, value ),size=0.7)+
                scale_x_continuous(breaks =seq(0,3600*24*3,by=4*3600),
                                   labels=snum2str(seq(0,3600*24*3,by=4*3600),"hour"))+
                geom_line(aes(delta_tnum, value))
usfzkie_corfac
ggsave("usfzkie_corfac.png",plot=usfzkie_corfac ,dpi="retina")

#tchz_vs_angwid
tchzangwid_corfac<-ggplot(tchz_vs_angwid$Corr_factor)+
                   xlab("Delta t [h]")+
                   ylab("Correlation factor \"r\"")+
                   theme(panel.background = element_rect(fill = "white", colour = "black"),
                   panel.grid.major = element_line(colour = "gray90"),
                   panel.grid.minor = element_line(colour = "gray90"))+
                   geom_point(aes(delta_tnum, value ),size=0.7)+
                   scale_x_continuous(breaks =seq(0,3600*24*3,by=4*3600),
                                      labels=snum2str(seq(0,3600*24*3,by=4*3600),"hour"))+
                   geom_line(aes(delta_tnum, value))
tchzangwid_corfac
ggsave("tchzangwid_corfac.png",plot=tchzangwid_corfac ,dpi="retina")

#tchz_vs_kie
tchzkie_corfac<-ggplot(tchz_vs_kie$Corr_factor)+
                xlab("Delta t [h]")+
                ylab("Correlation factor \"r\"")+
                theme(panel.background = element_rect(fill = "white", colour = "black"),
                panel.grid.major = element_line(colour = "gray90"),
                panel.grid.minor = element_line(colour = "gray90"))+
                geom_point(aes(delta_tnum, value ),size=0.7)+
                scale_x_continuous(breaks =seq(0,3600*24*3,by=4*3600),
                                   labels=snum2str(seq(0,3600*24*3,by=4*3600),"hour"))+
                geom_line(aes(delta_tnum, value))
tchzkie_corfac
ggsave("tchzkie_corfac.png",plot=tchzkie_corfac ,dpi="retina")

#MAke a color for each event.
ev_col<-as.vector(rainbow(36,alpha = 0.7))
ev_label<-as.character(c(1:36))

#Scatter_plots for swpvsmdp
#"We want to animate then since plotting 4x73 scatterplots is useless"
#We want to see the evolution of the scatter plots!
#Thus in order to do it we have to tidy the data.

#usfz_vs_angwid
usfzvsangwid<-matrix(nrow=36*73,ncol=4)
colnames(usfzvsangwid)<-c("time","Angular_width","usfz","delta_t")
usfzvsangwid<-as_tibble(usfzvsangwid)
usfzvsangwid
for (i in c(1:73)){
  usfzvsangwid[]
}

usfzangwid<-ggplot()+xlab("Unsigned Magnetic Flux [Mx]")+
                     ylab("CME's angular width [deg]")+
                     theme(panel.background = element_rect(fill = "white", colour = "black"),
                     panel.grid.major = element_line(colour = "gray90"),
                     panel.grid.minor = element_line(colour = "gray90"))+
                     geom_point(aes(swp,Angular_Width),size=0.7)+
                     geom_smooth(aes(swp,Angular_Width),size=0.7,se=F)



#####################################################################################################

#SHARP parameters computed for entire region!!!
#no correlation... yet. Just for testing purposes 
sh_dat<-list.files(pattern = "^all_swp_(.*)csv$", full.names = TRUE)
sharp<-map(sh_dat,read_csv)

for (i in seq_along(sharp))
{
  sharp[[i]]$T_REC<-as_datetime(sharp[[i]]$T_REC)
}

#Leemos data sobre CMES identificadas para la región!
#Corregimos la fecha y la transformamos como objeto datetime
cmes_fer<-read_csv("cmeinfo_fer.csv") %>%
  mutate(Date=as.Date(Date,'%d/%m/%Y')) %>%
  mutate(Date=str_c(as.character(Date),' ',Time),Time=NULL) %>%
  mutate(Date=as_datetime(Date))
cmes_fsel<-read_csv("5rotcme_locationsFV.csv") %>%
  mutate(Date=str_c(as.character(Date),' ',Hour),Hour=NULL) %>%
  mutate(Date=as_datetime(Date))

#De todas las cmes identificadas de la región seleccionamos las que podemos 
cme_fer_fv<-inner_join(cmes_fsel,cmes_fer)

##Ploteamos para diferentes parametros

#Unsigned Flux 
gg_usflux<-ggplot()+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "gray90"))+
  xlab("Date")+
  ggtitle("USFLUX vs  CME's Angular Widths for 5 Rot")
for (i in c(1:5))
{
  gg_usflux<-gg_usflux+
    geom_line(data = sharp[[i]], aes(T_REC,USFLUX),colour="black",size=.25)
  # assign each one data from each event
}

gg_usflux<-gg_usflux +
           geom_point(data=cme_fer_fv,aes(Date,`Angular Width`*3e+20), colour="Red",size=0.8)+
           scale_y_continuous(sec.axis = sec_axis(~./3e+20, name = "Angular width [Deg]"))
plot(gg_usflux)
ggsave("usflux_vs_Angwid_5rot.png",dpi="retina")#,width= 100,height = 50,units = "mm")

#Mean Current Density 
gg_mjzd<-ggplot()+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "gray90"))+
  xlab("Date")+
  ggtitle("MEANJZD vs  CME's Angular Widths for 5 Rot")
for (i in c(1:5))
{
  gg_mjzd<-gg_mjzd+
    geom_line(data = sharp[[i]], aes(T_REC,MEANJZD),colour="black",size=.25)
  # assign each one data from each event
}
gg_mjzd
gg_mjzd<-gg_mjzd +
  geom_point(data=cmes_fer,aes(Date,`Angular Width`/90), colour="Red",size=0.8)+
  scale_y_continuous(sec.axis = sec_axis(~.*90, name = "Angular width [Deg]"))
plot(gg_mjzd)
ggsave("meanjzd_vs_Angwid_5rot.png",dpi="retina")#,width= 100,height = 50,units = "mm")

#Mean Current Helicity
gg_mjzh<-ggplot()+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "gray90"))+
  xlab("Date")+
  ggtitle("MEANJZH vs  CME's Angular Widths for 5 Rot")
for (i in c(1:5))
{
  gg_mjzh<-gg_mjzh+
    geom_line(data = sharp[[i]], aes(T_REC,MEANJZH),colour="black",size=.25)
  # assign each one data from each event
}
gg_mjzh
gg_mjzh<-gg_mjzh +
  geom_point(data=cmes_fer,aes(Date,`Angular Width`/9000), colour="Red",size=0.8)+
  scale_y_continuous(sec.axis = sec_axis(~.*9000, name = "Angular width [Deg]"))
plot(gg_mjzh)
ggsave("meanjzh_vs_Angwid_5rot.png",dpi="retina")#,width= 100,height = 50,units = "mm")


#Total unsigned Current Helicity
gg_totusjh<-ggplot()+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "gray90"))+
  xlab("Date")+
  ggtitle("MEANJZH vs  CME's Angular Widths for 5 Rot")
for (i in c(1:5))
{
  gg_totusjh<-gg_totusjh+
    geom_line(data = sharp[[i]], aes(T_REC,TOTUSJH),colour="black",size=.25)
}
gg_totusjh
gg_totusjh<-gg_totusjh +
  geom_point(data=cmes_fer,aes(Date,`Angular Width`*25), colour="Red",size=0.8)+
  scale_y_continuous(sec.axis = sec_axis(~./25, name = "Angular width [Deg]"))
plot(gg_totusjh)
ggsave("totusjh_vs_Angwid_5rot.png",dpi="retina")


###PLotting of all events on one graph for each quantity
####Analysis for each computed box that corresponds to a single CME in the cmds_fsel ^ cme_fer_fv.
######## 30x30 pixel Boxes Data #########

#search & read 36 events 30x30_pix boxes computed data
dat_30x30<-list.files(pattern = "^swp_(.*)csv$", full.names = TRUE)
evs_30x30<-map(dat_30x30,read_csv)

#event dates we need them to normalize dates in graphs!
ev_dat<-read_csv("5rotcme_locationsFV.csv") %>%
  mutate(Date=str_c(as.character(Date),' ',Hour),Hour=NULL) %>%
  mutate(Date=as_datetime(Date))
summary((ev_dat$Date[1]-ev_dat$Date[2])/86400)

#converting original date format in every tibble 
#and normalize it (in days units) to the event time
for (i in seq_along(evs_30x30))
{
  evs_30x30[[i]]$Date_obs<-as.duration(as_datetime(evs_30x30[[i]]$Date_obs)-ev_dat$Date[i])/duration(86400,units="seconds")
}

#We create a ggplot per swp quantity
#First we create a color palette for each event
ev_col<-as.vector(rainbow(36,alpha = 0.7))
ev_label<-as.character(c(1:36))
#MCDZ
gg_30x30_MCDZ<-ggplot()+
  scale_x_continuous(limits=c(-3,1),breaks=seq(-3,1,by=0.5))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "gray90"))+
  geom_vline(xintercept = 0)+
  xlab("Days from the event time")+
  ggtitle("TCHz for all events vs days from the event-time")

for (i in c(1:36))
{
  gg_30x30_MCDZ<-gg_30x30_MCDZ+
    geom_line(data = evs_30x30[[i]], aes(Date_obs,TCHz),colour=ev_col[i],size=.25)
  # assign each one data from each event
}
gg_30x30_MCDZ<-gg_30x30_MCDZ+scale_colour_manual(name="Event #",values=ev_col,labels=ev_label)
plot(gg_30x30_MCDZ)

stop()

######## 90x90 pixel Boxes Data #########

#search & read 36 events 90x90_pix boxes computed data
dat_90x90<-list.files(pattern = "^swp45_(.*)csv$", full.names = TRUE)
evs_90x90<-map(dat_90x90,read_csv)

#converting original date format in every tibble
for (i in seq_along(evs_90x90))
{
  evs_90x90[[i]]$Date_obs<-as_datetime(evs_90x90[[i]]$Date_obs)
}

#create a lot of gg objects
gg_90x90<-vector(mode="list", length = length(evs_90x90))
for (i in seq_along(gg_90x90))
{
  gg_90x90[[i]]<-ggplot(evs_90x90[[i]]) # assign each one, data from each event
}
