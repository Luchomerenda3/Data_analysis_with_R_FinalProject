###################
#Author:          # 
#LucianoA.Merenda #
#2019             #
###################

swp_mp<-function(mordyn_ts,swp_ts,del_t,erup_time)
     {
      #Esta Función recibe un time-serie de alguna caracteristica morfologico y/o dinamica de la cme "mordyn_ts",
      #un arreglo de 36 time-series de algún parámetro fotosférico "swp_ts",
      #un valor de tiempo en formato string (e.g. "1 hour", ver ?ts_lag) correspondiente al tiempo a atrasar
      #o adelantar cada time-serie("swp_ts") "del_t",
      #y un vector de largo 36 con los tiempos de cada evento eyectivo "erup_time"
      #Luego se usa en conjunto con "sw_md" aunque puede utilizarse sola.
      
      nsw<-length(swp_ts) #take number of events! to generalize the function.
      
      #We create the swp_time_series to correlate with the mordyn_ts
      swp<-vector(mode="numeric",nsw)
      #assign its values
      for_swp<-swp_ts #local variable created to maintain the original swp_ts safe while we work with it
      for (i in c(1:nsw))          
        {
         print(str_c("Working on event ",as.character(i)," of ",as.character(nsw)))
         len=nrow(for_swp[[i]])
         
         #First we lag or lead the input space_weather_parameter time-serie by delta_t ("del_t") 
         for_swp[[i]]<-slice(ts_lag(for_swp[[i]], by=del_t),1:len)
         
         #Then we have to generate the time series for the value of the swp at the eruption_time
         swp[i]<-filter(for_swp[[i]],abs(as.numeric(difftime(time,erup_time[i],units="secs")))<361) %>%
                          summarise(avg=mean(value)) %>%
                          as.numeric() 
         
         #IMPORTANT_NOTE!: there's 720 seconds between each observation in the space_weather_quantities (e.g. usfz),
         #When we are comparing which corresponds to the eruption time, it may fall between 2 observations! 
         #so the difference between the erup_time and the observation is at least 360s (right in the middle)
         #so we take 361 and compute the mean value if the erup falls right in the middle, otherwise it'll fall
         # in either one of the 2 neareast observation! Hence the 361 comparison value.
         ### los NAN aparecen acá! fixing->>>
         #check if swp[i] is NaN, this may happen whether the cadence isnt 720s, or the eruption-time doesnt fall
         #in the leaded/lagged time-series.
         if (is.nan(swp[i])) 
           {
            #First we check if cadence is fucked up so we change the 360s in intervals of 360s until we found the nearest.
            k=1
            aux=361
            while(k<37 && is.nan(swp[i])) 
               {
                print("Correcting NaN since no value has been found.")
                swp[i]<-filter(for_swp[[i]],abs(as.numeric(difftime(time,erup_time[i],units="secs")))<aux) %>%
                summarise(avg=mean(value)) %>%
                as.numeric()
                aux<-aux+360
               }
            if (is.nan(swp[i]))
               {
                print("Correcting NaN since last correction wasnt able to found one,
                      this means the value of erup_time is no longer in the time-series")
                #If this checks true that means the eruption times is no longer in the time series  
                #and the closest value for swp[i] will be the first value for that time-series or we could assign 0
                #since we are checking if previous values may correlate with the morphologic parameter, 
                #we will assign the closest value which is the first of for_swp[[i]]
                swp[i]<-as.numeric(for_swp[[i]]$value[1])
                print(str_c("WARNING!!! replacing NaN value with closest value for that eruption time,
                            Difference in time between values is ",
                            as.character(as.numeric(difftime(for_swp[[i]]$time[i],erup_time[i],units="secs"))/3600),
                            " hours"))
                #swp[i]<-0.
               }
            }
        }
      
      
      #Final tibble to compute the linear regression 
      lr_ts<-add_column(mordyn_ts,swp)
      
      attr(lr_ts,"delta_t_used")<-del_t #Here we add the delta_t used to compute the time-series, 
                                        #to track it down if needed
      
      #with cor() we can compute correlation coefficient.
      #betwenn swp our "x" and the morpho/dynamic our "y"

      cor_fac<-cor(lr_ts[[3]],as.numeric(lr_ts[[2]]))
      
      #linear regression between the swp and mordyn_ts
      #same here swp our "x" and the morpho/dynamic our "y"
      lr_swpmor<-lm(as.numeric(lr_ts[[2]])~lr_ts[[3]])
      
      #Define the final object to return containing the final tibble, correlation factor, and linear regression.
      flm_swpmor<-list(lr_ts,cor_fac,lr_swpmor)
      
      #Return -> linear regression coefficient? entire linear regression object so I can plot it? both? yeahhh 
      return(flm_swpmor)
     }

sw_md<-function(mordyn_ts,swp_ts,delta_t,erup_time,delta_tnum)
       {
        #Esta Función recibe un time-serie de alguna caracteristica morfologico y/o dinamica de la cme "mordyn_ts",
        #un arreglo de 36 time-series de algún parámetro fotosférico "swp_ts",
        #un array de valor de tiempo en formato string (e.g. "1 hour", ver ?ts_lag) correspondiente al tiempo a atrasar
        #o adelantar cada time-serie("swp_ts") "delta_t",
        #y un vector de largo 36 con los tiempos de cada evento eyectivo "erup_time"
        #Tiene como simple objetivo extender la función "swp_mp" a una cantidad length(delta_t) de adelantos.
        
        n_delt<-length(delta_t) # number of delta_t
        
        #Final object to return!
        #IMPORTANT NOTE : This will be a list of n_delt elements,each one containing
        #a tibble, a correlation factor for variables in the tibble, and a linearmodel for those variables.
        swp_vs_mp<-vector(mode = "list",n_delt)
        
        for (i in c(1:n_delt))
           {
            print(str_c("###########WORKING ON DELTA_T ",as.character(i)," OF ",as.character(n_delt),"#################"))
            swp_vs_mp[[i]]<-swp_mp(mordyn_ts,swp_ts,delta_t[i],erup_time)
           }
        
        #We will now change the structure of the returned object so its easier to manipulate
        tibbles<-vector(mode = "list",n_delt)
        cor_facs<-vector(mode = "numeric",n_delt)
        lr_models<-vector(mode = "list",n_delt)
        
        for (i in c(1:n_delt)) 
           {
           tibbles[[i]]<-swp_vs_mp[[i]][[1]]
           cor_facs[[i]]<-swp_vs_mp[[i]][[2]]
           lr_models[[i]]<-swp_vs_mp[[i]][[3]]
        }
        
        tibbles<-list("swpvsmdp"=tibbles,"delta_tnum"=delta_tnum)
        cor_facs<-add_column(as_tibble(cor_facs),delta_tnum)
        lr_models<-list("lrmodels"=lr_models,"delta_tnum"=delta_tnum)
        
        fswp_vs_mp<-list("swp_vs_mdp"=tibbles,"Corr_factor"=cor_facs,"LinReg_models"=lr_models)
        
        return(fswp_vs_mp)
        
       }

snum2str<-function(numsec,unitt)
     {
      #Esta función toma de argumento tipo numeric y se asume que son segundos.
      #luego transforma este número y lo traduce a un string en unidades especificadas por units.
      #Units can be one of one of "sec", "min", "hour", "day", "week", "month", "quarter" or "year" *from ts_lag
      switch(unitt,
             sec=str_c(as.character(numsec)," sec"),
             min=str_c(as.character(numsec/60)," min"),
             hour=str_c(as.character(numsec/3600),"h"),
             day=str_c(as.character(numsec/(3600*24))," d"),
             week=str_c(as.character(numsec/(3600*24*7))," w")
             )
      
     }
