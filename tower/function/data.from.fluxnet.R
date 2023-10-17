#------------------------------------------------------------------------------------------#
# Calculate absolute_humidity with  temperature  and relative_humidity                     #
#------------------------------------------------------------------------------------------#

absolute_humidity <- function(temperature = 25,relative_humidity = 50){
  R = 287.05 # Specific gas constant for dry air (J/(kg*K))
  eps = 0.622 # Ratio of the molecular weight of water vapor to dry air
  (eps * (relative_humidity / 100) * 6.112 * exp((17.67 * temperature) / (temperature + 243.5))) /
    (R * (temperature + 273.15))
}


#------------------------------------------------------------------------------------------#
# Calculate rain in kg/m2/s by unit mm per half hour                                       #
#------------------------------------------------------------------------------------------#
rainfall_kg_m2_s1 <- function(rain){
  rainfall <- (rain * 0.000001) / 1800
  return(rainfall)
}

#------------------------------------------------------------------------------------------#
# Calculate absolute_humidity with temperature  and relative_humidity                     #
#------------------------------------------------------------------------------------------#
#     need variables:                                                                   #
# year — Year (time in UTC)
# month — Month (time in UTC)
# day — Day (time in UTC)
# hour — Hour (time in UTC)
# min — Minute (time in UTC)
# sec — Second (time in UTC)
# atm.tmp — Air temperature (K)
# atm.prss — Atmospheric pressure (Pa)
# atm.shv — Specific humidity (kgH2O kgAir-1)
# atm.vels — Wind speed (m s-1)
# atm.vdir — Wind direction (°). This is not needed by ED-2. add a dummy column with zeroes.
# rain — Precipitation rate (kg m-2 s-1)
# rlong.in — Incoming longwave radiation (W m-2)
# rshort.in — Incoming shortwave radiation (W m-2)
#---------------------------------------------------------------------------------------#
fluxnet_vesion_data <- function(file.in){
  fluxnet_data = read.csv(file.in)
  colnames(fluxnet_data)
  fluxnet_data$time = ymd_hm(fluxnet_data$TIMESTAMP_START)
  fluxnet_data$year <- year(fluxnet_data$time)
  fluxnet_data$month <- month(fluxnet_data$time)
  fluxnet_data$day <- day(fluxnet_data$time)
  fluxnet_data$hour <- hour(fluxnet_data$time)
  fluxnet_data$min <- minute(fluxnet_data$time)
  fluxnet_data$sec <- second(fluxnet_data$time)

  fluxnet_data$atm.tmp <- fluxnet_data$TA_F + 273.15 # unit: K
  fluxnet_data$atm.prss <- fluxnet_data$PA_F*1000 #unit: Pa
  fluxnet_data$atm.shv <- absolute_humidity(fluxnet_data$TA_F,fluxnet_data$RH) #unit: "kgH2O kgAir⁻¹"
  fluxnet_data$atm.vels <- fluxnet_data$WS_F
  fluxnet_data$atm.vdir <- fluxnet_data$WD
  fluxnet_data$rain <- rainfall_kg_m2_s1(fluxnet_data$P_F) # unit: kg/m2*s1
  fluxnet_data$rlong.in <- fluxnet_data$LW_IN_F
  fluxnet_data$rshort.in <- fluxnet_data$SW_IN_F
  datum <- fluxnet_data[c('year','month','day','hour','min','sec',
                   'atm.tmp','atm.prss','atm.shv','atm.vels',
                   'atm.vdir','rain','rlong.in','rshort.in')]

  return(datum)
}



