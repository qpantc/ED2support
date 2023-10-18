#------------------------------------------------------------------------------------------#
# Calculate absolute_humidity with  temperature  and relative_humidity                     #
#------------------------------------------------------------------------------------------#

absolute_humidity <- function(atmospheric_pressure_kpa = 101.3, temperature_celsius = 25,relative_humidity = 50){
  # Atmospheric pressure in kPa

  # Convert temperature to Kelvin
  temperature_kelvin <- temperature_celsius + 273.15

  # Calculate the saturation vapor pressure (es) in kPa
  es <- 0.611 * exp((17.27 * temperature_celsius) / (temperature_celsius + 237.3))

  # Calculate the actual vapor pressure (ea) in kPa
  ea <- (relative_humidity / 100) * es

  # Calculate the specific humidity (water content per kilogram of dry air) in Kg H2O/Kg dry air
  specific_humidity <- (0.622 * ea) / (atmospheric_pressure_kpa - (0.378 * ea))

  return(specific_humidity)
}



#------------------------------------------------------------------------------------------#
# Calculate rain in kg/m2/s by unit mm per half hour                                       #
#------------------------------------------------------------------------------------------#
rainfall_kg_m2_s1 <- function(rain_mm){
  rainfall <- rain_mm / 1800
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
  fluxnet_data$atm.shv <- absolute_humidity(fluxnet_data$PA_F,fluxnet_data$TA_F,fluxnet_data$RH) #unit: "kgH2O kgAir⁻¹"
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



