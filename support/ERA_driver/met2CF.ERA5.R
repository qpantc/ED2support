###############################################################################
# met2cf.ERA5
met2CF.ERA5<- function(lat,
                       long,
                       start_date,
                       end_date,
                       sitename,
                       outfolder,
                       out.xts,
                       overwrite = FALSE,
                       verbose = TRUE) {

  years <- seq(lubridate::year(start_date),
               lubridate::year(end_date),
               1)

  ensemblesN <- seq(1, 1)

  start_date <- paste0(lubridate::year(start_date),"-01-01")  %>% as.Date()
  end_date <- paste0(lubridate::year(end_date),"-12-31") %>% as.Date()
  # adding RH and converting rain

  out.new <- ensemblesN %>%
    purrr::map(function(ensi) {
      tryCatch({

        ens <- out.xts[[ensi]]
        # Solar radation conversions
        #https://confluence.ecmwf.int/pages/viewpage.action?pageId=104241513
        #For ERA5 daily ensemble data, the accumulation period is 3 hours. Hence to convert to W/m2:

        # Reanalysis have hourly time-step!

        ens[, "ssrd"] <- ens[, "ssrd"] / (1*3600)
        ens[, "strd"] <- ens[, "strd"] / (1*3600)
        #precipitation it's originaly in meters. Meters times the density will give us the kg/m2
        ens[, "tp"] <-
          ens[, "tp"] * 1000 / 1 # divided by 3 because we have 1 hours data --> mm/h
        ens[, "tp"] <-
          udunits2::ud.convert(ens[, "tp"], "kg m-2 hr-1", "kg m-2 s-1")  #There are 21600 seconds in 6 hours??
        #RH
        #Adopted from weathermetrics/R/moisture_conversions.R
        t <-
          udunits2::ud.convert(ens[, "t2m"] %>% as.numeric(), "K", "degC")
        dewpoint  <-
          udunits2::ud.convert(ens[, "d2m"] %>% as.numeric(), "K", "degC")
        beta <- (112 - (0.1 * t) + dewpoint) / (112 + (0.9 * t))
        relative.humidity <- beta ^ 8
        #specific humidity
        specific_humidity <-
          PEcAn.data.atmosphere::rh2qair(relative.humidity,
                                         ens[, "t2m"] %>% as.numeric(),
                                         ens[, "sp"] %>% as.numeric()) # Pressure in Pa
      },
      error = function(e) {
        PEcAn.logger::logger.severe("Something went wrong during the unit conversion in met2cf ERA5.",
                                    conditionMessage(e))
      })


      #adding humidity
      xts::merge.xts(ens[, -c(3)], (specific_humidity)) %>%
        `colnames<-`(
          c(
            "air_temperature",
            "air_pressure",
            "precipitation_flux",
            "eastward_wind",
            "northward_wind",
            "surface_downwelling_shortwave_flux_in_air",
            "surface_downwelling_longwave_flux_in_air",
            "specific_humidity"
          )
        )

    })


  #These are the cf standard names
  cf_var_names = colnames(out.new[[1]])
  cf_var_units = c("K", "Pa", "kg m-2 s-1", "m s-1", "m s-1", "W m-2", "W m-2", "1")  #Negative numbers indicate negative exponents


  results_list <-  ensemblesN %>%
    purrr::map(function(i) {

      start_date <- min(zoo::index(out.new[[i]]))
      end_date <- max(zoo::index(out.new[[i]]))
      # Create a data frame with information about the file.  This data frame's format is an internal PEcAn standard, and is stored in the BETY database to
      # locate the data file.
      results <- data.frame(
        file = "",
        #Path to the file (added in loop below).
        host = PEcAn.remote::fqdn(),
        mimetype = "application/x-netcdf",
        formatname = "CF Meteorology",
        startdate = paste0(format(
          start_date , "%Y-%m-%dT%H:%M:00 %z"
        )),
        enddate = paste0(format(
          end_date , "%Y-%m-%dT%H:%M:00 %z"
        )),
        dbfile.name = paste0("ERA5.", i),
        stringsAsFactors = FALSE
      )

      # i is the ensemble number
      #Generating a unique identifier string that characterizes a particular data set.
      identifier <- paste("ERA5", sitename, i, sep = "_")

      identifier.file <- paste("ERA5",
                               i,
                               lubridate::year(start_date),
                               sep = ".")

      ensemble_folder <- file.path(outfolder, identifier)

      #Each file will go in its own folder.
      if (!dir.exists(ensemble_folder)) {
        dir.create(ensemble_folder,
                   recursive = TRUE,
                   showWarnings = FALSE)
      }

      flname <-file.path(ensemble_folder, paste(identifier.file, "nc", sep = "."))

      #Each ensemble member gets its own unique data frame, which is stored in results_list
      results$file <- flname

      years %>%
        purrr::map(function(year) {
          #
          identifier.file <- paste("ERA5",
                                   i,
                                   year,
                                   sep = ".")

          flname <-file.path(ensemble_folder, paste(identifier.file, "nc", sep = "."))
          # Spliting it for this year
          data.for.this.year.ens <- out.new[[i]]
          data.for.this.year.ens <- data.for.this.year.ens[year %>% as.character]


          #Each ensemble gets its own file.
          time_dim = ncdf4::ncdim_def(
            name = "time",
            paste(units = "hours since", format(start_date, "%Y-%m-%dT%H:%M")),
            seq(0, (length(zoo::index(
              data.for.this.year.ens
            )) * 3) - 1 , length.out = length(zoo::index(data.for.this.year.ens))),
            create_dimvar = TRUE
          )
          lat_dim = ncdf4::ncdim_def("latitude", "degree_north", lat, create_dimvar = TRUE)
          lon_dim = ncdf4::ncdim_def("longitude", "degree_east", long, create_dimvar = TRUE)

          #create a list of all ens
          nc_var_list <- purrr::map2(cf_var_names,
                                     cf_var_units,
                                     ~ ncdf4::ncvar_def(.x, .y, list(time_dim, lat_dim, lon_dim), missval = NA_real_))

          #results$dbfile.name <- flname


          if (!file.exists(flname) || overwrite) {
            tryCatch({
              nc_flptr <- ncdf4::nc_create(flname, nc_var_list, verbose = verbose)

              #For each variable associated with that ensemble
              for (j in seq_along(cf_var_names)) {
                # "j" is the variable number.  "i" is the ensemble number.
                ncdf4::ncvar_put(nc_flptr,
                                 nc_var_list[[j]],
                                 zoo::coredata(data.for.this.year.ens)[, nc_var_list[[j]]$name])
              }

              ncdf4::nc_close(nc_flptr)  #Write to the disk/storage
            },
            error = function(e) {
              PEcAn.logger::logger.severe("Something went wrong during the writing of the nc file.",
                                          conditionMessage(e))
            })

          } else {
            PEcAn.logger::logger.info(paste0(
              "The file ",
              flname,
              " already exists.  It was not overwritten."
            ))
          }


        })

      return(results)
    })
  #For each ensemble
  return(results_list )
}


met2model.ED2 <- function (in.path, in.prefix, outfolder, start_date, end_date,
          lst = 0, lat = NA, lon = NA, overwrite = FALSE, verbose = FALSE,time_interval=8,
          leap_year = TRUE, ...)
{
  overwrite <- as.logical(overwrite)
  start_date <- as.POSIXlt(start_date, tz = "UTC")
  end_date <- as.POSIXlt(end_date, tz = "UTC")
  met_folder <- outfolder
  met_header_file <- file.path(met_folder, "ED_MET_DRIVER_HEADER")
  results <- data.frame(file = met_header_file, host = PEcAn.remote::fqdn(),
                        mimetype = "text/plain", formatname = "ed.met_driver_header files format",
                        startdate = start_date, enddate = end_date, dbfile.name = "ED_MET_DRIVER_HEADER",
                        stringsAsFactors = FALSE)
  dir.create(met_folder, recursive = TRUE, showWarnings = FALSE)
  month <- c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL",
             "AUG", "SEP", "OCT", "NOV", "DEC")
  day2mo <- function(year, day, leap_year) {
    dm <- c(0, 32, 60, 91, 121, 152, 182, 213, 244, 274,
            305, 335, 366)
    dl <- c(0, 32, 61, 92, 122, 153, 183, 214, 245, 275,
            306, 336, 367)
    mo <- rep(NA, length(day))
    if (!leap_year) {
      mo <- findInterval(day, dm)
      return(mo)
    }
    else {
      leap <- lubridate::leap_year(year)
      mo[leap] <- findInterval(day[leap], dl)
      mo[!leap] <- findInterval(day[!leap], dm)
      return(mo)
    }
  }
  start_year <- lubridate::year(start_date)
  end_year <- lubridate::year(end_date)
  year_seq <- seq(start_year, end_year)
  day_secs <- PEcAn.utils::ud_convert(1, "day", "seconds")
  need_input_files <- file.path(in.path, paste(in.prefix,
                                               year_seq, "nc", sep = "."))
  have_input_files <- file.exists(need_input_files)
  if (!all(have_input_files)) {
    PEcAn.logger::logger.severe("Missing the following required input files: ",
                                paste(sprintf("'%s'", need_input_files[!have_input_files]),
                                      collapse = ", "))
  }
  month_seq <- seq(lubridate::floor_date(start_date, "month"),
                   lubridate::floor_date(end_date, "month"), by = "1 month")
  target_fnames <- paste0(toupper(strftime(month_seq, "%Y%b",
                                           tz = "UTC")), ".h5")
  target_out_files <- file.path(met_folder, target_fnames)
  have_target_out_files <- file.exists(target_out_files)
  if (any(have_target_out_files)) {
    if (overwrite) {
      PEcAn.logger::logger.warn("The following existing target output files will be overwritten:",
                                paste(sprintf("'%s'", target_out_files[have_target_out_files]),
                                      collapse = ", "))
    }
    else {
      have_output_byyear <- split(have_target_out_files,
                                  lubridate::year(month_seq))
      complete_years <- vapply(have_output_byyear, all,
                               logical(1))
      skip_years <- tryCatch(as.numeric(names(complete_years[complete_years])),
                             warning = function(e) PEcAn.logger::logger.severe(e))
      PEcAn.logger::logger.warn("The following output files already exist:",
                                paste(target_out_files[have_target_out_files]),
                                ". This means the following complete years will be skipped: ",
                                skip_years)
      year_seq <- setdiff(year_seq, skip_years)
    }
  }
  for (year in year_seq) {
    ncfile <- file.path(in.path, paste(in.prefix, year,
                                       "nc", sep = "."))
    nc <- ncdf4::nc_open(ncfile)
    flat <- try(ncdf4::ncvar_get(nc, "latitude"), silent = TRUE)
    if (!is.numeric(flat)) {
      flat <- nc$dim[[1]]$vals[1]
    }
    if (is.na(lat)) {
      lat <- drop(flat)
    }
    else if (lat != flat) {
      PEcAn.logger::logger.warn("Latitude does not match that of file",
                                lat, "!=", flat)
    }
    flon <- try(ncdf4::ncvar_get(nc, "longitude"), silent = TRUE)
    if (!is.numeric(flon)) {
      flon <- nc$dim[[2]]$vals[1]
    }
    if (is.na(lon)) {
      lon <- drop(flon)
    }
    else if (lon != flon) {
      PEcAn.logger::logger.warn("Longitude does not match that of file",
                                lon, "!=", flon)
    }
    tdays <- nc$dim$time$vals
    Tair <- ncdf4::ncvar_get(nc, "air_temperature")
    Qair <- ncdf4::ncvar_get(nc, "specific_humidity")
    U <- try(ncdf4::ncvar_get(nc, "eastward_wind"), silent = TRUE)
    V <- try(ncdf4::ncvar_get(nc, "northward_wind"), silent = TRUE)
    Rain <- ncdf4::ncvar_get(nc, "precipitation_flux")
    pres <- ncdf4::ncvar_get(nc, "air_pressure")
    SW <- ncdf4::ncvar_get(nc, "surface_downwelling_shortwave_flux_in_air")
    LW <- ncdf4::ncvar_get(nc, "surface_downwelling_longwave_flux_in_air")
    CO2 <- try(ncdf4::ncvar_get(nc, "mole_fraction_of_carbon_dioxide_in_air"),
               silent = TRUE)
    use_UV <- is.numeric(U) & is.numeric(V)
    if (!use_UV) {
      U <- try(ncdf4::ncvar_get(nc, "wind_speed"), silent = TRUE)
      if (is.numeric(U)) {
        PEcAn.logger::logger.info("eastward_wind and northward_wind are absent, using wind_speed to approximate eastward_wind")
        V <- rep(0, length(U))
      }
      else {
        PEcAn.logger::logger.severe("No eastward_wind and northward_wind or wind_speed in the met data")
      }
    }
    useCO2 <- is.numeric(CO2)
    sec <- PEcAn.utils::ud_convert(tdays, unlist(strsplit(nc$dim$time$units,
                                                          " "))[1], "seconds")
    ncdf4::nc_close(nc)
    dt <- drop(unique(diff(sec)))
    if (length(dt) > 1) {
      dt_old <- dt
      dt <- drop(round(mean(diff(sec))))
      PEcAn.logger::logger.warn(paste0("Time step (`dt`) is not uniform! Identified ",
                                       length(dt_old), " unique time steps. ", "`head(dt)` (in seconds): ",
                                       paste(utils::head(dt_old), collapse = ", "),
                                       " Using the rounded mean difference as the time step: ",
                                       dt))
    }
    toff <- -as.numeric(lst) * 3600/dt
    slen <- seq_along(sec)
    Tair <- c(rep(Tair[1], toff), Tair)[slen]
    Qair <- c(rep(Qair[1], toff), Qair)[slen]
    U <- c(rep(U[1], toff), U)[slen]
    V <- c(rep(V[1], toff), V)[slen]
    Rain <- c(rep(Rain[1], toff), Rain)[slen]
    pres <- c(rep(pres[1], toff), pres)[slen]
    SW <- c(rep(SW[1], toff), SW)[slen]
    LW <- c(rep(LW[1], toff), LW)[slen]
    if (useCO2) {
      CO2 <- c(rep(CO2[1], toff), CO2)[slen]
    }
    doy <- floor(tdays/time_interval) + 1
    invalid_doy <- doy < 1 | doy > PEcAn.utils::days_in_year(year,
                                                             leap_year)
    if (any(invalid_doy)) {
      PEcAn.logger::logger.severe(paste0("Identified at least one invalid day-of-year (`doy`). ",
                                         "PEcAn met standard uses days since start of year as its time unit, ",
                                         "so this suggests a problem with the input met file. ",
                                         "Invalid values are: ", paste(doy[invalid_doy],
                                                                       collapse = ", "), ". ", "Source file is: ",
                                         normalizePath(ncfile)))
    }
    hr <- ((tdays/time_interval)%%1) * 24
    cosz <- PEcAn.data.atmosphere::cos_solar_zenith_angle(doy,
                                                          lat, lon, dt, hr)
    rpot <- 1366 * cosz
    rpot <- rpot[seq_along(tdays)]
    SW[rpot < SW] <- rpot[rpot < SW]
    frac <- SW/rpot
    frac[frac > 0.9] <- 0.9
    frac[frac < 0] <- 0
    frac[is.na(frac)] <- 0
    frac[is.nan(frac)] <- 0
    SWd <- SW * (1 - frac)
    n <- length(Tair)
    nbdsfA <- (SW - SWd) * 0.57
    nddsfA <- SWd * 0.48
    vbdsfA <- (SW - SWd) * 0.43
    vddsfA <- SWd * 0.52
    prateA <- Rain
    dlwrfA <- LW
    presA <- pres
    hgtA <- rep(50, n)
    ugrdA <- U
    vgrdA <- V
    shA <- Qair
    tmpA <- Tair
    if (useCO2) {
      co2A <- CO2 * 1e+06
    }
    mo <- day2mo(year, doy, leap_year)
    for (m in unique(mo)) {
      selm <- which(mo == m)
      mout <- file.path(met_folder, paste0(year, month[m],
                                           ".h5"))
      if (file.exists(mout)) {
        if (overwrite) {
          file.remove(mout)
          ed_met_h5 <- hdf5r::H5File$new(mout)
        }
        else {
          PEcAn.logger::logger.warn("The file already exists! Moving to next month!")
          next
        }
      }
      else {
        ed_met_h5 <- hdf5r::H5File$new(mout)
      }
      dims <- c(length(selm), 1, 1)
      nbdsf <- array(nbdsfA[selm], dim = dims)
      nddsf <- array(nddsfA[selm], dim = dims)
      vbdsf <- array(vbdsfA[selm], dim = dims)
      vddsf <- array(vddsfA[selm], dim = dims)
      prate <- array(prateA[selm], dim = dims)
      dlwrf <- array(dlwrfA[selm], dim = dims)
      pres <- array(presA[selm], dim = dims)
      hgt <- array(hgtA[selm], dim = dims)
      ugrd <- array(ugrdA[selm], dim = dims)
      vgrd <- array(vgrdA[selm], dim = dims)
      sh <- array(shA[selm], dim = dims)
      tmp <- array(tmpA[selm], dim = dims)
      if (useCO2) {
        co2 <- array(co2A[selm], dim = dims)
      }
      ed_met_h5[["nbdsf"]] <- nbdsf
      ed_met_h5[["nddsf"]] <- nddsf
      ed_met_h5[["vbdsf"]] <- vbdsf
      ed_met_h5[["vddsf"]] <- vddsf
      ed_met_h5[["prate"]] <- prate
      ed_met_h5[["dlwrf"]] <- dlwrf
      ed_met_h5[["pres"]] <- pres
      ed_met_h5[["hgt"]] <- hgt
      ed_met_h5[["ugrd"]] <- ugrd
      ed_met_h5[["vgrd"]] <- vgrd
      ed_met_h5[["sh"]] <- sh
      ed_met_h5[["tmp"]] <- tmp
      if (useCO2) {
        ed_met_h5[["co2"]] <- co2
      }
      ed_met_h5$close_all()
    }
    metvar <- c("nbdsf", "nddsf", "vbdsf", "vddsf", "prate",
                "dlwrf", "pres", "hgt", "ugrd", "vgrd", "sh", "tmp",
                "co2")
    metvar_table <- data.frame(variable = metvar, update_frequency = dt,
                               flag = 1)
    if (!useCO2) {
      metvar_table_vars <- metvar_table[metvar_table$variable !=
                                          "co2", ]
    }
    else {
      metvar_table_vars <- metvar_table
    }
    ed_metheader <- list(list(path_prefix = met_folder,
                              nlon = 1, nlat = 1, dx = 1, dy = 1, xmin = lon,
                              ymin = lat, variables = metvar_table_vars))
    PEcAn.ED2::check_ed_metheader(ed_metheader)
    PEcAn.ED2::write_ed_metheader(ed_metheader, met_header_file, header_line = shQuote("Made_by_PEcAn_met2model.ED2"))
  }
  PEcAn.logger::logger.info("Done with met2model.ED2")
  return(invisible(results))
}
