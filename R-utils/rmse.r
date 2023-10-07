#==========================================================================================#
#==========================================================================================#
#     Generic RMSE.                                                                        #
#------------------------------------------------------------------------------------------#
rmse.gen <<- function(yobs,yhat,np=1,na.rm=TRUE,unitless=TRUE){

   #------ Discard data in case na.rm = TRUE. ---------------------------------------------#
   if (na.rm){
      keep = is.finite(yobs) & is.finite(yhat)
      yobs = yobs[keep]
      yhat = yhat[keep]
   }#end if (na.rm)
   yres = yobs - yhat
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     In case length is less than or equal to np, return NA.                            #
   #---------------------------------------------------------------------------------------#
   ny = length(yobs)
   if (ny <= (np+1)){
      ans = NA
      return(ans)
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Get the number of degrees of freedom, sum of squares, and root mean square       #
   # error.                                                                                #
   #---------------------------------------------------------------------------------------#
   df.tot   = ny - 1
   df.exp   = ny - np
   mse      = sum(yres^2) / df.exp
   if (unitless){
      ans   = sqrt(mse) / sd(yobs)
   }else{
      ans   = sqrt(mse)
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end rmse.gen
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Generic bias.                                                                        #
#------------------------------------------------------------------------------------------#
bias.gen <<- function(yobs,yhat,na.rm=TRUE,unitless=TRUE){

   #------ Discard data in case na.rm = TRUE. ---------------------------------------------#
   if (na.rm){
      keep = is.finite(yobs) & is.finite(yhat)
      yobs = yobs[keep]
      yhat = yhat[keep]
   }#end if (na.rm)
   yres = yobs - yhat
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     In case length is less than or equal to np, return NA.                            #
   #---------------------------------------------------------------------------------------#
   ny = length(yobs)
   if (ny <= 0){
      ans = NA
      return(ans)
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Get the mean bias.                                                               #
   #---------------------------------------------------------------------------------------#
   if (unitless){
      ans   = - mean(yres) / sd(yobs)
   }else{
      ans   = - mean(yres)
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end bias.gen
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Generic sigma of residuals.                                                          #
#------------------------------------------------------------------------------------------#
sigres.gen <<- function(yobs,yhat,np=1,na.rm=TRUE,unitless=TRUE){

   #------ Discard data in case na.rm = TRUE. ---------------------------------------------#
   if (na.rm){
      keep = is.finite(yobs) & is.finite(yhat)
      yobs = yobs[keep]
      yhat = yhat[keep]
   }#end if (na.rm)
   yres = yobs - yhat
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     In case length is less than or equal to np, return NA.                            #
   #---------------------------------------------------------------------------------------#
   ny = length(yobs)
   if (ny <= 0){
      ans = NA
      return(ans)
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Variance of the residuals.                                                       #
   #---------------------------------------------------------------------------------------#
   df.tot   = ny - 1
   df.exp   = ny - np
   varres   = sum((yres-mean(yres))^2) / df.exp
   varobs   = sum((yobs-mean(yobs))^2) / df.tot
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Get the mean sigma of residuals.                                                 #
   #---------------------------------------------------------------------------------------#
   if (unitless){
      ans   = sqrt(varres) / sqrt(varobs)
   }else{
      ans   = sqrt(varres)
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end sigres.gen
#==========================================================================================#
#==========================================================================================#
