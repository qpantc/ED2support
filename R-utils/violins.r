#==========================================================================================#
#==========================================================================================#
#    This function is very similar to vioplot, except that it allows different violin      #
# widths to be scaled by a single KDE, and uses function density instead of sm.density.    #
#  Also, all elements must be provided as a list.                                          #
#------------------------------------------------------------------------------------------#
violins <<- function( x
                    , at         = seq_along(x)
                    , xnames     = if(is.null(names(x))){seq_along(x)}else{names(x)}
                    , box.range  = 1.5
                    , ylim       = NULL
                    , horizontal = FALSE
                    , col        = "grey66"
                    , border     = "black"
                    , lty        = 1
                    , lwd        = 1
                    , rectCol    = "black"
                    , colMed     = "white"
                    , bgMed      = colMed
                    , pchMed     = 19
                    , add        = FALSE
                    , wex        = 1
                    , drawRect   = TRUE
                    , plog       = FALSE
                    , from
                    , to
                    , pwidth     = FALSE
                    , n          = 512
                    , na.rm      = FALSE
                    , ...
                    ){

   #----- Initialise several vectors that will be used to generate the violins. -----------#
   nx        = length(x)
   x.upper   = numeric(length=nx)
   x.lower   = numeric(length=nx)
   x.q250    = numeric(length=nx)
   x.q500    = numeric(length=nx)
   x.q750    = numeric(length=nx)
   x.width   = numeric(length=nx)
   base      = replicate(n=nx,list())
   height    = replicate(n=nx,list())
   baserange = NULL
   #---------------------------------------------------------------------------------------#



   #----- Make sure the violin settings can vary. -----------------------------------------#
   seqx    = sequence(nx)
   col     = rep(x=col    ,times=nx)[seqx]
   border  = rep(x=border ,times=nx)[seqx]
   lty     = rep(x=lty    ,times=nx)[seqx]
   lwd     = rep(x=lwd    ,times=nx)[seqx]
   rectCol = rep(x=rectCol,times=nx)[seqx]
   colMed  = rep(x=colMed ,times=nx)[seqx]
   pchMed  = rep(x=pchMed ,times=nx)[seqx]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Define the common axis.                                                           #
   #---------------------------------------------------------------------------------------#
   if (missing(from) | missing(to)){
      x.range = range(unlist(x),finite=TRUE)
      if (missing(from)) from = x.range[1]
      if (missing(to)  ) to   = x.range[2]
      dbase = mean(diff(seq(from=from,to=to,length.out=n)))
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Go through each list element, and estimate the density function.                  #
   #---------------------------------------------------------------------------------------#
   for (ix in seqx){
       #----- Load data, remove missing values in case na.rm is TRUE. ---------------------#
       xnow = x[[ix]]
       if (na.rm) xnow = xnow[! is.na(xnow)]
       #-----------------------------------------------------------------------------------#


       #----- Find the boxplot limits. ----------------------------------------------------#
       x.min       = min(xnow)
       x.max       = max(xnow)
       x.q250[ix]  = quantile(xnow, 0.25)
       x.q750[ix]  = quantile(xnow, 0.75)
       x.q500[ix]  = median(xnow)
       x.iqd       = x.q750[ix] - x.q250[ix]
       x.upper[ix] = min(x.q750[ix] + box.range * x.iqd, x.max)
       x.lower[ix] = max(x.q250[ix] - box.range * x.iqd, x.min)
       #-----------------------------------------------------------------------------------#


       #----- Estimate kernel density, keep only values within range. ---------------------#
       densout      = density.safe(xnow,from=from,to=to,n=n,...)
       keep         = densout$x %wr% c(x.min-0.5*dbase,x.max+0.5*dbase)
       base  [[ix]] = densout$x[keep]
       height[[ix]] = densout$y[keep]
       if (any(is.finite(height[[ix]]))){
          x.width[ix] = max(height[[ix]],na.rm=TRUE)
       }else{
          x.width[ix] = NA_real_
       }#end if
       #-----------------------------------------------------------------------------------#

       #----- Update candidate limit. -----------------------------------------------------#
       baserange   = range(c(baserange,base[[ix]]),finite=TRUE)
       #-----------------------------------------------------------------------------------#
   }#end for (ix in seqx)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Define x and y limits in case they aren't provided.                              #
   #---------------------------------------------------------------------------------------#
   if (! add){
      if (nx == 1){
         xlim = c(at-0.5,at+0.5)
      }else{
         xlim = range(at) + min(diff(at))/2 * c(-1, 1)
      }#end if
      if (is.null(ylim)) ylim = baserange
   }#end if (! add)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Standardise the scaling factor.                                                   #
   #---------------------------------------------------------------------------------------#
   for (i in sequence(nx)){
      if (pwidth){
         x.scale = 0.4 / max(x.width,na.rm=TRUE)
      }else if (is.finite(x.width[i])){
         x.scale = 0.4 / x.width[i]
      }else{
         x.scale = 0.4
      }#end if (pwidth)

      height[[i]] = height[[i]] * x.scale
   }#end for (i in sequence(nx))
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #     Define the box width.                                                             #
   #---------------------------------------------------------------------------------------#
   boxwidth = 0.05 * wex
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Start device.                                                                     #
   #---------------------------------------------------------------------------------------#
   if (! add){
      plot.new()
      if (horizontal){
         plot.window(xlim=ylim,ylim=xlim,log=ifelse(plog,"x",""))
         axis(1)
         axis(2,at=at,labels=xnames)
      }else{
         plot.window(xlim=xlim,ylim=ylim,log=ifelse(plog,"y",""))
         axis(1)
         axis(2,at=at,labels=xnames)
      }#end if (horizontal)
   }#end if (! add)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Go through the list elements.                                                     #
   #---------------------------------------------------------------------------------------#
   for (ix in seqx){

      #----- Polygon coordinates. ---------------------------------------------------------#
      poly.height = c(at[ix] - height[[ix]],rev(at[ix] + height[[ix]]))
      poly.base   = c(base[[ix]],rev(base[[ix]]))
      #------------------------------------------------------------------------------------#


      if (horizontal){
         #----- Plot violin. --------------------------------------------------------------#
         epolygon( x      = poly.base
                 , y      = poly.height
                 , col    = col   [ix]
                 , border = border[ix]
                 , lty    = lty   [ix]
                 , lwd    = lwd   [ix]
                 )#end epolygon
         #---------------------------------------------------------------------------------#


         #----- Plot box-and-whisker plot. ------------------------------------------------#
         if (drawRect){
             lines (x=c(x.lower[ix], x.upper[ix]),y=at[c(ix,ix)],lwd=lwd[ix],lty=lty[ix])
             rect  ( xleft   = x.q250[ix]
                   , ybottom = at[ix] - boxwidth/2
                   , xright  = x.q750[ix]
                   , ytop    = at[ix] + boxwidth/2
                   , col     = rectCol
                   )#end rect
             points(x=x.q500[ix],y=at[ix],pch=pchMed[ix],col=colMed[ix],bg=bgMed[ix])
         }#end if (drawRect)
         #---------------------------------------------------------------------------------#
      }else{
         #----- Plot violin. --------------------------------------------------------------#
         epolygon( x      = poly.height
                 , y      = poly.base
                 , col    = col   [ix]
                 , border = border[ix]
                 , lty    = lty   [ix]
                 , lwd    = lwd   [ix]
                 )#end epolygon
         #---------------------------------------------------------------------------------#


         #----- Plot box-and-whisker plot. ------------------------------------------------#
         if (drawRect){
             lines (x=at[c(ix,ix)],y=c(x.lower[ix],x.upper[ix]),lwd=lwd[ix],lty=lty[ix])
             rect  ( xleft   = at[ix] - boxwidth/2
                   , ybottom = x.q250[ix]
                   , xright  = at[ix] + boxwidth/2
                   , ytop    = x.q750[ix]
                   , col     = rectCol
                   )#end rect
             points(x=at[ix],y=x.q500[ix],pch=pchMed[ix],col=colMed[ix],bg=bgMed[ix])
         }#end if (drawRect)
         #---------------------------------------------------------------------------------#
      }#end if (horizontal)
      #------------------------------------------------------------------------------------#
   }#end for (ix in seqx)
   #---------------------------------------------------------------------------------------#


   #----- Make invisible output. ----------------------------------------------------------#
   ans = list( upper  = x.upper
             , lower  = x.lower
             , median = x.q500
             , qlwr   = x.q250
             , qupr   = x.q750
             , width  = x.width
             )#end list
   invisible(ans)
   #---------------------------------------------------------------------------------------#
}#end function violins
#==========================================================================================#
#==========================================================================================#
