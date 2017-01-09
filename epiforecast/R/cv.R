##' Compare two cv.sim objects
cv.compare = function(cv.sim1, cv.sim2){
    print("A helpful message like this: setting 1 (wiggle, holiday effect, etc.) is better than setting 2( wiggle, etc)")
}



##' Class for cv objects; contains
##' (1) control list
##' (2) for each forecasting time and for each left-out season (fold), what are the densities in the 52 by n.grid block of the 2d plane?\
##' (3) several things pre-calculated; the prediction scores (negative log-likelihood) for 4 target forecasts.
##' The idea is that, instead of 52 by nsim curves, we can store 52 by
##' n.grid values that store the density estimates, where n.grid can be
##' hundreds, while n.sim may be 10,000's.
##' It should return /all/ quantities required for doing 
cv.sim = function(){
    print("Not written yet")
return()    
}

##' Holiday effect:
