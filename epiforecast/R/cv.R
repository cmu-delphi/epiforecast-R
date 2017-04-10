## author_header begin
## Copyright (C) 2016 Sangwon Hyun
##
## This file is part of epiforecast.  Algorithms included in epiforecast were developed by Logan C. Brooks, David C. Farrow, Sangwon Hyun, Shannon Gallagher, Ryan J. Tibshirani, Roni Rosenfeld, and Rob Tibshirani (Stanford University), members of the Delphi group at Carnegie Mellon University.
##
## Research reported in this publication was supported by the National Institute Of General Medical Sciences of the National Institutes of Health under Award Number U54 GM088491. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health. This material is based upon work supported by the National Science Foundation Graduate Research Fellowship Program under Grant No. DGE-1252522. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the authors and do not necessarily reflect the views of the National Science Foundation. David C. Farrow was a predoctoral trainee supported by NIH T32 training grant T32 EB009403 as part of the HHMI-NIBIB Interfaces Initiative. Ryan J. Tibshirani was supported by NSF grant DMS-1309174.
## author_header end
## license_header begin
## epiforecast is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, version 2 of the License.
##
## epiforecast is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with epiforecast.  If not, see <http://www.gnu.org/licenses/>.
## license_header end

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
