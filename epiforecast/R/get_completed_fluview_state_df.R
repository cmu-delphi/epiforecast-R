## author_header begin
## Copyright (C) 2017 Aaron Rumack
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

get_completed_fluview_state_df = function(epigroup, first.week.of.season=31L) {
  st = stringi::stri_trans_tolower(epigroup)
  st.dat = NULL
  if (st != "fl" && st != "la") {
    st.dat = fetchEpidataDF("fluview", st, first.week.of.season = first.week.of.season,
                            cache.file.prefix=sprintf("fluview_%s_fetch", st))
    if (st == "pr") {
      ## Estimate seasons 2010-2013 using other epigroups in Region 2
      njili = get_completed_fluview_state_df("nj",first.week.of.season = first.week.of.season)$ili
      nyili = get_completed_fluview_state_df("ny",first.week.of.season = first.week.of.season)$ili
      jfkili = get_completed_fluview_state_df("jfk",first.week.of.season = first.week.of.season)$ili
      pr.len = length(st.dat$ili)
      pr.fit = lm(st.dat$ili ~ tail(njili,pr.len) + tail(nyili,pr.len) + tail(jfkili,pr.len))
      pr.est = pr.fit$coefficients[1] + pr.fit$coefficients[2]*njili + pr.fit$coefficients[3]*nyili + pr.fit$coefficients[4]*jfkili
      pr.dat = get_completed_fluview_state_df("nj") # Placeholder for values such as epiweek, issue, etc.
      vars = colnames(pr.dat)[startsWith(colnames(pr.dat),"num_")] # Flush numeric variables in data frame
      firstepiweek = min(st.dat$epiweek[!is.na(st.dat$issue)])
      mask = pr.dat$epiweek >= firstepiweek # Have data from here and later
      pr.dat[mask,vars] = st.dat[st.dat$epiweek >= firstepiweek,vars]
      pr.dat$ili[!mask] = pr.est[!mask]
      pr.dat$ili[mask] = st.dat$ili[st.dat$epiweek >= firstepiweek]
      pr.dat$wili = pr.dat$ili
      pr.dat$region = st
      st.dat = pr.dat
    } else if (st == "vi") {
      ## Estimate seasons 2010-2015 using other epigroups in Region 2
      njili = get_completed_fluview_state_df("nj",first.week.of.season = first.week.of.season)$ili
      nyili = get_completed_fluview_state_df("ny",first.week.of.season = first.week.of.season)$ili
      jfkili = get_completed_fluview_state_df("jfk",first.week.of.season = first.week.of.season)$ili
      prili = get_completed_fluview_state_df("pr",first.week.of.season = first.week.of.season)$ili
      vi.len = length(st.dat$ili)
      vi.fit = lm(st.dat$ili ~ tail(njili,vi.len) + tail(nyili,vi.len) + tail(jfkili,vi.len) + tail(prili,vi.len))
      vi.est = vi.fit$coefficients[1] + vi.fit$coefficients[2]*njili + vi.fit$coefficients[3]*nyili + vi.fit$coefficients[4]*jfkili + vi.fit$coefficients[5]*prili
      vi.dat = get_completed_fluview_state_df("nj") # Placeholder for values such as epiweek, issue, etc.
      vars = colnames(vi.dat)[startsWith(colnames(vi.dat),"num_")] # Flush numeric variables in data frame
      firstepiweek = min(st.dat$epiweek[!is.na(st.dat$issue)])
      mask = vi.dat$epiweek >= firstepiweek # Have data from here and later
      vi.dat[mask,vars] = st.dat[st.dat$epiweek >= firstepiweek,vars]
      vi.dat$ili[!mask] = vi.est[!mask]
      vi.dat$ili[mask] = st.dat$ili[st.dat$epiweek >= firstepiweek]
      vi.dat$ili[vi.dat$ili < 0] = 0
      vi.dat$wili = vi.dat$ili
      vi.dat$region = st
      st.dat = vi.dat
    }
  } else if (st == "la" || st == "fl") {
    if (st == "la") {
      r = c("ok","ar","tx","nm","la")
      rdat = get_completed_fluview_state_df("hhs6",first.week.of.season = first.week.of.season)
    } else {
      r = c("ga","al","ms","tn","ky","nc","sc","fl")
      rdat = get_completed_fluview_state_df("hhs4",first.week.of.season = first.week.of.season)
    }
    st1 = get_completed_fluview_state_df(r[1],first.week.of.season = first.week.of.season)
    st2 = get_completed_fluview_state_df(r[2],first.week.of.season = first.week.of.season)
    x = merge(st1,st2,by="epiweek",suffixes=c("",""))
    for (s in r[3:length(r)-1]) {
      st3 = get_completed_fluview_state_df(s,first.week.of.season = first.week.of.season)
      x = merge(x,st3,by="epiweek",suffixes=c("",""))
    }
    vars = colnames(st1)[startsWith(colnames(st1),"num_")] # List of numeric variables in data frame
    firstepiweek = min(st1$epiweek)
    st.dat = cbind(rdat[(rdat["epiweek"] >= firstepiweek),]) # Default values are that of the region
    st.dat$region = st
    for (v in vars) {
      st.dat[[v]] = get(v,st.dat) - rowSums(x[colnames(x)==v],na.rm=FALSE) # Subtract the sum of the other states in the region
    }
    st.dat$wili = st.dat$num_ili * 100 / st.dat$num_patients # Reset %ILI and wILI
    st.dat$ili = st.dat$wili
  }
  return (st.dat)
}
