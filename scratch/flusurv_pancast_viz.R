
library("pipeR")

## shadowtext does not work without library'ing in ggplot2 & scales
library("ggplot2")
library("scales")

source("../sample/hospitalization-config.R")

epiproject.cache.dir = "~/files/nosync/epiforecast-epiproject/flusurv-pancast-viz"

## s.viz.seasons = s.retro.seasons["2017/2018",drop=FALSE]
s.viz.seasons = s.retro.seasons["2016/2017",drop=FALSE]
## w.viz.model.weeks = w.retro.model.weeks["MW50",drop=FALSE]
## w.viz.model.weeks = w.retro.model.weeks["MW54",drop=FALSE]
w.viz.model.weeks = w.retro.model.weeks["MW58",drop=FALSE]
## w.viz.model.weeks = w.retro.model.weeks["MW60",drop=FALSE]
b.viz.backcasters = b.backcasters["quantile_arx_pancast",drop=FALSE]
f.viz.forecasters = f.forecasters["Delphi_MarkovianDeltaDensity_PackageDefaults",drop=FALSE]

swg.viz.voxel.data =
    map_join(
        get_voxel_data,
        s.viz.seasons, w.viz.model.weeks, g.epigroups,
        last.losocv.issue,
        cache.prefix=file.path(epiproject.cache.dir,"swg.retro.voxel.data"),
        use.proxy=TRUE
    )

## Version of the above grouped by season and model week (an array of arrays of objects):
sw.g.viz.voxel.data = map_join(
    function(swg.array, s,w) {
        swg.array[s,w,,drop=FALSE] %>>%
            select_dims(1:2, "drop") # drop s & w dimensions, but keep g dimension even if size 1
    },
    no_join(swg.viz.voxel.data),
    named_arrayvec_to_name_arrayvec(s.viz.seasons),
    named_arrayvec_to_name_arrayvec(w.viz.model.weeks),
    lapply_variant=lapply, show.progress=FALSE
)

swgb.viz.full.dats = map_join(
    get_backcast,
    swg.viz.voxel.data, sw.g.viz.voxel.data, source.name, signal.name, b.viz.backcasters,
    use.proxy=TRUE,
    cache.prefix=file.path(epiproject.cache.dir,"swgb.viz.full.dats")
)

n.regions = 5L
probs = seq(1/n.regions,by=0.2,length.out=n.regions-1L)

voxel.data = swg.viz.voxel.data[[1L,1L,"Overall"]]
full.dat = swgb.viz.full.dats[[1L,1L,"Overall",1L]]
## voxel.data = swg.viz.voxel.data[[1L,1L,"65+ yr"]]
## full.dat = swgb.viz.full.dats[[1L,1L,"65+ yr",1L]]

provisional.display.traj =
    voxel.data %>>%
    {dplyr::filter(.[["epidata.dfs"]][[source.name]],
                   season == .[["season"]])
    } %>>%
    `[[`(signal.name) %>>%
    flusurv2017_target_trajectory_preprocessor() %>>%
    {.}

previous.issue.display.trajs =
    add_epiweek_integer(voxel.data[["issue"]],-1:-4) %>>%
    setNames(.) %>>%
    lapply(function(issue) {
        g.flusurv.network_all.history.dfs[[voxel.data[["epigroup"]]]] %>>%
            mimicPastEpidataDF(issue) %>>%
            dplyr::filter(!dplyr::between(epiweek%%100L, 18L,39L)) %>>%
            {.[[signal.name]][.[["season"]]==voxel.data[["season"]]]} %>>%
            flusurv2017_target_trajectory_preprocessor() %>>%
            {.}
    })

finalized.display.traj =
    get_observed_trajectory(voxel.data[["season"]], voxel.data[["epigroup"]]) %>>%
    flusurv2017_target_trajectory_preprocessor() %>>%
    {.}

previous.display.traj =
    get_observed_trajectory(voxel.data[["season"]]-1L, voxel.data[["epigroup"]]) %>>%
    flusurv2017_target_trajectory_preprocessor() %>>%
    {.}

pancast.quantile.display.trajs =
    full.dat %>>%
    `[[`(length(.)) %>>%
    {
        sapply(seq_len(nrow(.[["ys"]])), function(time.i) {
            wtd.quantile.or.na(.[["ys"]][time.i,,drop=TRUE], .[["weights"]], probs=probs)
        })
    } %>>%
    apply(1L,flusurv2017_target_trajectory_preprocessor) %>>%
    {.}

## matplot(type="l", full.dat%>>%`[[`(length(.))%>>%(ys)%>>%apply(2L,flusurv2017_target_trajectory_preprocessor))
matplot(type="l", pancast.quantile.display.trajs, ylim=c(0,
                                                         full.dat %>>%
                                                         `[[`(length(.)) %>>%
                                                         {
                                                             sapply(seq_len(nrow(.[["ys"]])), function(time.i) {
                                                                 wtd.quantile.or.na(.[["ys"]][time.i,,drop=TRUE], .[["weights"]], probs=0.95)
                                                             })
                                                         } %>>%
                                                         max(na.rm=TRUE) %>>%
                                                         `*`(1.5)
                                                         ))
lines(provisional.display.traj, lwd=5)
lines(finalized.display.traj, lwd=5)




month_label_locations = function(seasons, start_month_int=10L, end_month_int=5L) {
  seasons %>>%
    lapply(function(season) {
      data.frame(date=seq(
                   as.Date(paste0(season,sprintf("-%02d-01",start_month_int))),
                   as.Date(paste0(season+1L,sprintf("-%02d-31",end_month_int))),
                   by=as.difftime(1L, units="days"))) %>>%
        cbind(rweek={
          n = nrow(.)
          DateToYearWeekWdayDF(.[["date"]][[1L]], 0L,3L) %>>%
            {.[["week"]]+.[["wday"]]/7} %>>%
            magrittr::add(seq(from=0, by=1/7, length.out=n)) %>>%
            magrittr::subtract(0.5) # align to Sunday midnight (from Wednesday mid-day)
        }) %>>%
        tibble::as_tibble() %>>%
        {.}
    }) %>>%
    dplyr::bind_rows() %>>%
    dplyr::mutate(month=lubridate::month(date)) %>>%
    dplyr::group_by(month) %>>%
    dplyr::summarize(rweek=mean(rweek)) %>>%
    dplyr::ungroup() %>>%
    {.}
}
month_sep_locations = function(seasons, start_month_int=10L, end_month_int=5L) {
  seasons %>>%
    lapply(function(loop_season) {
      tibble::tibble(date=
                         paste0(loop_season,sprintf("-%02d-01", start_month_int)) %>>%
                         as.POSIXlt() %>>%
                         {.[["mon"]] <- .[["mon"]] + seq_len(12L+end_month_int+1L-start_month_int+1L)-1L; .} %>>%
                         as.Date()
                     ) %>>%
        cbind(DateToYearWeekWdayDF(.[["date"]], 0L,3L)) %>>%
        dplyr::mutate(dateyear=lubridate::year(date)) %>>%
        dplyr::mutate(relyear=dateyear-loop_season) %>>%
        cbind(yearWeekToSeasonModelWeekDF(.$year, .$week, 31L,3L)) %>>%
        tibble::as_tibble() %>>%
        dplyr::mutate(rweek=(week + wday/7 - 0.5)[1L] + as.integer(date - date[1L])/7) %>>%
        {.}
    }) %>>%
    dplyr::bind_rows() %>>%
    dplyr::mutate(month=lubridate::month(date)) %>>%
    dplyr::group_by(relyear, month) %>>%
    dplyr::summarize(rweek=mean(rweek)) %>>%
    dplyr::ungroup() %>>%
    {.}
}

traj.times = voxel.data %>>%
    {DatesOfSeason(.[["season"]], .[["first.week.of.season"]], 0L,3L)[[1L]]} %>>%
    DateToYearWeekWdayDF(0L,3L) %>>%
    yearWeekDFToSeasonModelWeekDF(31L,3L) %>>%
    (model.week) %>>%
    {.}
forecast.time.t =
    max(traj.times[!is.na(provisional.display.traj)])
forecast.time.i =
    which(traj.times==forecast.time.t)
season.end.i =
    max(which(!is.na(finalized.display.traj)))
month.label.locations = month_label_locations(voxel.data[["season"]])
month.sep.locations = month_sep_locations(voxel.data[["season"]])

data.ymin = 0
## ymax = max(flusurv2017.age.group.percentage.bin.infos[[swg.viz.voxel.data[[1L]][["epigroup"]]]][["breaks"]])
## ymax = 2*max(pancast.quantile.display.trajs, na.rm=TRUE)
ymax = 1.5*max(pancast.quantile.display.trajs, na.rm=TRUE)
## ymax = 1.75*max(pancast.quantile.display.trajs, na.rm=TRUE)
## annot.ymin = data.ymin-0.05*ymax
annot.ymin = data.ymin
width.in = 7.5
height.in = 3.75
xmin = min(traj.times[!is.na(finalized.display.traj)])
## xmax = traj.times[[forecast.time.i]]+diff(traj.times)[[1L]]*30L
## xmax = traj.times[season.end.i]+(2L+15L)*diff(traj.times)[[1L]]
xmax = traj.times[season.end.i]
ribbon.tbl =
    pancast.quantile.display.trajs %>>%
    {list(Low=cbind(data.ymin,.), High=cbind(.,Inf))} %>>%
    ## {list(Low=cbind(data.ymin,.), High=cbind(.,ymax))} %>>%
    lapply(function(mat) {
        mat %>>%
            ## magrittr::set_colnames(LETTERS[seq_len(ncol(.))]) %>>%
            magrittr::set_colnames(sprintf("%04d",seq_len(ncol(.)))) %>>%
            tibble::as_tibble() %>>%
            dplyr::mutate(Time=traj.times) %>>%
            tidyr::gather(Region,Rate,-Time) %>>%
            {.}
    }) %>>%
    dplyr::bind_rows(.id="Division") %>>%
    tidyr::spread(Division, Rate) %>>%
    na.omit() %>>%
    ## dplyr::filter(match(Time,traj.times) <= forecast.time.i+4L) %>>%
    ## dplyr::filter(Time <= traj.times[forecast.time.i+8L]) %>>%
    ## dplyr::bind_rows(
    ##            . %>>%
    ##            dplyr::filter(Time == traj.times[season.end.i]) %>>%
    ##            dplyr::mutate(Time = traj.times[season.end.i]+0.5*diff(traj.times)[[1L]])
    ##           ,
    ##            tibble::tibble(
    ##                        Region = sprintf("%04d",seq_len(n.regions)),
    ##                        ## Time = traj.times[forecast.time.i+8L]+4L*diff(traj.times)[[1L]],
    ##                        Time = traj.times[season.end.i]+2L*diff(traj.times)[[1L]],
    ##                        High = data.ymin + (ymax-data.ymin)*c(seq_len(n.regions-1L),Inf)/n.regions,
    ##                        Low = data.ymin + (ymax-data.ymin)*(seq_len(n.regions)-1L)/n.regions
    ##                    ) %>>%
    ##            dplyr::bind_rows(. %>>% dplyr::mutate(Time=xmax))
    ##        ) %>>%
    {.}
ribbon.annot.tbl = ribbon.tbl %>>%
    dplyr::filter(Time == traj.times[season.end.i]) %>>%
    dplyr::bind_rows(
               . %>>%
               dplyr::mutate(Time = traj.times[season.end.i]+0.5*diff(traj.times)[[1L]])
              ,
               tibble::tibble(
                           Region = sprintf("%04d",seq_len(n.regions)),
                           ## Time = traj.times[forecast.time.i+8L]+4L*diff(traj.times)[[1L]],
                           Time = traj.times[season.end.i]+2L*diff(traj.times)[[1L]],
                           High = data.ymin + (ymax-data.ymin)*c(seq_len(n.regions-1L),Inf)/n.regions,
                           Low = data.ymin + (ymax-data.ymin)*(seq_len(n.regions)-1L)/n.regions
                       ) %>>%
               dplyr::bind_rows(. %>>% dplyr::mutate(Time=xmax))
           ) %>>%
    {.}
point.tbl =
    tibble::tibble(
                Time=traj.times,
                provisional=provisional.display.traj,
                finalized=finalized.display.traj
                ## previous=previous.display.traj
            ) %>>%
    tidyr::gather(Color, Rate, -Time) %>>%
    dplyr::bind_rows(
               previous.issue.display.trajs %>>%
               lapply(function(traj) tibble::tibble(Rate=traj,Time=traj.times)) %>>%
               dplyr::bind_rows(.id="Color")
           ) %>>%
    dplyr::mutate_at(dplyr::vars(Color), dplyr::recode_factor,
                     provisional="Provisional data (as of time prediction was made)",
                     finalized="Finalized data (available after season end)",
                     .ordered=TRUE) %>>%
    na.omit() %>>%
    {.}
line.tbl =
    tibble::tibble(
                Time=traj.times,
                ## provisional=provisional.display.traj,
                ## finalized=finalized.display.traj,
                previous=previous.display.traj
            ) %>>%
    tidyr::gather(Color,Rate,-Time) %>>%
    dplyr::filter(Color!="previous") %>>%
    na.omit() %>>%
    {.}
fill.scale =
    viridis::scale_fill_viridis(
                 ## option="A",
                 ## option="B",
                 ## option="C",
                 option="D",
                 discrete=TRUE)
ribbon.tbl %>>%
    ## dplyr::mutate(alpha=1/(1+High-Low)%>>%magrittr::divide_by(max(.))) %>>%
    ggplot2::ggplot(ggplot2::aes(Time)) %>>%
    `+`(fill.scale) %>>%
    `+`(ggplot2::geom_ribbon(ggplot2::aes(ymin=Low,ymax=High,
                                          ## alpha=alpha,
                                          fill=Region),
                             alpha=0.5
                         )) %>>%
    `+`(ggplot2::geom_line(ggplot2::aes(y=Low, group=Region),
                           ribbon.tbl %>>% dplyr::filter(Region!=sprintf("%04d",1L)),
                           ## linetype="21",
                           ## colour="black"
                           colour="#2A2A2A"
                           )) %>>%
    `+`(ggplot2::geom_line(ggplot2::aes(y=High, group=Region),
                           ribbon.tbl %>>% dplyr::filter(Region!=sprintf("%04d",n.regions)),
                           ## linetype="21",
                           colour="#2A2A2A"
                           )) %>>%
    ## `+`(ggplot2::geom_ribbon(ggplot2::aes(ymin=Low,ymax=High,
    ##                                       fill=Region),
    ##                          ribbon.annot.tbl,
    ##                          alpha=0.5
    ##                          )) %>>%
    ## `+`(ggplot2::geom_line(ggplot2::aes(y=Low, group=Region),
    ##                        ribbon.annot.tbl %>>% dplyr::filter(Region!=sprintf("%04d",1L)),
    ##                        linetype="32",
    ##                        colour="#2A2A2A"
    ##                        )) %>>%
    ## `+`(ggplot2::geom_line(ggplot2::aes(y=High, group=Region),
    ##                        ribbon.annot.tbl %>>% dplyr::filter(Region!=sprintf("%04d",n.regions)),
    ##                        linetype="32",
    ##                        colour="#2A2A2A"
    ##                        )) %>>%
    ## `+`(ggplot2::geom_hline(yintercept=data.ymin)) %>>%
    ## `+`(ggplot2::geom_hline(yintercept=Inf)) %>>%
    `+`(ggplot2::annotate("segment", x=xmin,xend=xmax,y=data.ymin,yend=data.ymin)) %>>%
    `+`(ggplot2::annotate("segment", x=xmin,xend=xmax,y=Inf,      yend=Inf)) %>>%
    ## `+`(ggplot2::geom_vline(
    ##              xintercept=forecast.time.t,
    ##              linetype="32"
    ##          )) %>>%
    `+`(ggplot2::annotate("path",
                     x=c(forecast.time.t+0.5*diff(traj.times)[[1L]],
                         forecast.time.t+0.5*diff(traj.times)[[1L]],
                         forecast.time.t-3*diff(traj.times)[[1L]],
                         forecast.time.t-3*diff(traj.times)[[1L]]
                         ),
                     y=c(annot.ymin,
                         ymax*0.90,
                         ymax*0.95,
                         Inf
                         ),
                     linetype="32"
                 )) %>>%
    `+`(ggplot2::annotate("path",
                     x=c(forecast.time.t+1.5*diff(traj.times)[[1L]],
                         forecast.time.t+1.5*diff(traj.times)[[1L]],
                         forecast.time.t+5*diff(traj.times)[[1L]],
                         forecast.time.t+5*diff(traj.times)[[1L]]
                         ),
                     y=c(annot.ymin,
                         ymax*0.90,
                         ymax*0.95,
                         Inf
                         ),
                     linetype="32"
                 )) %>>%
    ## `+`(ggplot2::geom_hline(
    ##                  yintercept=ymax*0.90,
    ##                  linetype="32"
    ##              )) %>>%
    ## `+`(ggplot2::geom_hline(
    ##                  yintercept=ymax*0.95,
    ##                  ## linetype="32"
    ##                  linetype="solid"
    ##              )) %>>%
    ## `+`(ggplot2::annotate("rect",
    ##                       ymin=ymax*0.90,
    ##                       ymax=ymax*0.95,
    ##                       ## xmin=-Inf,
    ##                       xmin=xmin,
    ##                       ## xmax=Inf,
    ##                       ## xmax=traj.times[season.end.i]+diff(traj.times)[[1L]],
    ##                       ## xmax=traj.times[season.end.i]+0.5*diff(traj.times)[[1L]],
    ##                       xmax=traj.times[season.end.i],
    ##                       fill="black",
    ##                       ## alpha=0.05
    ##                       alpha=0.1
    ##              )) %>>%
    `+`(ggplot2::geom_point(ggplot2::aes(y=Rate),
                            colour="black", size=3.5,
                            point.tbl)) %>>%
    `+`(ggplot2::geom_line(ggplot2::aes(y=Rate, group=Color),
                           colour="black", size=2,
                           point.tbl)) %>>%
    ## ggplot2::geom_point(ggplot2::aes(y=Rate),
    ##                     colour="black", size=3.5,
    ##                     point.tbl)) %>>%
    `+`(ggplot2::geom_point(ggplot2::aes(y=Rate, colour=Color),
                        size=2,
                        point.tbl)) %>>%
    `+`(ggplot2::geom_line(ggplot2::aes(y=Rate, colour=Color),
                       size=1.25,
                       point.tbl)) %>>%
    ## `+`(ggplot2::geom_line(ggplot2::aes(y=Rate, group=Color),
    ##                    colour="black",
    ##                    linetype="41",
    ##                    line.tbl)) %>>%
    `+`(ggplot2::geom_line(ggplot2::aes(y=Rate, colour=Color),
                           line.tbl)) %>>%
    ## ## Revision arrow
    ## `+`(ggplot2::geom_segment(
    ##              arrow=ggplot2::arrow(angle=20, length=ggplot2::unit(0.03,"npc")),
    ##              colour="darkred",
    ##              size=2,
    ##              x=forecast.time.t, xend=forecast.time.t,
    ##              y=provisional.display.traj[[forecast.time.i]]+grid::unit(3.5/2, "mm")%>>%grid::convertY("npc")%>>%as.numeric()*(ymax-annot.ymin)*
    ##                  3,
    ##              yend=finalized.display.traj[[forecast.time.i]]-grid::unit(3.5/2+2/2, "mm")%>>%grid::convertY("npc")%>>%as.numeric()*(ymax-annot.ymin)*
    ##                  3
    ##          )) %>>%
    ## `+`(ggplot2::geom_segment(
    ##              arrow=ggplot2::arrow(angle=20, length=ggplot2::unit(0.03,"npc")),
    ##              colour="#FFAAAA",
    ##              x=forecast.time.t, xend=forecast.time.t,
    ##              y=provisional.display.traj[[forecast.time.i]]+grid::unit(3.5/2, "mm")%>>%grid::convertY("npc")%>>%as.numeric()*(ymax-annot.ymin)*
    ##                  3,
    ##              yend=finalized.display.traj[[forecast.time.i]]-grid::unit(3.5/2+2/2, "mm")%>>%grid::convertY("npc")%>>%as.numeric()*(ymax-annot.ymin)*
    ##                  3
    ##          )) %>>%
    ## `+`(ggplot2::geom_segment(
    ##              arrow=ggplot2::arrow(angle=20, length=ggplot2::unit(0.03,"npc")),
    ##              colour="darkblue",
    ##              size=2,
    ##              x=traj.times[[20L]], xend=traj.times[[20L]],
    ##              y=previous.display.traj[[20L]],
    ##              yend=finalized.display.traj[[20L]]-grid::unit(3.5/2, "mm")%>>%grid::convertY("npc")%>>%as.numeric()*(ymax-annot.ymin)*
    ##                  3
    ##              )) %>>%
    ## `+`(ggplot2::geom_segment(
    ##              arrow=ggplot2::arrow(angle=20, length=ggplot2::unit(0.03,"npc")),
    ##              colour="skyblue",
    ##              x=traj.times[[20L]], xend=traj.times[[20L]],
    ##              y=previous.display.traj[[20L]],
    ##              yend=finalized.display.traj[[20L]]-grid::unit(3.5/2, "mm")%>>%grid::convertY("npc")%>>%as.numeric()*(ymax-annot.ymin)*
    ##                  3
    ##              )) %>>%
    `+`(ggplot2::annotate("text",
                          label="Nowcast",
                          ## forecast.time.t, ymax)) %>>%
                          forecast.time.t+diff(traj.times)[[1L]], ymax)) %>>%
    `+`(ggplot2::annotate("text",
                      ## label="Past",
                      ## hjust=1.3,
                      label="Backcasts",
                      hjust=1.2,
                      forecast.time.t-3*diff(traj.times)[[1L]], ymax)) %>>%
    ## `+`(ggplot2::geom_segment(
    ##              arrow=ggplot2::arrow(angle=20, length=ggplot2::unit(0.03,"npc")),
    ##              ## x=forecast.time.t-1.6*as.numeric(grid::convertX(grid::unit(strwidth("Past","inches"),"inches"),"npc"))*(xmax-xmin),
    ##              ## xend=forecast.time.t-3.6*as.numeric(grid::convertX(grid::unit(strwidth("Past","inches"),"inches"),"npc"))*(xmax-xmin),
    ##              x=forecast.time.t-2.7*1.6*as.numeric(grid::convertX(grid::unit(strwidth("Past","inches"),"inches"),"npc"))*(xmax-xmin),
    ##              xend=forecast.time.t-4.7*1.6*as.numeric(grid::convertX(grid::unit(strwidth("Past","inches"),"inches"),"npc"))*(xmax-xmin),
    ##              y=ymax, yend=ymax
    ##          )) %>>%
    `+`(ggplot2::annotate("text",
                      ## label="Future",
                      ## hjust=-0.2,
                      label="Forecasts",
                      hjust=-0.2,
                      forecast.time.t+5*diff(traj.times)[[1L]], ymax)) %>>%
    ## `+`(ggplot2::geom_segment(
    ##              arrow=ggplot2::arrow(angle=20, length=ggplot2::unit(0.03,"npc")),
    ##              ## x=forecast.time.t+1.4*as.numeric(grid::convertX(grid::unit(strwidth("Future","inches"),"inches"),"npc"))*(xmax-xmin),
    ##              ## xend=forecast.time.t+3.4*as.numeric(grid::convertX(grid::unit(strwidth("Future","inches"),"inches"),"npc"))*(xmax-xmin),
    ##              ## x=forecast.time.t+3*1.4*as.numeric(grid::convertX(grid::unit(strwidth("Future","inches"),"inches"),"npc"))*(xmax-xmin),
    ##              ## xend=forecast.time.t+5*1.4*as.numeric(grid::convertX(grid::unit(strwidth("Future","inches"),"inches"),"npc"))*(xmax-xmin),
    ##              x=forecast.time.t+2.7*1.4*as.numeric(grid::convertX(grid::unit(strwidth("Future","inches"),"inches"),"npc"))*(xmax-xmin),
    ##              xend=forecast.time.t+(2.7*1.4*as.numeric(grid::convertX(grid::unit(strwidth("Future","inches"),"inches"),"npc")) + (4.7-2.7)*1.6*as.numeric(grid::convertX(grid::unit(strwidth("Past","inches"),"inches"),"npc")))*(xmax-xmin),
    ##              y=ymax, yend=ymax
    ##          )) %>>%
    ## `+`(ggplot2::scale_colour_manual(values=c("grey","green"))) %>>%
    ## `+`(ggplot2::annotate("text",
    ## `+`(shadowtext::geom_shadowtext(
    ##                 ggplot2::aes(x=NULL),
    ##                 tibble::tibble(arbitrary=1), # ensure text is only added once
    ##                 colour="skyblue",
    ##                 label="Change from last season",
    ##                 ## todo better wording than "change"? "up" or "down" based on delta amount?
    ##                 hjust=-0.05,
    ##                 ## hjust=1.05, colour="#FFAAAA",
    ##                 x=traj.times[[20L]],
    ##                 y=(previous.display.traj[[20L]]+finalized.display.traj[[20L]])/2)) %>>%
    ## `+`(ggplot2::annotate("text",
    ##                   label="Actual revision",
    ##                   hjust=-0.05,
    ##                   traj.times[[forecast.time.i]],
    ##                   (provisional.display.traj[[forecast.time.i]]+finalized.display.traj[[forecast.time.i]])/2)) %>>%
    ## `+`(shadowtext::geom_shadowtext(
    ##                 ggplot2::aes(x=NULL),
    ##                 tibble::tibble(arbitrary=1), # ensure text is only added once
    ##                 colour="#FFAAAA",
    ##                 label="Actual revision",
    ##                 hjust=-0.05,
    ##                 x=traj.times[[forecast.time.i]],
    ##                 y=(provisional.display.traj[[forecast.time.i]]+finalized.display.traj[[forecast.time.i]])/2,
    ##                 )) %>>%
    `+`(shadowtext::geom_shadowtext(
                    ggplot2::aes(x=NULL),
                    tibble::tibble(arbitrary=1), # ensure text is only added once
                    colour="grey",
                    label="Provisional",
                    vjust=1.5,
                    x=traj.times[[forecast.time.i]],
                    y=provisional.display.traj[[forecast.time.i]]
                    )) %>>%
    `+`(shadowtext::geom_shadowtext(
                    ggplot2::aes(x=NULL),
                    tibble::tibble(arbitrary=1), # ensure text is only added once
                    colour="green",
                    label="Finalized",
                    vjust=-0.5,
                    x=traj.times[[forecast.time.i]],
                    y=finalized.display.traj[[forecast.time.i+2L]]
                )) %>>%
    ## `+`(shadowtext::geom_shadowtext(
    ## `+`(ggplot2::geom_label(
    ##              ggplot2::aes(x=Time,y=(Low+High)/2),
    ##              ## ggplot2::aes(x=Time,y=(Low+High)/2,fill=Region),
    ##              ribbon.tbl %>>%
    ##              dplyr::filter(Time==traj.times[forecast.time.i+8L]) %>>%
    ##              dplyr::mutate(High=dplyr::if_else(is.infinite(High),ymax,High)) %>>%
    ##              dplyr::arrange(Region),
    ##              ## colour=fill.scale[["palette"]](n.regions),
    ##              ## bg.colour=fill.scale[["palette"]](n.regions),
    ##              ## alpha=0.2,
    ##              label=sprintf("Est. 1 in %d chance that finalized point falls here",n.regions),
    ##              ## hjust=-0.01
    ##              hjust=+0.01
    ##              ## hjust=0
    ##          )) %>>%
    ## `+`(ggplot2::geom_label(
    ##                 ## ggplot2::aes(x=Time,y=(Low+High)/2),
    ##                 ggplot2::aes(x=Time,y=(Low+High)/2,fill=Region),
    ##                 ribbon.tbl %>>%
    ##                 dplyr::filter(Time==traj.times[forecast.time.i+8L]) %>>%
    ##                 dplyr::mutate(High=dplyr::if_else(is.infinite(High),ymax,High)) %>>%
    ##                 dplyr::arrange(Region),
    ##                 colour="transparent",
    ##                 ## colour=fill.scale[["palette"]](n.regions),
    ##                 ## bg.colour=fill.scale[["palette"]](n.regions),
    ##                 alpha=0.2,
    ##                 label=sprintf("Est. 1 in %d chance that finalized point falls here",n.regions),
    ##                 ## hjust=-0.01
    ##                 hjust=+0.01
    ##                 ## hjust=0
    ##                 )) %>>%
    ## `+`(ggplot2::geom_text(
    ##              ggplot2::aes(x=Time,y=(Low+High)/2),
    ##              ## ribbon.tbl %>>%
    ##              ribbon.annot.tbl %>>%
    ##              dplyr::filter(Time==traj.times[season.end.i]+2L*diff(traj.times)[[1L]]) %>>%
    ##              dplyr::mutate(High=dplyr::if_else(is.infinite(High),ymax,High)) %>>%
    ##              dplyr::arrange(Region),
    ##              label=sprintf("Estimated 1 in %d chance that\nfinalized point falls here",n.regions),
    ##              hjust=-0.02
    ##          )) %>>%
    ## `+`(ggplot2::geom_text(
    ##     ggplot2::aes(x=rweek, y=0, label=month.abb[month], group=NULL, colour=NULL),
    ##     month.label.locations,
    ##     hjust=0.5,
    ##     vjust=1.4,
    ##     show.legend=FALSE
    ## )) %>>%
    {
        plt = .
        for (row.i in seq_len(nrow(month.label.locations))) {
            ## plt <- plt + ggplot2::annotate("text",
            plt <- plt + ggplot2::annotation_custom(
                                      grob=grid::textGrob(
                                                     hjust=0.5,
                                                     vjust=1.4,
                                                     label=month.abb[[month.label.locations[[row.i,"month"]]]],
                                                     ),
                                      xmin=month.label.locations[[row.i,"rweek"]]-10,
                                      xmax=month.label.locations[[row.i,"rweek"]]+10,
                                      ymin=0-10,
                                      ymax=0+10
                                  )
        }
        plt
    } %>>%
    ## `+`(ggplot2::geom_vline(xintercept=traj.times[season.end.i]+diff(traj.times)[[1L]],
    ##                         linetype="32")) %>>%
    ## `+`(ggplot2::geom_vline(xintercept=traj.times[season.end.i]+2L*diff(traj.times)[[1L]],
    ##                         ## linetype="32")) %>>%
    ##                         linetype="solid")) %>>%
    ## `+`(ggplot2::annotate("rect",
    ##                       xmin=traj.times[season.end.i]+0.5*diff(traj.times)[[1L]],
    ##                       xmax=traj.times[season.end.i]+2L*diff(traj.times)[[1L]],
    ##                       ymin=data.ymin,
    ##                       ymax=Inf,
    ##                       fill="black",
    ##                       ## alpha=0.05
    ##                       alpha=0.1
    ##                       )) %>>%
    ## `+`(ggplot2::expand_limits(y=c(annot.ymin,ymax))) %>>%
    `+`(ggplot2::coord_cartesian(clip="off")) %>>%
    `+`(ggplot2::expand_limits(y=c(data.ymin,ymax))) %>>%
    `+`(ggplot2::expand_limits(x=c(xmin,xmax))) %>>%
    `+`(ggplot2::theme(legend.position="bottom")) %>>%
    `+`(ggplot2::guides(fill=FALSE)) %>>%
    `+`(ggplot2::scale_x_continuous(
                 labels=function(breaks) model_week_to_epi_week(breaks,31L,lastWeekNumber(voxel.data[["season"]], 3L)),
                 breaks=round(month.sep.locations[["rweek"]]),
                 minor_breaks=round(month.sep.locations[["rweek"]])
    )) %>>%
    `+`(ggplot2::ylab("Rate")) %>>%
    `+`(ggplot2::ggsave("pancast_banded.png", width=width.in, height=height.in, units="in", dpi=300))
    ## `+`(ggplot2::ggsave("pancast_banded.png", width=width.in, height=height.in, units="in", dpi=300, type="cairo"))




## todo filter to top group
## todo rearrange legend to match visual order
## todo region labels on plot
## todo line labels on plot
## todo month labels, epiweek axis, ...
## todo if do black line outlines, adjust legend to match.  See if there is a built-in command to do this outlining
## todo better color scale
## todo region explanations
## todo add individual point annotations? state quantile level somehow?
## todo relating any/all of the above and/or something else to historical levels?
## todo add historical median?
## todo annotated arrow from provisional up to revised
## todo arrows left and right near past & future labels





## todo series label on graph
## todo to introduce nowcast label, potentially allow full forecast trajectory display, and clean up right-hand side by reducing overlap and perhaps allowing multi-line labels: have transition from narrow colored/outlined regions to wider labeling areas through some linear transition


## change to something unbiased, but with large revisions
## last four issues in different colors + finalized in black
## put in terms of issues
## bracket on top
##       Pancasting
##   _________^__________
##   |                  |
##   back,now,forecasting
## remove RHS labeling
## discussion on bottom
