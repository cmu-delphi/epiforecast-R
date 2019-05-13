
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
## n.regions = 4L
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

n.older.lags = 4L
older.lags = seq_len(n.older.lags)
previous.issue.display.trajs =
    add_epiweek_integer(voxel.data[["issue"]], -older.lags) %>>%
    setNames(sprintf("%d issue (%d week%s old)",
                     ., older.lags, dplyr::recode(older.lags,`1`="",.default="s"))) %>>%
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
fill.ymax = scales::pretty_breaks()(c(data.ymin,0.8*ymax)) %>>% `[[`(length(.))
## ymax = 1.75*max(pancast.quantile.display.trajs, na.rm=TRUE)
## annot.ymin = data.ymin-0.05*ymax
annot.ymin = data.ymin
width.in = 7.5
## height.in = 3.75
height.in = 4
## height.in = 5
xmin = min(traj.times[!is.na(finalized.display.traj)])
## xmax = traj.times[[forecast.time.i]]+diff(traj.times)[[1L]]*30L
## xmax = traj.times[season.end.i]+(2L+15L)*diff(traj.times)[[1L]]
xmax = traj.times[season.end.i]
ribbon.tbl =
    pancast.quantile.display.trajs %>>%
    ## {list(Low=cbind(data.ymin,.), High=cbind(.,Inf))} %>>%
    ## {list(Low=cbind(data.ymin,.), High=cbind(.,ymax))} %>>%
    {list(Low=cbind(data.ymin,.), High=cbind(.,fill.ymax))} %>>%
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
previous.label = paste0(voxel.data[["issue"]]," issue (latest available at prediction time)")
## finalized.label = "Finalized data (available after season end)"
finalized.label = "Finalized data (unavailable until after season end)"
traj.tbl =
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
    dplyr::mutate(Type = dplyr::recode_factor(c("point","thickline","line", Color),
                                              point="point", thickline="thickline",line="line",
                                              finalized="point",
                                              provisional="thickline",
                                              .default="thickline")[-(1:3),drop=FALSE]
                  ) %>>%
    dplyr::mutate(shortname=Color) %>>%
    dplyr::mutate_at(dplyr::vars(Color), ordered,
                     c(rev(names(previous.issue.display.trajs)),"provisional","finalized")) %>>%
    dplyr::mutate_at(dplyr::vars(Color), dplyr::recode,
                     ## provisional="Provisional data (as of time prediction was made)",
                     provisional=previous.label,
                     finalized=finalized.label) %>>%
    na.omit() %>>%
    ## dplyr::arrange(Color) %>>%
    ## dplyr::mutate_at(dplyr::vars(Color), as.character) %>>%
    {.}
traj.tbls = split(traj.tbl, traj.tbl[["Type"]], drop=FALSE)
point.tbl = traj.tbls[["point"]]
thickline.tbl = traj.tbls[["thickline"]]
line.tbl = traj.tbls[["line"]]
n.thickline.trajs = length(unique(thickline.tbl[["Color"]]))
colour.fill.hue.sep = 0.2
colour.scale =
    ## viridis::scale_colour_viridis(
    ##              ## option="A",
    ##              ## option="B",
    ##              option="C",
    ##              ## option="D",
    ##              ## option="E",
    ##              ## begin=0,end=(1-colour.fill.hue.sep)*(n.thickline.trajs)/(n.thickline.trajs+n.regions+1L),
    ##              begin=(1-colour.fill.hue.sep)*(n.thickline.trajs)/(n.thickline.trajs+n.regions+1L),end=0,
    ##              breaks=levels(traj.tbl[["Color"]]),
    ##              ## labels=identity,
    ##              ## labels=as.character,
    ##              discrete=TRUE)
    ggplot2::scale_colour_manual(
                 values=viridisLite::viridis(
                                         ## option="A",
                                         ## option="B",
                                         ## option="C",
                                         option="D",
                                         ## option="E",
                                         ## begin=0,end=(1-colour.fill.hue.sep)*(n.thickline.trajs)/(n.thickline.trajs+n.regions+1L),
                                         begin=(1-colour.fill.hue.sep)*(n.thickline.trajs)/(n.thickline.trajs+n.regions+1L),end=0,
                                         n.thickline.trajs
                                     ) %>>% c("black"),
                 breaks=levels(traj.tbl[["Color"]])
    )
fill.scale =
    viridis::scale_fill_viridis(
                 ## option="A",
                 ## option="B",
                 ## option="C",
                 option="D",
                 ## option="E",
                 begin=colour.fill.hue.sep+(1-colour.fill.hue.sep)*(n.thickline.trajs+1L)/(n.thickline.trajs+n.regions+1L),end=1,
                 ## begin=1,end=(n.thickline.trajs+1L)/(n.thickline.trajs+n.regions+1L),
                 labels=sprintf(
                     ## "Estimated %s percentile to %s percentile",
                     ## "Estimated %s to %s percentiles",
                     ## "Estimated %s to %s percentile range",
                     ## "Predicted %s to %s percentile range",
                     "%s to %s percentile range",
                     scales::ordinal(100/n.regions*(seq_len(n.regions)-1L)),
                     scales::ordinal(100/n.regions*seq_len(n.regions))
                 ),
                 discrete=TRUE)
plt = ribbon.tbl %>>%
    ## dplyr::mutate(alpha=1/(1+High-Low)%>>%magrittr::divide_by(max(.))) %>>%
    ggplot2::ggplot(ggplot2::aes(Time)) %>>%
    `+`(colour.scale) %>>%
    `+`(fill.scale) %>>%
    `+`(ggplot2::annotate("rect", xmin=-Inf, xmax=Inf, ymin=fill.ymax, ymax=Inf, fill="white")) %>>%
    `+`(ggplot2::geom_ribbon(ggplot2::aes(ymin=Low,ymax=High,
                                          ## alpha=alpha,
                                          fill=Region),
                             colour="#2A2A2A",
                             alpha=0.5
                         )) %>>%
    ## `+`(ggplot2::geom_line(ggplot2::aes(y=Low, group=Region),
    ##                        ribbon.tbl %>>% dplyr::filter(Region!=sprintf("%04d",1L)),
    ##                        ## linetype="21",
    ##                        ## colour="black"
    ##                        colour="#2A2A2A"
    ##                        )) %>>%
    ## `+`(ggplot2::geom_line(ggplot2::aes(y=High, group=Region),
    ##                        ribbon.tbl %>>% dplyr::filter(Region!=sprintf("%04d",n.regions)),
    ##                        ## linetype="21",
    ##                        colour="#2A2A2A"
    ##                        )) %>>%
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
    ## `+`(ggplot2::annotate("segment", x=xmin,xend=xmax,y=Inf,      yend=Inf)) %>>%
    `+`(ggplot2::annotate("segment", x=xmin,xend=xmax,y=fill.ymax,yend=fill.ymax)) %>>%
    ## `+`(ggplot2::geom_vline(
    ##              xintercept=forecast.time.t,
    ##              linetype="32"
    ##          )) %>>%
    {
        ## left.path.xs = c(forecast.time.t+0.5*diff(traj.times)[[1L]],
        ##                  forecast.time.t+0.5*diff(traj.times)[[1L]],
        ##                  ## forecast.time.t-3*diff(traj.times)[[1L]],
        ##                  ## forecast.time.t-3*diff(traj.times)[[1L]]
        ##                  forecast.time.t-1*diff(traj.times)[[1L]],
        ##                  forecast.time.t-1*diff(traj.times)[[1L]]
        ##                  )
        ## right.path.xs = c(forecast.time.t+1.5*diff(traj.times)[[1L]],
        ##                   forecast.time.t+1.5*diff(traj.times)[[1L]],
        ##                   ## forecast.time.t+5*diff(traj.times)[[1L]],
        ##                   ## forecast.time.t+5*diff(traj.times)[[1L]]
        ##                   forecast.time.t+3*diff(traj.times)[[1L]],
        ##                   forecast.time.t+3*diff(traj.times)[[1L]]
        ##                   )
        ## path.ys = c(annot.ymin,
        ##             ## ymax*0.90,
        ##             ## ymax*0.95,
        ##             ## Inf
        ##             ## ymax*0.80,
        ##             ## ymax*0.85,
        ##             ## ymax*0.95
        ##             fill.ymax,
        ##             fill.ymax+0.03*ymax,
        ##             fill.ymax+0.11*ymax
        ##             )
        left.path.xs = c(forecast.time.t+0.5*diff(traj.times)[[1L]],
                         forecast.time.t+0.5*diff(traj.times)[[1L]]
                         )
        right.path.xs = c(forecast.time.t+1.5*diff(traj.times)[[1L]],
                          forecast.time.t+1.5*diff(traj.times)[[1L]]
                          )
        path.ys = c(annot.ymin,
                    fill.ymax
                    )
        . %>>%
            `+`(ggplot2::annotate("path", x=left.path.xs, y=path.ys, linetype="32")) %>>%
            `+`(ggplot2::annotate("path", x=right.path.xs, y=path.ys, linetype="32")) %>>%
            ## `+`(ggplot2::annotate("polygon",
            ##                       x=c(left.path.xs[-1L]+0.1,rev(right.path.xs[-1L])-0.1),
            ##                       y=c(path.ys[-1L],rev(path.ys[-1L]))+ymax*0.005,
            ##                       fill="gray",
            ##                       linetype="32")) %>>%
            ## `+`(ggforce::geom_shape(
            ##                  data=data.frame(arbitrary=rep(1,2L*(length(left.path.xs)-1L))),
            ##                  x=c(left.path.xs[-1L],rev(right.path.xs[-1L])),
            ##                  y=c(path.ys[-1L],rev(path.ys[-1L])),
            ##                  expand=grid::unit(-0.002,"npc"),
            ##                  fill="gray",
            ##                  linetype="32")) %>>%
            ## `+`(ggforce::geom_shape(
            ##                  data=data.frame(arbitrary=rep(1,length(left.path.xs)-1L+2L)),
            ##                  x=c(left.path.xs[-1L],xmin,xmin),
            ##                  y=c(path.ys[-1L],path.ys[c(4L,2L)]),
            ##                  expand=grid::unit(-0.002,"npc"),
            ##                  fill="gray",
            ##                  linetype="32")) %>>%
            ## `+`(ggforce::geom_shape(
            ##                  data=data.frame(arbitrary=rep(1,length(right.path.xs)-1L+2L)),
            ##                  x=c(right.path.xs[-1L],xmax,xmax),
            ##                  y=c(path.ys[-1L],path.ys[c(4L,2L)]),
            ##                  expand=grid::unit(-0.002,"npc"),
            ##                  fill="gray",
            ##                  linetype="32")) %>>%
            ## `+`(ggforce::geom_shape(
            ##                  data=data.frame(arbitrary=rep(1,4L)),
            ##                  x=c(xmin,xmax,xmax,xmin),
            ##                  y=c(path.ys[c(4L,4L)],Inf,Inf),
            ##                  expand=grid::unit(-0.002,"npc"),
            ##                  fill="gray",
            ##                  linetype="32")) %>>%
            {.}
    } %>>%
    {
        plus_bracket = function(plt, leftx, rightx, ybase, midx=0.5*(leftx+rightx)) {
            edge.dx = 0.2
            mid.dx = 0.1
            plt + ggplot2::annotate("path",
                                    x=c(leftx,leftx+edge.dx,midx-mid.dx,midx,midx+mid.dx,rightx-edge.dx,rightx),
                                    ## y=ybase+c(0,3,3,5,3,3,0)*0.005*ymax
                                    y=ybase+c(0,3,3,5,3,3,0)*0.003*ymax
                                    )
        }
        pancast.ybase = fill.ymax+ymax-(fill.ymax+0.065*ymax)
        . %>>%
            ## `+`(ggplot2::annotate("path",
            ##                       x=c(xmin, xmin+0.1,
            ##                           forecast.time.t+1-0.5-0.1, forecast.time.t+1-0.5
            ##                           ),
            ##                       y=c(fill.ymax, fill.ymax+0.01*ymax, fill.ymax+0.01*ymax, fill.ymax)
            ##                       )) %>>%
            ## `+`(ggplot2::annotate("path",
            ##                       x=c(forecast.time.t+1-0.5,forecast.time.t+1-0.5+0.1,forecast.time.t+1+0.5-0.1,forecast.time.t+1+0.5),
            ##                       y=c(fill.ymax, fill.ymax+0.01*ymax, fill.ymax+0.01*ymax, fill.ymax)
            ##                       )) %>>%
            ## `+`(ggplot2::annotate("path",
            ##                       x=c(forecast.time.t+1+0.5,forecast.time.t+1+0.5+0.1,xmax-0.1,xmax),
            ##                       y=c(fill.ymax, fill.ymax+0.01*ymax, fill.ymax+0.01*ymax, fill.ymax)
            ##                       )) %>>%
            plus_bracket(xmin, forecast.time.t+1-0.5, fill.ymax) %>>%
            plus_bracket(forecast.time.t+1-0.5, forecast.time.t+1+0.5, fill.ymax) %>>%
            plus_bracket(forecast.time.t+1+0.5, xmax, fill.ymax) %>>%
            ## plus_bracket(xmin, xmax, pancast.ybase, forecast.time.t+1) %>>%
            plus_bracket(xmin, xmax, pancast.ybase, 0.5*(xmin+xmax)) %>>%
            `+`(ggplot2::annotate("segment",
                                  x=xmin, xend=xmin,
                                  y=fill.ymax, yend=pancast.ybase,
                                  linetype="32")) %>>%
            `+`(ggplot2::annotate("segment",
                                  x=xmax, xend=xmax,
                                  y=fill.ymax, yend=pancast.ybase,
                                  linetype="32")) %>>%
            {.}
    } %>>%
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
    ## `+`(ggplot2::geom_line(ggplot2::aes(y=Rate, group=Color),
    ##                        ## colour="black",
    ##                        colour="#AAAAAA",
    ##                        size=2,
    ##                        thickline.tbl)) %>>%
    `+`(ggplot2::geom_line(ggplot2::aes(y=Rate, colour=Color),
                           size=1.25,
                           thickline.tbl)) %>>%
    ## `+`(ggplot2::geom_point(ggplot2::aes(y=Rate),
    ##                         ## colour="black",
    ##                         colour="#AAAAAA",
    ##                         size=3.5,
    ##                         point.tbl)) %>>%
    ## `+`(ggplot2::geom_line(ggplot2::aes(y=Rate, group=Color),
    ##                        ## colour="black",
    ##                        colour="#AAAAAA",
    ##                        size=2,
    ##                        point.tbl)) %>>%
    ## ggplot2::geom_point(ggplot2::aes(y=Rate),
    ##                     colour="black", size=3.5,
    ##                     point.tbl)) %>>%
    ## `+`(ggplot2::geom_point(ggplot2::aes(y=Rate, colour=Color),
    ##                     size=2,
    ##                     point.tbl)) %>>%
    `+`(ggplot2::geom_line(ggplot2::aes(y=Rate, colour=Color),
                       size=1.25,
                       point.tbl)) %>>%
    ## `+`(ggplot2::geom_line(ggplot2::aes(y=Rate, group=Color),
    ##                    colour="black",
    ##                    linetype="41",
    ##                    line.tbl)) %>>%
    ## `+`(ggplot2::geom_line(ggplot2::aes(y=Rate, colour=Color),
    ##                        line.tbl)) %>>%
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
                          label="Pancasts",
                          ## forecast.time.t, ymax)) %>>%
                          ## forecast.time.t+diff(traj.times)[[1L]],
                          0.5*(xmin+xmax),
                          ## ymax
                          ymax-0.015*ymax
                          )) %>>%
    `+`(ggplot2::annotate("text",
                          label="Nowcast",
                          ## forecast.time.t, ymax)) %>>%
                          forecast.time.t+diff(traj.times)[[1L]],
                          hjust=0.7,
                          ## ymax*0.9
                          ## fill.ymax+0.065*ymax
                          fill.ymax+0.050*ymax
                          )) %>>%
    `+`(ggplot2::annotate("text",
                      ## label="Past",
                      ## hjust=1.3,
                      label="Backcasts",
                      ## hjust=1.2,
                      ## forecast.time.t-3*diff(traj.times)[[1L]],
                      (xmin+forecast.time.t+1-0.5)/2,
                      ## ymax*0.9
                      ## fill.ymax+0.065*ymax
                      fill.ymax+0.050*ymax
                      )) %>>%
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
                      ## hjust=-0.2,
                      ## forecast.time.t+5*diff(traj.times)[[1L]],
                      (forecast.time.t+1+0.5+xmax)/2,
                      ## ymax*0.9
                      ## fill.ymax+0.065*ymax
                      fill.ymax+0.050*ymax
                      )) %>>%
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
                        bg.color="#AAAAAA",
                        ## `+`(ggplot2::geom_text(
                        ## ggplot2::aes(x=NULL,colour=Color),
                        ggplot2::aes(x=NULL),
                        tibble::tibble(Color=names(previous.issue.display.trajs)[[4L]]),
                        ## xxx or just use colour=shortname + guide naming in scale
                        ## colour="grey",
                        colour= colour.scale[["palette"]](length(colour.scale[["breaks"]]))[match(names(previous.issue.display.trajs)[[4L]],colour.scale[["breaks"]])],
                        label=sprintf('bold("%s week%s old")',n.older.lags,if(n.older.lags==1L)""else"s"), parse=TRUE,
                        ## vjust=1.5,
                        vjust=1.3,
                        hjust=0.2,
                        ## hjust=-0.1,
                        x=traj.times[[forecast.time.i-n.older.lags]],
                        y=previous.issue.display.trajs[[n.older.lags]][[forecast.time.i-n.older.lags]]
                    )) %>>%
    `+`(shadowtext::geom_shadowtext(
                        bg.color="#AAAAAA",
    ## `+`(ggplot2::geom_text(
                    ## ggplot2::aes(x=NULL,colour=Color),
                    ggplot2::aes(x=NULL),
                    tibble::tibble(Color=previous.label),
                    ## xxx or just use colour=shortname + guide naming in scale
                    ## colour="grey",
                    colour= colour.scale[["palette"]](length(colour.scale[["breaks"]]))[match(previous.label,colour.scale[["breaks"]])],
                    ## label="Latest provisional",
                    ## label="Latest",
                    label='bold("Latest")', parse=TRUE,
                    ## vjust=1.5,
                    vjust=1.3,
                    hjust=0.2,
                    ## hjust=-0.1,
                    x=traj.times[[forecast.time.i]],
                    y=provisional.display.traj[[forecast.time.i]]
                    )) %>>%
    `+`(shadowtext::geom_shadowtext(
                        bg.color="#AAAAAA",
    ## `+`(ggplot2::geom_text(
                    ## ggplot2::aes(x=NULL,colour=Color),
                    ggplot2::aes(x=NULL),
                    tibble::tibble(Color=finalized.label),
                    ## colour="green",
                    colour= colour.scale[["palette"]](length(colour.scale[["breaks"]]))[match(finalized.label,colour.scale[["breaks"]])],
                    ## label="Finalized",
                    label='bold("Finalized")', parse=TRUE,
                    vjust=-0.5,
                    x=traj.times[[forecast.time.i+2L]],
                    y=finalized.display.traj[[forecast.time.i+2L]]
                )) %>>%
    ## `+`(ggplot2::annotate("text",
    {
        ## dt = 8L
        ## dt = 6L
        ## dt = 4L
        dt = 5L
        . %>>%
            `+`(ggplot2::annotate("text",
            ## `+`(ggrepel::geom_text_repel(ggplot2::aes(x=NULL), data.frame(arbitrary=1), force=0,
            ##                              arrow=arrow(length=unit(0.02, "npc")),
            ##                              ## nudge_y=0.2*ymax,
            ##                              ## nudge_y=0.1*ymax,
            ##                              nudge_y=0.14*ymax,
            ##                              ## nudge_x=2*diff(traj.times)[[1L]],
                                  ## label=sprintf("Est. %s percentile", scales::ordinal(100*probs[[length(probs)]])),
                                  ## label=sprintf("Pred. %s percentile", scales::ordinal(100*probs[[length(probs)]])),
                                  label=sprintf("Predicted\n%s percentile", scales::ordinal(100*probs[[length(probs)]])),
                                  lineheight=1.0,
                                  ## vjust=-0.05,
                                  ## vjust=-0.15,
                                  vjust=-0.2,
                                  ## hjust=-0.05,
                                  ## box.padding=0.25+1.00,
                                  ## hjust=0.40,
                                  ## hjust=0.30,
                                  ## hjust=0.20,
                                  ## hjust=0.45,
                                  hjust=0.35,
                                  colour="#2A2A2A",
                                  x=forecast.time.t+dt*diff(traj.times)[[1L]],
                                  ## y=pancast.quantile.display.trajs%>>%`[[`(forecast.time.i+dt,ncol(.))
                                  ## y=0.14*ymax+pancast.quantile.display.trajs%>>%`[[`(forecast.time.i+dt,ncol(.))
                                  y=0.09*ymax+pancast.quantile.display.trajs%>>%`[[`(forecast.time.i+dt,ncol(.))
                                  )) %>>%
            `+`(ggplot2::annotate("segment",
                                  xend=forecast.time.t+dt*diff(traj.times)[[1L]],
                                  x=forecast.time.t+dt*diff(traj.times)[[1L]],
                                  yend=pancast.quantile.display.trajs%>>%`[[`(forecast.time.i+dt,ncol(.)),
                                  ## y=0.14*ymax+pancast.quantile.display.trajs%>>%`[[`(forecast.time.i+dt,ncol(.)),
                                  y=0.09*ymax+pancast.quantile.display.trajs%>>%`[[`(forecast.time.i+dt,ncol(.)),
                                  arrow=arrow(length=unit(0.02, "npc")),
                                  )) %>>%
            `+`(ggplot2::annotate("text",
            ## `+`(ggrepel::geom_text_repel(ggplot2::aes(x=NULL), data.frame(arbitrary=1), force=0,
            ##                              arrow=arrow(length=unit(0.02, "npc")),
            ##                              ## nudge_y=-0.08*ymax,
            ##                              ## nudge_y=-0.18*ymax,
            ##                              nudge_y=-0.11*ymax,
            ##                              ## nudge_x=-diff(traj.times)[[1L]],
                                         ## label=sprintf("Est. %s percentile", scales::ordinal(100*probs[[1L]])),
                                  ## label=sprintf("Pred. %s percentile", scales::ordinal(100*probs[[1L]])),
                                  label=sprintf("Predicted %s percentile", scales::ordinal(100*probs[[1L]])),
                                         vjust=1.05,
                                         ## hjust=1.05,
                                         ## hjust=0.7,
                                         hjust=0.8,
                                         colour="#2A2A2A",
                                         x=forecast.time.t+dt*diff(traj.times)[[1L]],
                                         ## y=pancast.quantile.display.trajs%>>%`[[`(forecast.time.i+dt,1L)
                                         ## y=-0.11*ymax+pancast.quantile.display.trajs%>>%`[[`(forecast.time.i+dt,1L)
                                  y=-0.13*ymax+pancast.quantile.display.trajs%>>%`[[`(forecast.time.i+dt,1L)
                                         )) %>>%
            `+`(ggplot2::annotate("segment",
                                  xend=forecast.time.t+dt*diff(traj.times)[[1L]],
                                  x=forecast.time.t+dt*diff(traj.times)[[1L]],
                                  yend=pancast.quantile.display.trajs%>>%`[[`(forecast.time.i+dt,1L),
                                  ## y=-0.11*ymax+pancast.quantile.display.trajs%>>%`[[`(forecast.time.i+dt,1L),
                                  y=-0.13*ymax+pancast.quantile.display.trajs%>>%`[[`(forecast.time.i+dt,1L),
                                  arrow=arrow(length=unit(0.02, "npc")),
                                  )) %>>%
            {.}
    } %>>%
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
            if (
                ## TRUE
                month.label.locations[[row.i,"month"]]!=5L
            ) {
                ## plt <- plt + ggplot2::annotate("text",
                plt <- plt + ggplot2::annotation_custom(
                                          grob=grid::textGrob(
                                                         hjust=0.5,
                                                         ## vjust=1.4,
                                                         vjust=2,
                                                         label=month.abb[[month.label.locations[[row.i,"month"]]]],
                                                         ),
                                          xmin=month.label.locations[[row.i,"rweek"]]-10,
                                          xmax=month.label.locations[[row.i,"rweek"]]+10,
                                          ymin=0-10,
                                          ymax=0+10
                                      )
            }
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
    ## `+`(ggplot2::annotate("text",
    ##                       ## x=xmin+0.2*diff(traj.times)[[1L]], y=fill.ymax,
    ##                       x=xmin+1.0*diff(traj.times)[[1L]], y=fill.ymax,
    ##                       ## nudge_x=grid::unit(0.02,"npc"), nudge_y=grid::unit(-0.02,"npc"),
    ##                       ## label=sprintf("FluSurv-NET %s hospitalization rates:\nearly revisions are often large and upward;\nlater re-revisions are usually smaller."),
    ##                       ## label=sprintf("FluSurv-NET overall hospitalization rates:\nfirst revisions are often large and upward;\nlater re-revisions are usually smaller."),
    ##                       ## label=sprintf("FluSurv-NET overall hospitalization rates:\noften revised significantly upward\nmultiple times in a row, then\nre-revised either upward or downward."),
    ##                       ## label=sprintf("FluSurv-NET overall hospitalization rates:\noften revised significantly upward\nmultiple times in a row, then\nre-revised less strongly\nin either direction."),
    ##                       ## label=sprintf("FluSurv-NET overall hospitalization rates:\nrevisions and re-revisions are generally upward;\nlater re-revisions are ",
    ##                       ## label=sprintf("%s hospitalization rates as of issue %d:\nrecent revisions (and re-revisions) are mostly upward;\nbackcasts predict mostly upward future revisions;\nforecast shows large chance of\nhigh upcoming rates.",
    ##                                     label=sprintf("%s hospitalization rates as of issue %d:\nrecent issues contain mostly upward revisions;\nbackcasts predict mostly upward future revisions;\nforecast shows large chance of\nhigh upcoming rates.",
    ##                                     voxel.data[["epigroup"]],
    ##                                     ## season_to_Season(voxel.data[["season"]],31L),
    ##                                     voxel.data[["issue"]]
    ##                                     ),
    ##                       hjust=0,
    ##                       ## vjust=1.05
    ##                       vjust=1.1
    ##                       )) %>>%
    `+`(ggplot2::annotate("text",
                          x=traj.times[[forecast.time.i-6L]]-1.0*diff(traj.times)[[1L]],
                          y=finalized.display.traj[[forecast.time.i-6L]],
                          ## hjust=1.05,
                          ## vjust=-0.1,
                          hjust=1,
                          label="Recent revisions\nand re-revisions\nare large, upward.",
                          lineheight=1.1
                          )) %>>%
    `+`(ggplot2::coord_cartesian(clip="off")) %>>%
    `+`(ggplot2::expand_limits(y=c(data.ymin,ymax))) %>>%
    `+`(ggplot2::expand_limits(x=c(xmin,xmax))) %>>%
    ## `+`(ggplot2::theme(legend.position="bottom")) %>>%
    `+`(ggplot2::theme(legend.position="right")) %>>%
    ## `+`(ggplot2::guides(fill=FALSE)) %>>%
    ## `+`(ggplot2::guides(fill=ggplot2::guide_legend(override.aes=list(colour="black")))) %>>%
    `+`(ggplot2::guides(colour=ggplot2::guide_legend(direction="vertical", reverse=TRUE, order=1L,
                                                     title="FluSurv-NET data",
                                                     ## labels=levels(traj.tbl[["Color"]]),
                                                     override.aes=list(shape=c(16L,rep(26L,n.older.lags+1L)))),
                        ## fill=ggplot2::guide_legend(direction="vertical", reverse=TRUE, order=2L, title="Pancasts"))) %>>%
                        fill=ggplot2::guide_legend(direction="vertical", reverse=TRUE, order=2L, title="Pancasts (percentiles for finalized data)"))) %>>%
    `+`(ggplot2::scale_x_continuous(
                 labels=function(breaks) model_week_to_epi_week(breaks,31L,lastWeekNumber(voxel.data[["season"]], 3L)),
                 breaks=round(month.sep.locations[["rweek"]]),
                 minor_breaks=round(month.sep.locations[["rweek"]])
    )) %>>%
    `+`(ggplot2::scale_y_continuous(
                     breaks=scales::pretty_breaks()(c(data.ymin, fill.ymax))
                 )) %>>%
    `+`(ggplot2::xlab("Epiweek")) %>>%
    `+`(ggplot2::ylab("Rate")) %>>%
    {.}
    ## `+`(ggplot2::ggsave("pancast_banded.png", width=width.in, height=height.in, units="in", dpi=300))
    ## `+`(ggplot2::ggsave("pancast_banded.png", width=width.in, height=height.in, units="in", dpi=300, type="cairo"))
## print(plt)
ggplot2::ggsave("pancast_banded.png", plt, width=width.in, height=height.in, units="in", dpi=300, type="cairo")
ggplot2::ggsave("pancast_banded.pdf", plt, width=width.in, height=height.in, units="in", dpi=300)


## grid::grid.newpage()
## print(plt)
## pBrackets::grid.brackets(xmin, fill.ymax+0.11*ymax, xmax, fill.ymax+0.13*ymax)
## pBrackets::grid.brackets(grid::unit(44, "native"), grid::unit(3, "native"), grid::unit(46,"native"), grid::unit(6, "native"), lwd=2)

## xxx try lemon package?

plt.gtbl = ggplot2::ggplotGrob(plt)
guide.box.gtbl = plt.gtbl$grobs[[which(plt.gtbl$layout$name=="guide-box")]]
fill.guide.gtbl = guide.box.gtbl$grobs[[2L]]

## label.gtree = fill.guide.gtbl$grobs[[which(fill.guide.gtbl$layout$name=="label-3-3")]]
## label.title.grob = label.gtree$children[[1L]]
## label.text.grob = label.title.grob$children[[1L]]
## label.text.grob$label

ribbon.guide.region.tbl =
    tibble::tibble(
                Region = sprintf("%04d",seq_len(n.regions)),
                Low = seq_len(n.regions)-1L,
                High = seq_len(n.regions)
            ) %>>%
    {
        dplyr::bind_rows(
                   . %>>% dplyr::mutate(Time=-1),
                   . %>>% dplyr::mutate(Time=1)
               )
    }
ribbon.guide.quantile.tbl =
    ## tibble::tibble(Percentile = 100*probs) %>>%
    tibble::tibble(Percentile = sprintf(
                       "Estimated %s percentile",
                       scales::ordinal(100/n.regions*(seq_len(n.regions+1L)-1L))
                   ),
                   Rate = seq_len(n.regions+1L)-1L
                   ) %>>%
    {
        dplyr::bind_rows(
                   . %>>% dplyr::mutate(Time=-1),
                   . %>>% dplyr::mutate(Time=1)
               )
    }
ribbon.guide.gtbl =
    ggplot2::ggplot(,ggplot2::aes(x=Time)) %>>%
    `+`(ggplot2::geom_ribbon(ggplot2::aes(ymin=Low,ymax=High,
                                          ## alpha=alpha,
                                          fill=Region),
                             ribbon.guide.region.tbl,
                             alpha=0.5
                             )) %>>%
    `+`(ggplot2::geom_line(ggplot2::aes(y=Rate, group=Percentile),
                           ribbon.guide.quantile.tbl,
                           colour="#2A2A2A"
                           )) %>>%
    `+`(ggplot2::annotate("text",
                          x=0,
                          y=seq_len(n.regions)-0.5,
                          ## label=sprintf("Estimated 1 in %d chance that\nfinalized point falls here",n.regions)
                          label=sprintf("Est. 1 in %d chance finalized is here",n.regions)
                          )) %>>%
    ggplot2::ggplotGrob() %>>%
    {.}

fill.guide.gtbl.p =
    fill.guide.gtbl %>>%
    gtable::gtable_filter("(?:background)|(?:title)") %>>%
    gtable::gtable_add_grob(ribbon.guide.gtbl$grobs[[which(ribbon.guide.gtbl$layout$name=="panel")]],
                            t=
                                fill.guide.gtbl$layout %>>%
                                dplyr::filter(name=="title") %>>%
                                {.[["b"]]+1L},
                            l=1,
                            b=
                                fill.guide.gtbl$layout %>>%
                                dplyr::filter(name=="background") %>>%
                                {.[["b"]]-1L}
                            ,
                            r=2) %>>%
    {.}
fill.guide.gtbl.p$widths[[2L]] <- fill.guide.gtbl.p$widths[[4L]]
## xxx maybe do the other way around --- add the title & background onto the plot?
guide.box.gtbl.p = guide.box.gtbl
guide.box.gtbl.p$grobs[[2L]] <- fill.guide.gtbl.p
plt.gtbl.p = plt.gtbl
plt.gtbl.p$grobs[[which(plt.gtbl.p$layout$name=="guide-box")]] <- guide.box.gtbl.p

## plt.gtbl$grobs[[which(plt.gtbl$layout$name=="panel")]] <-

grid::grid.newpage()
grid::grid.draw(plt.gtbl.p)


## todo legend showing quantile lines and fill region interpretation?

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




## xxx lineheight aes instead of manually placing multiple lines due to large \n

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
