
backfill_ignorant_backsim = function(voxel.data, signal.name) {
  old.dat = voxel.data[["epidata.df"]] %>>%
    dplyr::filter(season != voxel.data[["season"]]) %>>%
    dplyr::group_by(season) %>>%
    {.[!any(is.na(.[[signal.name]])),]} %>>%
    dplyr::ungroup() %>>%
    {split(.[[signal.name]], .[["Season"]])}
  new.dat = voxel.data[["epidata.df"]] %>>%
    dplyr::filter(season == voxel.data[["season"]]) %>>%
    magrittr::extract2(signal.name)
  new.dat.sim = match.new.dat.sim(new.dat)
  voxel.Season = season_to_Season(
    voxel.data[["season"]], voxel.data[["first.week.of.season"]])
  ## concatenate new.dat.sim onto old.dat, setting the :
  full.dat = c(old.dat,
               setNames(list(new.dat.sim), voxel.Season))
  return (full.dat)
}