
gc()
s.prospective.seasons = current.issue.sw[["season"]] %>>%
  stats::setNames(paste0(.,"/",.+1L)) %>>%
  with_dimnamesnames("Season")
w.prospective.model.weeks = current.issue.sw[["model.week"]] %>>%
  stats::setNames(paste0("MW",.)) %>>%
  with_dimnamesnames("Model Week")

print("Current season: select available data")
swg.prospective.voxel.data = map_join(
  get_voxel_data,
  s.prospective.seasons, w.prospective.model.weeks, g.epigroups,
  last.losocv.issue)

## Version of the above grouped by season and model week (an array of arrays of objects):
sw.g.prospective.voxel.data = map_join(
  function(swg.array, s,w) {
    swg.array[s,w,,drop=FALSE] %>>%
      select_dims(1:2, "drop") # drop s & w dimensions, but keep g dimension even if size 1
  },
  no_join(swg.prospective.voxel.data),
  named_arrayvec_to_name_arrayvec(s.prospective.seasons),
  named_arrayvec_to_name_arrayvec(w.prospective.model.weeks),
  lapply_variant=lapply, show.progress=FALSE
)

print("Current season: generate backcasts")
swgb.prospective.full.dats = map_join(
  get_backcast,
  swg.prospective.voxel.data, sw.g.prospective.voxel.data, source.name, signal.name, b.backcasters,
  epidata_df_to_chopped_trajectory_df=epidata_df_to_chopped_trajectory_df
  , shuffle=FALSE
  , cache.prefix=file.path(epiproject.cache.dir,"swgb.prospective.full.dats")
)

print("Current season: generate component forecasts")
swgbf.prospective.component.target.multicasts = map_join(
  target_multicast,
  swg.prospective.voxel.data, swgb.prospective.full.dats, f.forecasters,
  target_trajectory_preprocessor,
  no_join(t.target.specs),
  no_join(m.forecast.types),
  full_dat_fixup=full_dat_fixup
  , cache.prefix=file.path(epiproject.cache.dir,"swgbf.prospective.component.target.multicasts")
)

swgtmbf.prospective.component.forecast.values =
  swgbf.prospective.component.target.multicasts %>>%
  ## first, get forecast.value's in swgbf.tm format:
  map_join(f=`[[`, "forecast.values") %>>%
  ## un-nest lists to get swgbftm format:
  map_join(f=`[[`,
           named_arrayvec_to_name_arrayvec(t.target.specs),
           named_arrayvec_to_name_arrayvec(m.forecast.types)
           ) %>>%
  ## permute dimension order to get desired swgtmbf format:
  aperm(c(1:3,6:7,4:5))
