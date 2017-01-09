library(devtools)
outputdir = "."
path_to_Rpkg_directory = "."
setwd(path_to_Rpkg_directory)
load_all()
build()
document()
check()

