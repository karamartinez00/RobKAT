# Nothing so far
.onAttach <- function(libname, pkgname) {
}

.onLoad <- function(libname, pkgname) {
}

# If we C++, this is good practice
.onUnload <- function (libpath) {
  library.dynam.unload("RobKAT", libpath)
}
