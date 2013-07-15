.First.lib <- function(lib, pkg) {
   library.dynam("QTLRel", pkg, lib)
}
.onLoad <- function(lib, pkg) cat("R/QTLRel is loaded\n")
.noGenerics <- TRUE
.onUnload <- function(libpath) library.dynam.unload("QTLRel", libpath)

