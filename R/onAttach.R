.onAttach <- function(libname, pkgname) {
    packageStartupMessage("Attaching chiptsne version ",
                          packageDescription("chiptsne")$Version, ".")
    #When adding new options here, also add them to the "names" setMethod below
    CT_DIMREDUCE_METHODS <<- valid_dimreduce_methods
}
