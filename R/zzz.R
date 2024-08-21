.onAttach <- function(libname, pkgname) {
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s", pkgname, "3.21")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
}
