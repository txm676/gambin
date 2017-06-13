#' @rdname logLik.gambin
#' @importFrom stats BIC
#' @export
BIC.gambin = function(object, ...)
{
  if (length(list(...)) > 0L) 
    warning("additional arguments ignored")
  stats::BIC(logLik(object))
}
