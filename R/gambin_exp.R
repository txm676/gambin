#' @rdname dgambin
#' @export
gambin_exp <-
function(alpha, maxoctave, total_species)
{
  exp <- dgambin(0:maxoctave, alpha, maxoctave)
  exp <- exp * total_species 
  return(exp)
}
