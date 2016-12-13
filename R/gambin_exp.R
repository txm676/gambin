#' @rdname dgambin
#' @export
gambin_exp = function(alpha, maxoctave, total_species)
{
  dgambin(0:maxoctave, alpha, maxoctave) * total_species
}
