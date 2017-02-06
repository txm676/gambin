#' @rdname dgambin
#' @export
gambin_exp = function(alpha, maxoctave, w = 1, total_species)
{
  dgambin(0:maxoctave, alpha, w, maxoctave) * total_species
}
