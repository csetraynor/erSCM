#' Generate reproducible per-replicate L'Ecuyer-CMRG RNG streams
#'
#' Derives \code{n} independent L'Ecuyer-CMRG random-number streams from a single
#' integer \code{seed}, so that a simulation split into \code{n} replicates
#' produces identical results whether the replicates are evaluated sequentially
#' or dispatched to parallel workers via
#' \code{future.apply::future_lapply(future.seed = streams)}.
#'
#' The caller's RNG kind and \code{.Random.seed} are captured up front; the
#' returned \code{restore} closure reinstates them exactly. Callers must register
#' this closure with \code{on.exit(<result>$restore(), add = TRUE)} immediately
#' after calling, so the internal switch to \code{"L'Ecuyer-CMRG"} does not leak
#' into the caller's session.
#'
#' @param n Integer; number of streams (replicates) to generate. Must be >= 1.
#' @param seed Integer; the master seed the streams are derived from.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{streams}}{A length-\code{n} list of \code{.Random.seed}
#'       vectors, one per replicate.}
#'     \item{\code{restore}}{A zero-argument function that restores the caller's
#'       original RNG kind and \code{.Random.seed}.}
#'   }
#'   As a side effect the global RNG kind is switched to \code{"L'Ecuyer-CMRG"}
#'   and remains active until \code{restore()} is invoked.
#'
#' @keywords internal
.make_rng_streams <- function(n, seed) {
  old_kind <- RNGkind()
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }

  restore <- function() {
    RNGkind(
      kind        = old_kind[1L],
      normal.kind = old_kind[2L],
      sample.kind = old_kind[3L]
    )
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  streams <- vector("list", n)
  streams[[1L]] <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (n > 1L) {
    for (i in 2:n) {
      streams[[i]] <- parallel::nextRNGStream(streams[[i - 1L]])
    }
  }

  list(streams = streams, restore = restore)
}
