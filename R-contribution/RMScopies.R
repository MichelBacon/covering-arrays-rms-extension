#' Strength 3 CAs from the Raaphorst-Moura-Stevens construction, extended by shifted column copies
#'
#' The Raaphorst-Moura-Stevens (2014) construction yields a CA(2q^3-1; 3, q^2+q+1, q)
#' for prime \code{q}. Function \code{RMScopies} extends it by \code{m} additional columns,
#' each of which is a shifted copy of an existing column, plus a few layers of
#' additional runs.
#'
#' @rdname RMScopies
#'
#' @aliases RMScopiesCA
#' @aliases RMScopies
#' @aliases N_RMScopiesCA
#'
#' @usage RMScopiesCA(k, v, ...)
#' @usage RMScopies(q, m=0)
#' @usage N_RMScopiesCA(t=3, k, v, ...)
#'
#' @param k integer: the target number of columns
#' @param v integer: the number of levels (must be a prime)
#' @param q integer: a prime
#' @param m integer: the number of added columns (shifted copies),
#'          \code{0 <= m <= q^2+q+1}
#' @param t integer: the target strength (always 3)
#' @param ... currently not used
#'
#' @details
#' \bold{Base array} (\code{m=0}, Raaphorst, Moura and Stevens 2014):
#' the k0 = q^2+q+1 columns correspond to the points \eqn{p_j}{p_j} of the
#' projective plane PG(2,q), ordered as the powers of a root of a primitive
#' cubic polynomial over GF(q) (Singer cycle).
#' The array consists of a row of zeroes, a block B1 with entries
#' \eqn{x \cdot p_j}{x.p_j} for all non-zero x in GF(q)^3, and a block B2
#' constructed in the same way with the column order reversed.
#' Three columns that are collinear in B1 are not collinear in B2, so that
#' all triples of columns are covered.
#'
#' \bold{Shifted copies} (\code{m>=1}): new column \eqn{n_i}{n_i} is a copy of
#' column \eqn{c_i = i}{c_i = i} in the zero row and in B1, and equals
#' \eqn{c_i + 1}{c_i + 1} (mod q) in B2. The only uncovered 3-tuples are then
#' found in column triples \eqn{\{c_i, n_i, w\}}{{c_i, n_i, w}}, for differences
#' \eqn{n_i - c_i \in \{2, \dots, q-1\}}{n_i - c_i in {2, ..., q-1}}.
#'
#' \bold{Layers} (\code{m>=2}): each layer consists of q^2 runs indexed by
#' y in GF(q)^2. The q+1 directions of PG(1,q) yield q+1 mutually orthogonal
#' linear functions \eqn{y \cdot P}{y.P}. In a layer, up to q copies are active:
#' active copy i receives its own direction (\eqn{c_i = y \cdot P_s}{c_i = y.P_s},
#' \eqn{n_i = y \cdot P_s + e}{n_i = y.P_s + e}),
#' all other columns receive direction \eqn{P_0}{P_0}. Each copy is active in
#' q-2 layers, with e = 2, ..., q-1. A final run (0 for the original columns,
#' 1 for the new columns) covers the image of x=0 in B2.
#' The number of layers is \eqn{L = \max(q-2, \lceil m(q-2)/q \rceil)}{L = max(q-2, ceiling(m(q-2)/q))}.
#'
#' For \code{m=1}, the q(q-2) layers' worth of runs are replaced by q^2(q-2) simple runs
#' (column \eqn{c_1=a}{c_1=a}, column \eqn{n_1 = a+e}{n_1 = a+e}, all other columns equal to s),
#' and no final run is needed.
#'
#' Run sizes: \eqn{2q^3-1}{2q^3-1} for \code{m=0},
#' \eqn{2q^3-1+q^2(q-2)}{2q^3-1+q^2(q-2)} for \code{m=1},
#' \eqn{2q^3+q^2 L}{2q^3+q^2 L} for \code{m>=2}.
#'
#' Compared to the Colbourn tables (status November 2024), this improves
#' 60 entries: v=5 (k=32,...,37), v=7 (k=58,...,65), v=11 (k=134,...,145),
#' v=13 (k=184,...,197), v=19 (k=382,...,401); for example CA(324; 3, 32, 5)
#' instead of 349 runs, and CA(350; 3, 37, 5) instead of 365 runs.
#'
#' @returns Function \code{RMScopiesCA} returns a matrix of class \code{ca}
#' with attributes, which is a strength 3 covering array.\cr
#' Function \code{RMScopies} returns a matrix without class and attributes,\cr
#' function \code{N_RMScopiesCA} returns the run size of the CA
#' (\code{NA} if the construction is not applicable).
#'
#' @references Raaphorst, S., Moura, L. and Stevens, B. (2014). A construction for
#' strength-3 covering arrays from linear feedback shift register sequences.
#' Designs, Codes and Cryptography 73, 949-968.
#'
#' @examples
#' D <- RMScopies(5)        ## the RMS array CA(249; 3, 31, 5)
#' dim(D)
#' E <- RMScopiesCA(32, 5)  ## CA(324; 3, 32, 5)
#' dim(E)
#' N_RMScopiesCA(3, 36, 5)  ## 325
#' N_RMScopiesCA(3, 37, 5)  ## 350
#' \dontrun{
#' caverify::ca_verify(E, 3)$covered
#' }
#'
#' @export
RMScopies <- function(q, m=0){
  stopifnot(length(q) == 1, q == round(q), q >= 2)
  if (!.RMS_isprime(q)) stop("q must be a prime")
  k0 <- q^2 + q + 1
  stopifnot(m == round(m), m >= 0, m <= k0)
  ## Singer cycle: powers of a root of a primitive cubic
  P <- .RMS_singer_points(q)[1:k0, , drop=FALSE]
  X <- as.matrix(expand.grid(0:(q-1), 0:(q-1), 0:(q-1)))[-1, , drop=FALSE]  ## non-zero vectors
  dimnames(X) <- NULL
  B1 <- (X %*% t(P)) %% q
  B2 <- (X %*% t(P[k0:1, , drop=FALSE])) %% q
  aus <- rbind(rep(0, k0), B1, B2)
  if (m == 0) return(aus)
  ## shifted copies of columns 1, ..., m
  inB2 <- (1 + nrow(B1) + 1):nrow(aus)
  neu <- aus[, 1:m, drop=FALSE]
  neu[inB2, ] <- (neu[inB2, ] + 1) %% q
  aus <- cbind(aus, neu)
  k <- k0 + m
  if (m == 1){
    extra <- do.call(rbind, lapply(0:(q-1), function(a)
      do.call(rbind, lapply(2:(q-1), function(e)
        t(sapply(0:(q-1), function(s){
          r <- rep(s, k); r[1] <- a; r[k] <- (a + e) %% q; r
        }))))))
    return(rbind(aus, extra))
  }
  ## layers
  act <- cbind(i=rep(1:m, times=q-2), e=rep(2:(q-1), each=m))
  if (m <= q) layers <- lapply(2:(q-1), function(e) cbind(i=1:m, e=e))
  else layers <- lapply(seq(1, nrow(act), by=q),
                        function(st) act[st:min(st + q - 1, nrow(act)), , drop=FALSE])
  dirs <- rbind(c(1, 0), c(0, 1), cbind(1, 1:(q-1)))   ## P_0, ..., P_q
  Y <- as.matrix(expand.grid(0:(q-1), 0:(q-1))); dimnames(Y) <- NULL
  extra <- do.call(rbind, lapply(layers, function(lay){
    stopifnot(!any(duplicated(lay[, "i"])))
    B <- matrix((Y %*% dirs[1, ]) %% q, nrow(Y), k)
    for (s in seq_len(nrow(lay))){
      f <- (Y %*% dirs[s + 1, ]) %% q
      B[, lay[s, "i"]] <- f
      B[, k0 + lay[s, "i"]] <- (f + lay[s, "e"]) %% q
    }
    B
  }))
  rbind(aus, extra, c(rep(0, k0), rep(1, m)))
}

#' @export
RMScopiesCA <- function(k, v, ...){
  ## t=3 always holds
  Call <- sys.call()
  if (!.RMS_isprime(v)) stop("v must be a prime")
  k0 <- v^2 + v + 1
  m <- max(0, k - k0)
  if (m > k0) stop("k must not exceed 2*(v^2+v+1)")
  aus <- RMScopies(v, m)[, 1:k, drop=FALSE]
  class(aus) <- c("ca", class(aus))
  attr(aus, "origin") <- if (m == 0) "Raaphorst-Moura-Stevens 2014"
                         else "Raaphorst-Moura-Stevens 2014, with shifted column copies"
  attr(aus, "Call") <- Call
  attr(aus, "t") <- 3
  aus
}

#' @export
N_RMScopiesCA <- function(t=3, k, v, ...){
  if (!t == 3) return(NA)
  if (!.RMS_isprime(v)) return(NA)
  q <- v; k0 <- q^2 + q + 1
  m <- max(0, k - k0)
  if (m > k0) return(NA)
  if (m == 0) return(2*q^3 - 1)
  if (m == 1) return(2*q^3 - 1 + q^2*(q - 2))
  2*q^3 + q^2*max(q - 2, ceiling(m*(q - 2)/q))
}

## internal helpers
.RMS_isprime <- function(n){
  if (length(n) != 1 || n != round(n) || n < 2) return(FALSE)
  if (n < 4) return(TRUE)
  all(n %% 2:floor(sqrt(n)) != 0)
}

.RMS_singer_points <- function(q){
  ## returns the q^3-1 powers of alpha as rows (coefficients of 1, x, x^2),
  ## alpha a root of a primitive cubic x^3 + f2 x^2 + f1 x + f0 over GF(q)
  mulx <- function(a, f){
    ## a: coefficients (1, x, x^2); multiply by x and reduce mod f (monic, degree 3)
    top <- a[3]
    (c(0, a[1], a[2]) - top*f) %% q
  }
  for (f0 in 1:(q-1)) for (f1 in 0:(q-1)) for (f2 in 0:(q-1)){
    f <- c(f0, f1, f2)
    cur <- c(1, 0, 0)
    pw <- matrix(0, q^3 - 1, 3)
    for (j in 1:(q^3 - 1)){
      pw[j, ] <- cur
      cur <- mulx(cur, f)
    }
    if (all(cur == c(1, 0, 0)) &&
        !any(duplicated(pw[, 1] + q*pw[, 2] + q^2*pw[, 3])))
      return(pw)
  }
  stop("no primitive polynomial found")
}
