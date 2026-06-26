# =============================================================================
# Host-conditioned point-process analysis (spatstat). Consumes per-block CSV grids from
# analysis/export_blocks_for_R.py (rasterio bridge) -> only spatstat needed (no terra/GDAL).
#
# Intensity model is fit GLOBALLY (pooled forest cells across all blocks) so it captures the
# regional host gradient (host-tracking is largely between-block); the fitted lambda is then
# applied PER BLOCK for the inhomogeneous L-function with host-conditioned simulation envelopes.
# Pre-mortality (2021) MS-NFI spruce + total volume + stand-edge.  Primary object: L_inhom.
#
# Usage:  Rscript host_pointprocess_spatstat.R [centroids|trees] [RMAX]   (default: centroids 1500)
#         Rscript host_pointprocess_spatstat.R selftest
#   unit  = analysis object: cluster centroids (n~370-4275/blk) or individual trees (subsampled, 20k cap)
#   RMAX  = max radius (m); block windows are 25-54 km so reliable r ~ window/4 (>= 3 km supported)
# =============================================================================
suppressPackageStartupMessages(library(spatstat))
RES <- 100; NSIM <- 39; RDIR <- "analysis/out/rblocks"
REPORT_R <- c(250, 500, 1000, 1500, 2000, 3000)

grid_im <- function(df, z) {                       # regular grid -> im
  xs <- sort(unique(df$x)); ys <- sort(unique(df$y))
  m <- matrix(NA_real_, length(ys), length(xs))
  m[cbind(match(df$y, ys), match(df$x, xs))] <- df[[z]]
  im(m, xcol = xs, yrow = ys)
}

load_block <- function(b, unit) {
  g  <- read.csv(sprintf("%s/block%d_grid.csv",  RDIR, b))
  pf <- if (unit == "trees") "trees" else "pts"
  pt <- read.csv(sprintf("%s/block%d_%s.csv",    RDIR, b, pf))
  ti <- read.csv(sprintf("%s/block%d_tiles.csv", RDIR, b))
  spr <- grid_im(g, "spr"); tot <- grid_im(g, "tot"); lc <- grid_im(g, "lc")
  forest <- solutionset(eval.im(round(lc) == 1))
  tw <- lapply(seq_len(nrow(ti)), function(i) owin(c(ti$x0[i], ti$x1[i]), c(ti$y0[i], ti$y1[i])))
  W  <- intersect.owin(do.call(union.owin, tw), forest)
  edge <- bdist.pixels(forest)
  fcell <- g[round(g$lc) == 1, ]
  fcell$edge <- lookup.im(edge, fcell$x, fcell$y, naok = TRUE, strict = FALSE)
  fcell <- fcell[is.finite(fcell$edge), ]
  if (unit == "trees") {
    fcell$count <- fcell$ntree                      # intensity model fit on FULL tree counts/cell
  } else {
    xs <- sort(unique(g$x)); ys <- sort(unique(g$y))
    ix <- round((pt$cx - min(xs)) / RES) + 1; iy <- round((pt$cy - min(ys)) / RES) + 1
    key <- paste(pmax(1, pmin(length(xs), ix)), pmax(1, pmin(length(ys), iy)))
    ckey <- paste(match(fcell$x, xs), match(fcell$y, ys))
    fcell$count <- as.integer(table(factor(key, levels = ckey)))
  }
  list(b = b, spr = spr, tot = tot, edge = edge, forest = forest, W = W,
       pp = ppp(pt$cx, pt$cy, window = W, checkdup = FALSE),
       cells = data.frame(b = b, x = fcell$x, y = fcell$y, spr = fcell$spr,
                          tot = fcell$tot, edge = fcell$edge, count = fcell$count))
}

main <- function(unit = "centroids", rmax = 1500) {
  radii <- REPORT_R[REPORT_R <= rmax]
  cat(sprintf("=== unit=%s  RMAX=%dm  radii=[%s] ===\n", unit, rmax, paste(radii, collapse=",")))
  blocks <- sort(as.integer(gsub("\\D", "", list.files(RDIR, "_grid.csv$"))))
  L <- lapply(blocks, load_block, unit = unit); names(L) <- blocks
  # ---- GLOBAL intensity model (pooled forest cells) ----
  allc <- do.call(rbind, lapply(L, `[[`, "cells"))
  mu <- colMeans(allc[c("spr","tot","edge")]); sdv <- sapply(allc[c("spr","tot","edge")], sd)
  zc <- scale(allc[c("spr","tot","edge")])
  glm_g <- glm(allc$count ~ zc, family = poisson())
  co <- coef(glm_g)
  cat(sprintf("GLOBAL intensity  lambda ~ spruce+total+edge  (standardised): spr=%+.3f tot=%+.3f edge=%+.3f\n",
              co[2], co[3], co[4]))
  cell_area <- RES * RES
  predict_lambda <- function(cells) {
    z <- scale(cells[c("spr","tot","edge")], center = mu, scale = sdv)
    exp(cbind(1, z) %*% co) / cell_area
  }
  res <- list()
  cat(sprintf("%3s %7s | L_inhom(r)-r vs GLOBAL-host envelope @ [%s] m\n","blk","n",paste(radii,collapse="/")))
  for (b in blocks) {
    lb <- L[[as.character(b)]]
    if (npoints(lb$pp) < 60) next
    lam_cells <- predict_lambda(lb$cells)
    lam <- grid_im(data.frame(x = lb$cells$x, y = lb$cells$y, z = lam_cells), "z")
    lam <- lam[lb$W, drop = FALSE]
    sf <- npoints(lb$pp) / integral.im(lam)                  # rescale null to observed count (handles subsample)
    lam <- eval.im(lam * sf)
    E <- envelope(lb$pp, Linhom, funargs = list(lambda = lam),
                  simulate = expression(rpoispp(lam)),
                  nsim = NSIM, correction = "translate", rmax = rmax, verbose = FALSE)
    at <- function(r0){ i <- which.min(abs(E$r-r0)); c(obs=E$obs[i]-r0, lo=E$lo[i]-r0, hi=E$hi[i]-r0) }
    f <- sapply(radii, at); colnames(f) <- radii; res[[as.character(b)]] <- list(b=b, n=npoints(lb$pp), E=E, f=f)
    v <- function(j) if (f["obs",j]>f["hi",j]) "ABOVE" else if (f["obs",j]<f["lo",j]) "below" else "withn"
    cat(sprintf("%3d %7d | %s\n", b, npoints(lb$pp),
        paste(sapply(seq_along(radii), function(j) sprintf("%s(%.0f)", v(j), f["obs",j])), collapse=" ")))
  }
  tag <- if (unit == "trees") "_trees" else ""
  saveRDS(list(unit=unit, rmax=rmax, coef=co, res=res),
          sprintf("analysis/out/spatstat_host%s.rds", tag))
  if (length(res)) {
    pdf(sprintf("figures/host_conditioned_spatstat%s.pdf", tag), width=10, height=7)
    par(mfrow=c(ceiling(length(res)/3),3), mar=c(4,4,2,1))
    for (nm in names(res)) plot(res[[nm]]$E, .-r~r, main=paste(unit,"block",nm), legend=FALSE)
    dev.off()
  }
  cat(sprintf("\nsaved analysis/out/spatstat_host%s.rds + figures/host_conditioned_spatstat%s.pdf\n", tag, tag))
}

selftest <- function() {
  set.seed(1); W <- owin(c(0,6000),c(0,6000))
  lam <- as.im(function(x,y) 50e-6*exp(2*x/6000), W); pp <- rpoispp(lam)
  E <- envelope(pp, Linhom, funargs=list(lambda=lam), simulate=expression(rpoispp(lam)),
                nsim=39, correction="translate", rmax=800, verbose=FALSE)
  i <- which.min(abs(E$r-500))
  cat(sprintf("[selftest] n=%d Linhom(500)-500=%.1f env=[%.1f,%.1f] -> %s\n",
      npoints(pp), E$obs[i]-500, E$lo[i]-500, E$hi[i]-500,
      ifelse(E$obs[i]>E$lo[i]&&E$obs[i]<E$hi[i],"PASS","CHECK")))
}

a <- commandArgs(trailingOnly=TRUE)
if (length(a) && a[1]=="selftest") {
  selftest()
} else {
  unit <- if (length(a) >= 1 && nzchar(a[1])) a[1] else "centroids"
  rmax <- if (length(a) >= 2) as.numeric(a[2]) else 1500
  main(unit, rmax)
}
