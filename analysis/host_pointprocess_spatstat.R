# =============================================================================
# Host-conditioned point-process analysis (spatstat). Consumes per-block CSV grids from
# analysis/export_blocks_for_R.py (rasterio bridge) -> only spatstat needed (no terra/GDAL).
#
# Intensity model is fit GLOBALLY (pooled forest cells across all blocks) so it captures the
# regional host gradient (host-tracking is largely between-block); the fitted lambda is then
# applied PER BLOCK for the inhomogeneous L-function with host-conditioned simulation envelopes.
# Pre-mortality (2021) MS-NFI spruce + total volume + stand-edge.  Primary object: L_inhom.
# =============================================================================
suppressPackageStartupMessages(library(spatstat))
RES <- 100; RMAX <- 1500; NSIM <- 39; RDIR <- "analysis/out/rblocks"

grid_im <- function(df, z) {                       # regular grid -> im
  xs <- sort(unique(df$x)); ys <- sort(unique(df$y))
  m <- matrix(NA_real_, length(ys), length(xs))
  m[cbind(match(df$y, ys), match(df$x, xs))] <- df[[z]]
  im(m, xcol = xs, yrow = ys)
}

load_block <- function(b) {
  g  <- read.csv(sprintf("%s/block%d_grid.csv",  RDIR, b))
  pt <- read.csv(sprintf("%s/block%d_pts.csv",   RDIR, b))
  ti <- read.csv(sprintf("%s/block%d_tiles.csv", RDIR, b))
  spr <- grid_im(g, "spr"); tot <- grid_im(g, "tot"); lc <- grid_im(g, "lc")
  forest <- solutionset(eval.im(round(lc) == 1))
  tw <- lapply(seq_len(nrow(ti)), function(i) owin(c(ti$x0[i], ti$x1[i]), c(ti$y0[i], ti$y1[i])))
  W  <- intersect.owin(do.call(union.owin, tw), forest)
  edge <- bdist.pixels(forest)
  # forest cells (centres) + covariates + cluster counts (for the global GLM)
  fcell <- g[round(g$lc) == 1, ]
  fcell$edge <- lookup.im(edge, fcell$x, fcell$y, naok = TRUE, strict = FALSE)
  fcell <- fcell[is.finite(fcell$edge), ]
  xs <- sort(unique(g$x)); ys <- sort(unique(g$y))
  ix <- round((pt$cx - min(xs)) / RES) + 1; iy <- round((pt$cy - min(ys)) / RES) + 1
  key <- paste(pmax(1, pmin(length(xs), ix)), pmax(1, pmin(length(ys), iy)))
  ckey <- paste(match(fcell$x, xs), match(fcell$y, ys))
  fcell$count <- as.integer(table(factor(key, levels = ckey)))
  list(b = b, spr = spr, tot = tot, edge = edge, forest = forest, W = W,
       pp = ppp(pt$cx, pt$cy, window = W, checkdup = FALSE),
       cells = data.frame(b = b, x = fcell$x, y = fcell$y, spr = fcell$spr,
                          tot = fcell$tot, edge = fcell$edge, count = fcell$count))
}

main <- function() {
  blocks <- sort(as.integer(gsub("\\D", "", list.files(RDIR, "_grid.csv$"))))
  L <- lapply(blocks, load_block); names(L) <- blocks
  # ---- GLOBAL intensity model (pooled forest cells) ----
  allc <- do.call(rbind, lapply(L, `[[`, "cells"))
  mu <- colMeans(allc[c("spr","tot","edge")]); sdv <- sapply(allc[c("spr","tot","edge")], sd)
  zc <- scale(allc[c("spr","tot","edge")])
  glm_g <- glm(allc$count ~ zc, family = poisson())
  co <- coef(glm_g)
  cat(sprintf("GLOBAL intensity  lambda ~ spruce+total+edge  (standardised): spr=%+.3f tot=%+.3f edge=%+.3f\n",
              co[2], co[3], co[4]))
  cell_area <- RES * RES
  predict_lambda <- function(cells) {              # fitted intensity per unit area at given cells
    z <- scale(cells[c("spr","tot","edge")], center = mu, scale = sdv)
    exp(cbind(1, z) %*% co) / cell_area
  }
  # ---- per-block L_inhom with the GLOBAL lambda + host-conditioned envelopes ----
  res <- list()
  cat(sprintf("%3s %6s | %s\n","blk","n","L_inhom(r)-r vs GLOBAL-host envelope  250/500/1000m"))
  for (b in blocks) {
    lb <- L[[as.character(b)]]
    lam_cells <- predict_lambda(lb$cells)
    lam <- grid_im(data.frame(x = lb$cells$x, y = lb$cells$y, z = lam_cells), "z")
    lam <- lam[lb$W, drop = FALSE]                 # restrict to the forest window
    if (npoints(lb$pp) < 60) next
    E <- envelope(lb$pp, Linhom, funargs = list(lambda = lam),
                  simulate = expression(rpoispp(lam)),
                  nsim = NSIM, correction = "translate", rmax = RMAX, verbose = FALSE)
    at <- function(r0){ i <- which.min(abs(E$r-r0)); c(obs=E$obs[i]-r0, lo=E$lo[i]-r0, hi=E$hi[i]-r0) }
    f <- sapply(c(250,500,1000), at); res[[as.character(b)]] <- list(b=b, n=npoints(lb$pp), E=E, f=f)
    v <- function(j) if (f["obs",j]>f["hi",j]) "ABOVE" else if (f["obs",j]<f["lo",j]) "below" else "within"
    cat(sprintf("%3d %6d | %s/%s/%s   (L500-r=%.0f, env=[%.0f,%.0f])\n",
        b, npoints(lb$pp), v(1), v(2), v(3), f["obs",2], f["lo",2], f["hi",2]))
  }
  saveRDS(list(coef=co, res=res), "analysis/out/spatstat_host.rds")
  if (length(res)) {
    pdf("figures/host_conditioned_spatstat.pdf", width=10, height=7)
    par(mfrow=c(ceiling(length(res)/3),3), mar=c(4,4,2,1))
    for (nm in names(res)) plot(res[[nm]]$E, .-r~r, main=paste("block",nm), legend=FALSE)
    dev.off()
  }
  cat("\nsaved analysis/out/spatstat_host.rds + figures/host_conditioned_spatstat.pdf\n")
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
if (length(a) && a[1]=="selftest") selftest() else main()
