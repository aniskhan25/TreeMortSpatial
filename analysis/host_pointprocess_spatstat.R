# =============================================================================
# Host-conditioned point-process analysis of boreal canopy mortality
# LOCKED DESIGN (Vasquez + Lindqvist). Run in R with spatstat + terra.
#
# NOTE: authored without an R environment to test against — validate interactively.
# The heavy steps (window union, tree-level envelopes) are flagged with scalability notes.
#
# Inputs (EPSG:3067 / ETRS-TM35FIN throughout) — MVMI 2021 "1721" product (PRE-mortality):
#   data/trees_clustered.csv          x, y, cluster_id, image_name   (698k dead-tree points)
#   data/kuusi_vmi1x_1721.tif         spruce volume      -> host covariate            [Lindqvist]
#   data/tilavuus_vmi1x_1721.tif      total volume       -> richer-model covariate
#   data/maaluokka_vmi1x_1721.tif     land class         -> forest mask (== forest land)
# (2023 spruce is ENDOGENOUS — do not use it as the host; see REVISION_PLAN.md.)
#
# Grab the rasters (rsync; ~3.5 GB total):
#   base=rsync://rsync.nic.funet.fi/ftp/index/geodata/luke/vmi/2021
#   rsync -P $base/kuusi_vmi1x_1721.tif $base/tilavuus_vmi1x_1721.tif $base/maaluokka_vmi1x_1721.tif data/
# =============================================================================
suppressPackageStartupMessages({ library(spatstat); library(terra); library(data.table) })

NODATA <- 32767
TILE_RES   <- 50     # m, resolution of covariate images / window mask (coarsen for speed)
RMAX       <- 1500   # m, max distance for K/L/pcf
NSIM       <- 99     # envelope simulations (drop to 39 while prototyping)

# ---- 1. points -------------------------------------------------------------
tr <- fread("data/trees_clustered.csv")[cluster_id != -1, .(x, y, cluster_id, image_name)]

# ---- 2. rasters -> spatstat im (host + forest mask) ------------------------
rast_to_im <- function(path, bbox, res = TILE_RES) {
  r <- rast(path)
  r <- crop(r, ext(bbox$xmin, bbox$xmax, bbox$ymin, bbox$ymax))
  r <- aggregate(r, fact = max(1, round(res / 16)), fun = "mean", na.rm = TRUE)
  r[r == NODATA] <- NA
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE); names(df) <- c("x", "y", "z")
  as.im(df)                              # spatstat im on a regular grid
}
bbox <- list(xmin = min(tr$x) - 3000, xmax = max(tr$x) + 3000,
             ymin = min(tr$y) - 3000, ymax = max(tr$y) + 3000)
spr <- rast_to_im("data/kuusi_vmi1x_1721.tif",    bbox)   # spruce (host) intensity covariate
tot <- rast_to_im("data/tilavuus_vmi1x_1721.tif", bbox)   # total volume (richer covariate)
mlk <- rast_to_im("data/maaluokka_vmi1x_1721.tif", bbox)  # land class (1 = metsämaa / forest land)

# ---- 3. observation window = forest-masked union of 6 km tile windows ------
# tile extent from tree bbox per image_name (~6 km), then keep only forest pixels.
tiles <- tr[, .(x0 = min(x), x1 = max(x), y0 = min(y), y1 = max(y), n = .N), by = image_name][n >= 10]
tile_win <- with(tiles, mapply(function(a,b,c,d) owin(c(a,b), c(c,d)),
                               x0, x1, y0, y1, SIMPLIFY = FALSE))
Wtiles <- do.call(union.owin, tile_win)                   # multi-rectangle survey footprint
forest <- solutionset(mlk == 1)                           # land-class forest mask (metsämaa); NOT spruce>0
W <- intersect.owin(Wtiles, forest)                       # << the correct observation window
# (Scalability: if W is too fine/large, set its mask resolution via as.mask(W, eps=TILE_RES).)

pp <- ppp(tr$x, tr$y, window = W, checkdup = FALSE)
cat(sprintf("points in window: %d ; window area: %.0f km^2\n", npoints(pp), area(W)/1e6))

# ---- 4. intensity models ---------------------------------------------------
# Primary host model: lambda(u) propto spruce(u)  (log-linear, Berman-Turner via ppm)
fit_spr  <- ppm(pp ~ spr, covariates = list(spr = spr))
# Richer discriminator: spruce + total volume + distance-to-forest-edge
edge <- distfun(as.psp(as.polygonal(W)))                  # crude edge covariate; refine as needed
fit_full <- ppm(pp ~ spr + tot + edge, covariates = list(spr = spr, tot = tot, edge = edge))
print(anova(fit_spr, fit_full, test = "LR"))              # does the richer model matter?

lam_spr  <- predict(fit_spr,  type = "trend")
lam_full <- predict(fit_full, type = "trend")

# ---- 5. PRIMARY inferential object: Linhom + host-conditioned envelopes -----
# Envelopes simulate from the FITTED INHOMOGENEOUS null (NOT CSR).
set.seed(42)
E_spr <- envelope(pp, Linhom, funargs = list(lambda = lam_spr),
                  simulate = expression(rpoispp(lam_spr)),
                  nsim = NSIM, correction = "translate", rmax = RMAX, savefuns = TRUE)
# Discriminator: does residual aggregation survive the richer covariate model?
E_full <- envelope(pp, Linhom, funargs = list(lambda = lam_full),
                   simulate = expression(rpoispp(lam_full)),
                   nsim = NSIM, correction = "translate", rmax = RMAX)
# Display PCF (Epanechnikov kernel; Stoyan bandwidth) -- scale localization only
g_spr <- pcfinhom(pp, lambda = lam_spr, kernel = "epanechnikov", correction = "translate")

pdf("figures/host_conditioned_spatstat.pdf", width = 11, height = 4)
par(mfrow = c(1, 3))
plot(E_spr,  main = "L_inhom vs spruce-host null")        # observed above band => residual foci
plot(E_full, main = "L_inhom vs spruce+total+edge null")  # collapse => host-structure, not foci
plot(g_spr,  main = "host-conditioned PCF (Epanechnikov)")
dev.off()

# ---- 6. per-block DESCRIPTIVE replicates (NO latitudinal trend; n=7 PSUs) ---
# assign each tile to a block via 15 km clustering of tile centres, then Linhom per block.
tiles[, cx := (x0 + x1)/2][, cy := (y0 + y1)/2]
tiles[, block := as.integer(factor(cutree(hclust(dist(cbind(cx, cy))), h = 15000)))]
tr <- merge(tr, tiles[, .(image_name, block)], by = "image_name")
for (b in sort(unique(tr$block))) {
  sub <- tr[block == b]
  if (nrow(sub) < 200) next
  Wb <- intersect.owin(do.call(union.owin, tile_win[tiles$block == b]), forest)
  ppb <- ppp(sub$x, sub$y, window = Wb, checkdup = FALSE)
  fb  <- ppm(ppb ~ spr, covariates = list(spr = spr))
  Lb  <- Linhom(ppb, lambda = predict(fb, type = "trend"), correction = "translate", rmax = RMAX)
  cat(sprintf("block %d: n=%d, L_inhom(500m)-500 = %.0f\n",
              b, npoints(ppb), Lb$trans[which.min(abs(Lb$r - 500))] - 500))
}
cat("\nReport per-block L_inhom as descriptive replicates with between-block spread.\n",
    "DO NOT fit a latitudinal trend (effective n = 7 blocks).\n")

# =============================================================================
# READING THE OUTPUT (locked interpretation):
#  - If observed L_inhom lies ABOVE the spruce-host envelope at 0.25-1.5 km AND stays above the
#    spruce+total+edge envelope  -> genuine residual mortality FOCI beyond host structure (the finding).
#  - If it collapses into the (richer) envelope -> mortality is host-tracking with little residual
#    clustering (a clean, publishable null).
#  - Clusters already sit on spruce-rich pixels (strong first-order host-tracking); this analysis
#    isolates the SECOND-ORDER residual.
# =============================================================================
