beta <- 0.6
p.vec <- c(100, 1000, 10000, 100000)
resolution <- 20
r.max <- 3
r.list <- seq(from = r.max/2/resolution, to = r.max - r.max/2/resolution, by = r.max/resolution)

load(file = "weak.boundary.one-vs-two-sided.Rdata")

for (p in p.vec) {
  setEPS()
  postscript(paste0("approx_recovery_one-vs-two-sided_beta0", beta*10, 
                    "_p", format(p, scientific=F),".eps"), height = 3.5, width = 3)
  par(mar = c(2.8,2.2,1.5,0.5))
  matplot(t(drop(result.array[,,as.character(p)])), type = 'n', 
          las = 1, ylim = c(0, 1), x = r.list, col = 1,
          pch = as.character(c(1,1,2,2)), lty = c(1,2,1,2))
  theoretical.boundary <- beta
  polygon(x = c(rep(theoretical.boundary, 2), rep(r.max, 2)),
          y = c(0, 1, 1, 0), col = grey(0.9), border = F)
  lines(x = rep(theoretical.boundary, 2), y = c(0, 1), lty = 3)
  abline(h = c(0, 1), lty = 3)
  matplot(t(drop(result.array[,,as.character(p)])), type = 'b', 
          las = 1, ylim = c(0, 1), x = r.list, col = 1, add = T,
          pch = as.character(c(1,1,2,2)), lty = c(1,2,1,2))
  mtext(text = bquote('signal size' ~ r), side = 1, line = 2, at = 2.5)
  mtext(text = 'Risk approx. recovery', side = 2, line = -5.8, at = 1 * 1.1, las = 1)
  text(x = 2, y = 0.80, pos = 4, labels = bquote(beta ~ '=' ~ .(beta)))
  # text(x = 2, y = 0.72, pos = 4, labels = paste0("p = ", format(p, scientific=F)))
  text(x = 2, y = 0.72, pos = 4, labels = bquote("p =" ~ 10^.(log10(p))))
  dev.off()
}

