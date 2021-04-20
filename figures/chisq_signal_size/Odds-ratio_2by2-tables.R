
OR <- 1.3^(-50:50)

# calculate signal size for a range of odds ratios, given marginal prob.
lam.finder <- function(p, q) {
  dmin <- max(-p*q, -(1-p)*(1-q), (1-p)*q-1, (1-q)*p-1)
  dmax <- min(1-p*q, 1-(1-p)*(1-q), (1-p)*q, (1-q)*p)
  d <- vector(mode = "numeric", length = 101)
  d.root.finder <- function(d, OR){(p*q+d)*((1-p)*(1-q)+d)/(p*(1-q)-d)/((1-p)*q-d) - OR}
  for (i in 1:length(OR)) {
    solution <- uniroot(f = d.root.finder, lower = dmin, upper = dmax, OR = OR[i])
    d[i] <- solution$root
  }
  lam <- d^2*(1/(p*q) + 1/(p*(1-q)) + 1/((1-p)*q) + 1/((1-p)*(1-q)))
  return(lam)
}


#### Figure 1: p = 1/2 ####

library(MASS)

setEPS()
postscript("singal-vs-odds-p05.eps", height = 4, width = 4)
p <- 1/2
q <- 1/2
lam <- lam.finder(p, q)
par(mar = c(3,3,2,1))
plot(x = OR, y = lam, log = 'x', type = 'n', las = 1, ylim = c(0,1), xaxt = 'n')
mtext(text = expression(paste("signal size ", w^2, '= ', lambda/n)), 
      side = 2, at = 1.1, las = 1, line = -5)
axis(side = 1, at = 10^seq(-6, 6, 2), 
     labels = parse(text = paste("10^", seq(-6, 6, 2), sep="")))
mtext(text = "odds ratio", side = 1, at = 3e5, las = 1, line = 2)

q.vec <- as.character(fractions(1/c(2, 3, 5, 10)))
for (i in 1:length(q.vec)) {
  q <- eval(parse(text = q.vec[i]))
  lam <- lam.finder(p, q)
  lines(x = OR, y = lam, col = 1, lty = i)
  text(x = 1e6, y = 0.95*lam[1] - 0.04, pos = 2,
       labels = bquote(paste(theta[1], " = ", .(q.vec[i]), ", ", phi[1], " = 1/2", sep="")))
}
dev.off()

#### Figure 2: p = 1/3 ####

setEPS()
postscript("singal-vs-odds-p0333.eps", height = 4, width = 4)
p <- 1/3
q <- 1/3
lam <- lam.finder(p, q)
par(mar = c(3,3,2,1))
plot(x = OR, y = lam, log = 'x', type = 'n', las = 1, ylim = c(0,1), xaxt = 'n')
mtext(text = expression(paste("signal size ", w^2, '= ', lambda/n)), 
      side = 2, at = 1.1, las = 1, line = -5)
axis(side = 1, at = 10^seq(-6, 6, 2), 
     labels = parse(text = paste("10^", seq(-6, 6, 2), sep="")))
mtext(text = "odds ratio", side = 1, at = 3e5, las = 1, line = 2)

q.vec <- as.character(fractions(1/c(3, 5, 10)))
for (i in 1:length(q.vec)) {
  q <- eval(parse(text = q.vec[i]))
  lam <- lam.finder(p, q)
  lines(x = OR, y = lam, col = 1, lty = i)
  text(x = 1e6, y = 0.95*lam[100] - 0.04, pos = 2,
       labels = bquote(paste(theta[1], " = ", .(q.vec[i]), ", ", phi[1], " = 1/3", sep="")))
}
dev.off()



#### analytical solution ####

signal.size.ana.sol <- function(p, q, R){
  A <- p*q*(1-p)*(1-q)
  B <- (p*q+(1-p)*(1-q)) 
  C <- (p*(1-q)+(1-p)*q)
  lam <- 1/(2*(R-1))^2 * 
    (B+C*R - sqrt((B+C*R)^2 - 4*A*(R-1)^2))^2 * 
    (1/(p*q)+1/((1-p)*q)+1/(p*(1-q))+1/((1-p)*(1-q)))
  ifelse(R==1, 0, lam)
}

plot(OR, signal.size.ana.sol(1/2,1/2,OR), xlim = c(1e-6, 1e6), ylim = c(0, 1),
     log = 'x', type = 'l', las = 1)
q <- 1/4
for (p in c(1/5, 1/4, 1/3, 1/2, 2/3, 3/4, 4/5)) {
  lines(OR, signal.size.ana.sol(p, q, OR), lty = 2, col = 1)
  abline(h = q*p/(1-p)/(1-q), lty = 2)
  abline(h = (1-p)*(1-q)/q/p, lty = 3)
}


plot(OR, signal.size.ana.sol(1/2,1/2,OR), xlim = c(1e-6, 1e6), ylim = c(0, 1),
     log = 'x', type = 'l', las = 1)
q <- 1/4
for (p in c(1/5, 1/4, 1/3, 1/2, 2/3, 3/4, 4/5)) {
  lines(OR, signal.size.ana.sol(p, q, OR), lty = 2, col = 1)
  abline(h = p*(1-q)/(1-p)/q, lty = 2)
  abline(h = q*(1-p)/(1-q)/p, lty = 3)
}


#### Figure 3: p = 1/100 ####

setEPS()
postscript("singal-vs-odds-p001.eps", height = 5, width = 5)
p <- 1/100
q <- 1/2
par(mar = c(3,3,2,1))
plot(OR, signal.size.ana.sol(p,q,OR), xlim = c(1e-6, 1e6), ylim = c(0, 1/(1/0.1+1)),
     log = 'x', type = 'l', las = 1,  xaxt = 'n')
mtext(text = expression(paste("signal size ", lambda/n)), 
      side = 2, at = 0.1, las = 1, line = -3.3)
axis(side = 1, at = 10^seq(-6, 6, 2), 
     labels = parse(text = paste("10^", seq(-6, 6, 2), sep="")))
mtext(text = "odds ratio", side = 1, at = 3e5, las = 1, line = 2)

q.vec <- as.character(fractions(1/c(2, 3, 5, 10)))
for (i in 1:length(q.vec)) {
  q <- eval(parse(text = q.vec[i]))
  lam <- signal.size.ana.sol(p, q, OR)
  lines(x = OR, y = lam, col = 1, lty = i)
  text(x = 1e6, y = 0.99*lam[100] - 0.005, pos = 2,
       labels = bquote(paste(theta[1], " = ", .(q.vec[i]), ", ", phi[1], " = 1/100", sep="")))
}
dev.off()


#### solve for ORs and RAFs at given signal size  ####

OR.finder <- function(q = 1/2, signal.size) {
  p.min <- (signal.size*q)/((1-q)+signal.size*q)
  p.mid <- q
  p.max <- (q)/(signal.size*(1-q)+q)
  OR.root.finder <- function(OR, p, q, signal.size){signal.size.ana.sol(p, q, OR) - signal.size}
  # left branch
  p.vec.left <- exp(log(p.min) + 1:100/100 * (log(p.mid) - log(p.min)))
  OR.vec.left <- vector(mode = "numeric", length = 100)
  for (i in 1:length(p.vec.left)) {
    solution <- uniroot(f = OR.root.finder, interval = c(1.01, 1e6),
                        p = p.vec.left[i], q = q, signal.size)
    OR.vec.left[i] <- solution$root
  }
  # right branch
  p.vec.right <- rev(1 - exp(log(1-p.max) + 1:100/100 * (log(1-p.mid) - log(1-p.max))))
  OR.vec.right <- vector(mode = "numeric", length = 100)
  for (i in 1:length(p.vec.right)) {
    solution <- uniroot(f = OR.root.finder, interval = c(1.01, 1e6),
                        p = p.vec.right[i], q = q, signal.size)
    OR.vec.right[i] <- solution$root
  }
  
  return(list(p = c(p.min, p.vec.left, p.vec.right, p.max),
              R = c(1e8, OR.vec.left, OR.vec.right, 1e8)))
}


#### plot "equi-power curve" ####

x.adj <- function(x) {ifelse(x < 0.5, log(x), -log(1-x) - log(4))}
# x.adj <- function(x) {x}
y.adj <- function(y) {y - 0.9}

q <- 1/2
solution.vec <- OR.finder(q = q, signal.size = 0.05)
plot(x = x.adj(solution.vec$p), y = y.adj(solution.vec$R), 
     log = 'y', type = 'n', las = 1, xaxt = 'n', yaxt = 'n', 
     xlim = x.adj(c(1e-4, 1-5e-4)), ylim = y.adj(c(1, y.plot.max)),
     xaxs="i", yaxs="i")
# abline(v = 1/(1/0.05+1), lty = 2)
axis(side = 1, at = x.adj(c(0.5, 10^(-1:-4), 1-10^(-1:-3))), 
     labels = c(0.5, 0.1, expression(10^{-2}), expression(10^{-3}), expression(10^{-4}),
                0.9, 0.99, 0.999))
axis(side = 2, at = y.adj(c(1, 2, 5, 10, 20, 50, 100)), 
     labels = c(1, 2, 5, 10, 20, 50, 100), las = 1)

for (signal.size in c(10^(-5:-1), 5*10^(-5:-2))) {
  solution.vec <- OR.finder(q = q, signal.size)
  lines(x = x.adj(solution.vec$p), y = y.adj(solution.vec$R), lty = 2)
}


#### save "equi-power curve" for p < 1/2, q = 1/2 ####

setEPS()
postscript("odds-ratio-allele-frequency-maxOR1000.eps", height = 5, width = 6)
par(mar = c(3,3,2,1))
options(scipen=0)
y.plot.max <- 1000
plot(x = solution.vec$p, y = solution.vec$R, 
     log = 'xy', type = 'n', las = 1, xaxt = 'n', xaxs="i", yaxs="i",
     xlim = c(1e-3, 1/2), ylim = c(1, y.plot.max))
mtext(text = "odds ratio", 
      side = 2, at = y.plot.max*1.6, las = 1, line = -1.5)
# axis(side = 1, at = 10^(-4:-1), 
#     labels = parse(text = paste("10^", seq(-4, -1, 1), sep="")))
# axis(side = 1, at = 5*10^(-3:-1), 
#      labels = expression(5%*%10^{-1}, 5%*%10^{-2}, 5%*%10^{-3}))
# axis(side = 1, at = c(5*10^(-1:-3), 10^(-3)), 
#     labels = expression(5%*%10^{-1}, 5%*%10^{-2}, 5%*%10^{-3}, 10^{-3}))
axis(side = 1, at = c(5*10^(-1:-3), 10^(-3)), 
     labels = c(0.5, 0.05, 0.005, 0.001))
mtext(text = "allele frequency", side = 1, at = 0.3, las = 1, line = 2)

for (signal.size in c(10^(-4:-1), 5*10^(-4:-2))) {
  solution.vec <- OR.finder(q = 1/2, signal.size)
  lines(x = solution.vec$p, y = solution.vec$R, lty = 2)
  if (signal.size >= 1e-3){
    if (signal.size == 1e-3) {
      text(x = solution.vec$p[1]*1.03, y = y.plot.max*0.7, pos = 4, 
           labels = bquote(paste(lambda/n, " = ", .(signal.size))))
    } else {
      text(x = solution.vec$p[1]*1.03, y = y.plot.max*0.7, pos = 4, labels = signal.size)
    }
  }
}
dev.off()


#### plotting both halves ####

# x.adj <- function(x) {tan((x-0.5)*pi)} # doesn't work
x.adj <- function(x) {ifelse(x < 0.5, log(x), -log(1-x) - log(4))}
# x.adj <- function(x) {x}
y.adj <- function(y) {y - 0.9}

par(mar = c(3,3,2,1))
options(scipen=0)
y.plot.max <- 100
plot(x = x.adj(solution.vec$p), y = y.adj(solution.vec$R), 
     log = 'y', type = 'n', las = 1, xaxt = 'n', yaxt = 'n', 
     xlim = x.adj(c(1e-4, 1-5e-3)), ylim = y.adj(c(1, y.plot.max)),
     xaxs="i", yaxs="i")
mtext(text = "odds ratio", 
      side = 2, at = y.plot.max*1.6, las = 1, line = -1.5)
axis(side = 1, at = x.adj(c(0.5, 10^(-1:-4), 1-10^(-1:-2))), 
     labels = c(0.5, 0.1, expression(10^{-2}), expression(10^{-3}), expression(10^{-4}),
                0.9, 0.99))
axis(side = 2, at = y.adj(c(1, 2, 5, 10, 20, 50, 100)), 
     labels = c(1, 2, 5, 10, 20, 50, 100), las = 1)
mtext(text = "risk allele frequency", side = 1, at = x.adj(0.97), las = 1, line = 2)

for (signal.size in c(10^(-4:-1), 5*10^(-4:-2))) {
  solution.vec <- OR.finder(q = 1/2, signal.size)
  lines(x = x.adj(solution.vec$p), y = y.adj(solution.vec$R), lty = 2)
  #lines(x = x.adj(1-solution.vec$p), y = y.adj(solution.vec$R), lty = 2)
  if (signal.size >= 1e-4){
    if (signal.size == 1e-1) {
      text(x = x.adj(solution.vec$p[1]*1.03), y = y.plot.max*0.7, pos = 4, 
           labels = bquote(paste(lambda/n, " = ", 10^{.(log10(signal.size))})))
    } else if (signal.size %in% 10^(-4:-2)) {
      text(x = x.adj(solution.vec$p[1]*0.9), y = y.plot.max*0.7, pos = 4, 
           labels = bquote(10^{.(log10(signal.size))}))
    }
  }
}


# rbPal <- colorRampPalette(c('red','blue'))
# Col <- rbPal(10)[as.numeric(cut(log((log(BC$PVALUE_MLOG))), breaks = 10))]
#text(x = x.adj(0.1), y = 5, labels = "Draft Figure", cex = 5, pos = 1, col = "grey70")
#text(x = x.adj(0.1), y = 1, labels = "Feb 13, 2019", cex = 5, pos = 1, col = "grey70")
