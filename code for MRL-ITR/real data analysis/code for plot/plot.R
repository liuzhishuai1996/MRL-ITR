rm(list=ls())
load("4-Result4.RData")
x <- c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800)
y1 <- Result[1, ]
y2 <- Result[2, ]
y3 <- Result[4, ]
sd1 <- Result[3, ]
sd2 <- Result[5, ]
plot(x, y1, 'l', ylim = c(100, 800), ylab = 'MRL', xlab = 'Time', lty = 2)
points(x, y1, pch = 20)
lines(x, y2, 'l', col = 'red', lty = 3)
points(x, y2, pch = 20, col = 'red')
# points (x, y2 + sd1, pch = 5)
for (i in 1:length(x)){
  segments(x0 = x[i], y0 = y2[i] + 1.96*sd1[i], y1 = y2[i] - 1.96*sd1[i], lty = 1, col = 'red', lwd = 1.5)
  segments(x0 = x[i] - 2, x1 = x[i] + 2, y0 = y2[i] + 1.96*sd1[i], lty = 1, col = 'red', lwd = 1.5)
  segments(x0 = x[i] - 2, x1 = x[i] + 2, y0 = y2[i] - 1.96*sd1[i], lty = 1, col = 'red', lwd = 1.5)
}

# plot(x, y1, 'l', ylim = c(100, 800), ylab = 'MRL', xlab = 'Time')
lines(x, y3, 'l', col = 'blue', lty = 4)
points(x, y3, pch = 20, col = 'blue')
for (i in 1:length(x)){
  segments(x0 = x[i], y0 = y3[i] + 1.96*sd2[i], y1 = y3[i] - 1.96*sd2[i], lty = 3, col = 'blue', lwd = 1.5)
  segments(x0 = x[i] - 2, x1 = x[i] + 2, y0 = y3[i] + 1.96*sd2[i], lty = 1, col = 'blue', lwd = 1.5)
  segments(x0 = x[i] - 2, x1 = x[i] + 2, y0 = y3[i] - 1.96*sd2[i], lty = 1, col = 'blue', lwd = 1.5)
}

legend(500, 800, c('IPW', 'S-IPW', 'Randomized'), lty = c(2, 4, 3), col = c('red', 'blue', 'black'))
# legend(400, 700, c('Randomized Treatment', 'estimated optimal regime'), col =2:3)
# plot(x = x, y = y2-y1, 'l', ylab = 'MRL', xlab = 'Time')
# for (i in 1:length(x)){
#   segments(x0 = x[i], y0 = y2[i] + sd1[i], y1 = y2[i] - sd1[i], lty = 3)
# }