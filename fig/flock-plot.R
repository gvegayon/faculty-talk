library(ape)

set.seed(12312)

sizes   <- sample.int(50, 50, TRUE) + 10
keytree <- sample.int(50, 10)
graphics.off()
pdf("playground/flock-plot.pdf", width = 10, height = 5)
op <- par(mfrow = c(5, 10), mai = rep(.15,4))
for (i in 1:length(sizes)) {
  atree <- rtree(sizes[i])
  plot(
    atree, show.tip.label = FALSE,
    edge.width = 1.5,
    edge.color = ifelse(i %in% keytree, "tomato", "gray")
    )
}
dev.off()
par(op)

