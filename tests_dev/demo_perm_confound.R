aggSum <- function(gg){ D <- model.matrix(~gg-1); colnames(D) <- levels(gg); t(D) %*% (v %*% D) }
aggAvg <- function(gg){ D <- model.matrix(~gg-1); colnames(D) <- levels(gg); S <- t(D)%*%(v%*%D); B <- t(D)%*%(1*(v>0))%*%D; S/B }

set.seed(1)
n <- 200L
pos <- c(runif(100, 0, 40), runif(100, 60, 100))   # two spatial domains
g  <- factor(c(rep("A",100), rep("B",100)))
v  <- exp(-as.matrix(dist(pos))/3); diag(v) <- 0    # pure distance decay, no LR specificity


X  <- aggSum(g); XA <- aggAvg(g)
set.seed(2)
R1s <- replicate(2000, aggSum(sample(g)))
R1a <- replicate(2000, aggAvg(sample(g)))
ord <- order(pos); lab <- as.character(g)[ord]
set.seed(4)
R2s <- replicate(2000, { s <- sample.int(n,1); aggSum(factor(lab[(seq_len(n)-s) %% n + 1L], levels=c("A","B"))) })

cat("== sum ==\n")
cat("P_AA obs =", signif(X["A","A"],3), " null(random) mean =", signif(mean(R1s["A","A",]),3),
    " p =", mean(R1s["A","A",] >= X["A","A"]), " p(shift-null) =", mean(R2s["A","A",] >= X["A","A"]), "\n")
cat("P_BB obs =", signif(X["B","B"],3), " null(random) mean =", signif(mean(R1s["B","B",]),3),
    " p =", mean(R1s["B","B",] >= X["B","B"]), "\n")
cat("P_AB obs =", signif(X["A","B"],3), " null(random) mean =", signif(mean(R1s["A","B",]),3),
    " p =", mean(R1s["A","B",] >= X["A","B"]), "\n")
cat("== avg ==\n")
cat("A_AA obs =", signif(XA["A","A"],3), " null(random) mean =", signif(mean(R1a["A","A",]),3),
    " p =", mean(R1a["A","A",] >= XA["A","A"]), "\n")
cat("A_AB obs =", signif(XA["A","B"],3), " null(random) mean =", signif(mean(R1a["A","B",], na.rm=TRUE),3),
    " p =", mean(R1a["A","B",] >= XA["A","B"], na.rm=TRUE), "\n")
