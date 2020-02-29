# Read in results of moretrees model -----------------------------------------
load(file = "./results/attempt4_northEast.Rdata")
sd_pm25 <- 7.48621
sd_tmmx <- 6.029194
sd_rmax <- 12.50401
n <- 3354364

# Create results table ------------------------------------------------------
OR_est <- ccs_table(root = "7", moretrees_results = moretrees_results,
                    digits = 3, mult = 10 / sd_pm25)
OR_est$long_label <- NULL
require(xtable)
row.names(OR_est) <- NULL
OR_xtable <- xtable(OR_est,align=c("l","l","p{8cm}","c","c"),digits=3,
                       display=c("d","d","s","d","f"))
names(OR_xtable) <- c("Group","CCS codes","#Outcomes","OR (95%CI)")

write(print(OR_xtable,floating=FALSE,include.rownames = FALSE),file="./figures/table1.tex")

# Create results tree plot ---------------------------------------------------
pdf(file = "./figures/tree1.pdf", width = 6, height = 3.2)
ccs_plot(root = "7", moretrees_results = moretrees_results)
dev.off()

