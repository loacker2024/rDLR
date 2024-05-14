library(ggplot2)
#setwd("/Volumes/PS2000/ubuntu/rDLR_final_code/simulation/")

res = readRDS("./illustration_fig_dat.rds")
names(res)

dat = res$part1
head(dat)
xlab="" 
ylab = "Coefficients"
gp11 = ggplot(dat, aes(X)) + geom_line(aes(y = X1, colour = "X1",size=0.5)) + 
  xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(a1)", x=min(dat$X),y=0.8*max(dat[,c(2,3,4,5,6)]),hjust=.2, color = "black", size = 26/.pt)

ylab = ""
gp12 = ggplot(dat, aes(X)) + geom_line(aes(y = X2, colour = "X2",size=0.5)) +
  xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(a2)", x=min(dat$X),y=0.8*max(dat[,c(2,3,4,5,6)]),hjust=.2, color = "black", size = 26/.pt)
gp13 = ggplot(dat, aes(X)) +  geom_line(aes(y = X3, colour = "X3",size=0.5))  +
  xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(a3)", x=min(dat$X),y=0.8*max(dat[,c(2,3,4,5,6)]),hjust=.2, color = "black", size = 26/.pt)
gp14 = ggplot(dat, aes(X)) +  geom_line(aes(y = X4, colour = "X4",size=0.5)) +
  xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(a4)", x=min(dat$X),y=0.8*max(dat[,c(2,3,4,5,6)]),hjust=.2, color = "black", size = 26/.pt)
gp15 = ggplot(dat, aes(X)) +  geom_line(aes(y = X5, colour = "X5",size=0.5)) +
  xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(a5)", x=min(dat$X),y=0.8*max(dat[,c(2,3,4,5,6)]),hjust=.2, color = "black", size = 26/.pt)

dat = res$part2
xlab=""; ylab = "Coefficients"
gp21 = ggplot(dat, aes(X)) + geom_line(aes(y = X1, colour = "X1",size=0.5)) + xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(b1)", x=min(dat$X),y=0.2*max(dat[,c(2,3,4,5,6)]),hjust=.2, color = "black", size = 26/.pt) + ylim(0,3)
ylab ="" #"Coefficients"
gp22 = ggplot(dat, aes(X)) +  geom_line(aes(y = X2, colour = "X2",size=0.5)) + xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(b2)", x=min(dat$X),y=0.8*max(dat[,c(2,3,4,5,6)]),hjust=.2, color = "black", size = 26/.pt)

gp23 = ggplot(dat, aes(X)) +   geom_line(aes(y = X3, colour = "X3",size=0.5)) + xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(b3)", x=min(dat$X),y=0.8*max(dat[,c(2,3,4,5,6)]),hjust=.2, color = "black", size = 26/.pt)

gp24 = ggplot(dat, aes(X)) +   geom_line(aes(y = X4, colour = "X4",size=0.5)) + xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(b4)", x=min(dat$X),y=0.8*max(dat[,c(2,3,4,5,6)]),hjust=.2, color = "black", size = 26/.pt)

gp25 = ggplot(dat, aes(X)) +   geom_line(aes(y = X5, colour = "X5",size=0.5)) +
  xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(b5)", x=min(dat$X),y=0.8*max(dat[,c(2,3,4,5,6)]),hjust=.2, color = "black", size = 26/.pt)

dat = res$part3
gp31 = ggplot(dat, aes(X)) + geom_line(aes(y = X1, colour = "X1",size=0.5)) + xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(c1)", x=min(dat$X),y=0.8,hjust=.2, color = "black", size = 26/.pt) + ylim(-1,1)
ylab ="" 
gp32 = ggplot(dat, aes(X)) +  geom_line(aes(y = X2, colour = "X2",size=0.5)) + xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(c2)", x=min(dat$X),y=0.8,hjust=.2, color = "black", size = 26/.pt) + ylim(-1,1)

gp33 = ggplot(dat, aes(X)) +   geom_line(aes(y = X3, colour = "X3",size=0.5)) + xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(c3)", x=min(dat$X),y=0.8,hjust=.2, color = "black", size = 26/.pt) + ylim(-1,1)

gp34 = ggplot(dat, aes(X)) +   geom_line(aes(y = X4, colour = "X4",size=0.5)) + xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(c4)", x=min(dat$X),y=0.8,hjust=.2, color = "black", size = 26/.pt) + ylim(-1,1)

gp35 = ggplot(dat, aes(X)) +   geom_line(aes(y = X5, colour = "X5",size=0.5)) +
  xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(c5)", x=min(dat$X),y=0.8,hjust=.2, color = "black", size = 26/.pt) + ylim(-1,1)

dat = res$part4
xlab = "Time"
gp41 = ggplot(dat, aes(X)) + geom_line(aes(y = X1, colour = "X1",size=0.5)) + xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(d1)", x=min(dat$X),y=0.2*max(dat[,c(2,3,4,5,6)]),hjust=.2, color = "black", size = 26/.pt) + ylim(0,3)
ylab ="" #"Coefficients"
gp42 = ggplot(dat, aes(X)) +  geom_line(aes(y = X2, colour = "X2",size=0.5)) + xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(d2)", x=min(dat$X),y=0.8*max(dat[,c(2,3,4,5,6)]),hjust=.2, color = "black", size = 26/.pt) + ylim(0,2.5)

gp43 = ggplot(dat, aes(X)) +   geom_line(aes(y = X3, colour = "X3",size=0.5)) + xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(d3)", x=min(dat$X),y=0.8*max(dat[,c(2,3,4,5,6)]),hjust=.2, color = "black", size = 26/.pt) + ylim(0,2.5)

gp44 = ggplot(dat, aes(X)) +   geom_line(aes(y = X4, colour = "X4",size=0.5)) + xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(d4)", x=min(dat$X),y=0.8,hjust=.2, color = "black", size = 26/.pt) + ylim(-1,1)

gp45 = ggplot(dat, aes(X)) +   geom_line(aes(y = X5, colour = "X5",size=0.5)) +
  xlab(xlab) + ylab(ylab) + theme(legend.position="none", text = element_text(size=26)) + theme(panel.background = element_rect(fill = "white", colour = "black")) +
  annotate("text", label = "(d5)", x=min(dat$X),y=0.8,hjust=.2, color = "black", size = 26/.pt) + ylim(-1,1)

pdf("illustration_simulated_coefficients_2024.pdf", width = 20, height = 24) 
gridExtra::grid.arrange(gp11, gp12, gp13, gp14, gp15, 
                        gp21, gp22, gp23, gp24, gp25,
                        gp31, gp32, gp33, gp34, gp35,
                        gp41, gp42, gp43, gp44, gp45, nrow = 4, ncol=5)
dev.off()




