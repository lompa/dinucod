#####
## 

par <- par()
pal <- palette()

## you may need:
# install.packages("RColorBrewer")

Capitalize <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

## choose working directory:
path <- "..."

setwd(path)

data <- read.csv("data.csv"); dim(data)

data$baltimore <- factor(data$baltimore,
                         levels = c("dsDNA","ssDNA","dsRNA","+ssRNA","-ssRNA","+ssRNA-RT","dsDNA-RT"))
sum(is.na(data$baltimore))
# levels = c("dsDNA","ssDNA","dsRNA","+ssRNA","-ssRNA","+ssRNA-RT","dsDNA-RT",""))
# levels(data$baltimore)[nlevels(df$baltimore)] <- "unknown"

# # RColorBrewer::display.brewer.all()
# # RColorBrewer::display.brewer.all(colorblindFriendly = TRUE)
# RColorBrewer::display.brewer.all(colorblindFriendly = TRUE, select = c("Set2","Paired","Dark2"))
# 
# baltimore <- c(RColorBrewer::brewer.pal(12, "Paired")[10],
#                RColorBrewer::brewer.pal(12, "Paired")[8],
#                RColorBrewer::brewer.pal(12, "Paired")[2],
#                RColorBrewer::brewer.pal(12, "Paired")[4],
#                RColorBrewer::brewer.pal(12, "Paired")[6],
#                RColorBrewer::brewer.pal(8, "Set2")[1],
#                RColorBrewer::brewer.pal(12, "Paired")[9])

data$host <- factor(rep("",nrow(data)),
                    levels = c("", "animals", "archaea", "bacteria", "fungi", "plants", "protists"))

sort(table(data$host_cell))
# sort(table(data$host_type))
data$host[data$host_cell=="animal; fungus; plant"] <- "fungi"
data$host[data$host_cell=="fungus; plant"] <- "fungi"
data$host[data$host_cell=="animal; protist"] <- "protists"
data$host[data$host_cell=="animal; plant"] <- "plants"
data$host[data$host_cell=="protist"] <- "protists"
data$host[data$host_cell=="fungus"] <- "fungi"
data$host[data$host_cell=="plant"] <- "plants"
data$host[data$host_cell=="animal"] <- "animals"
# data$host[data$host_cell=="prokaryote"] <- "prokaryote"
sort(table(data[data$host_cell=="prokaryote",]$host_type))
data$host[data$host_type=="archaea"] <- "archaea"
data$host[data$host_type=="bacteria"] <- "bacteria"
# data[data$host_cell=="prokaryote"&data$host_type=="",]
levels(data$host)

# # RColorBrewer::display.brewer.all()
# # RColorBrewer::display.brewer.all(colorblindFriendly = TRUE)
# RColorBrewer::display.brewer.all(colorblindFriendly = TRUE, select = c("Set2","Paired","Dark2"))
# 
# host <- c("gray",
#           RColorBrewer::brewer.pal(12, "Paired")[6],
#           RColorBrewer::brewer.pal(12, "Paired")[10],
#           RColorBrewer::brewer.pal(12, "Paired")[2],
#           RColorBrewer::brewer.pal(12, "Paired")[12],
#           RColorBrewer::brewer.pal(12, "Paired")[4],
#           RColorBrewer::brewer.pal(12, "Paired")[8])

#####

dir(pattern="*.fac*")

ass <- read.table("ass.tsv")
names(ass) <- c("taxid", "assembly")

# dim(read.table("genomicSums.fac2",sep=","))
dim(read.table("cdsSums.fac2",sep=","))
# dim(read.table("genomicSums.fac3",sep=","))
# dim(read.table("cdsSums.fac3",sep=","))
# f2 <- merge(ass,read.table("genomicSums.fac2",sep=","),by.x="assembly",by.y="V1")
f2 <- merge(ass,read.table("cdsSums.fac2",sep=","),by.x="assembly",by.y="V1")
# f3 <- merge(ass,read.table("genomicSums.fac3",sep=","),by.x="assembly",by.y="V1")
# f3 <- merge(ass,read.table("cdsSums.fac3",sep=","),by.x="assembly",by.y="V1")
names(f2) <- c("assembly","taxid",sort(gsub("[.]", "", levels(interaction(letters[c(1,3,7,20)],letters[c(1,3,7,20)])))))
# names(tn) <- c("assembly","taxid",sort(gsub("[.]", "", levels(interaction(letters[c(1,3,7,20)],letters[c(1,3,7,20)],letters[c(1,3,7,20)])))))

data0 <- data

data <- merge(data,f2,by="assembly"); dim(data)
# data <- merge(data,f2,by="taxid"); dim(data)
# data <- merge(data,f3,by="assembly"); dim(data)
# data <- merge(data,f3,by="taxid"); dim(data)

#####

par <- par()
par(family = "serif")
par(oma = c(0,0,0,0))
par(mar = c(4,4,2,3))# + 0.1)
par(mfrow = c(1,1))

names(data)[35] <- "gc."
names(data)[119] <- "gc"

#####
## Dinucleotide frequencies

names(data)[110:125]

data[110:125] <- data[110:125]/rowSums(data[110:125])

# with(data,boxplot(cg~baltimore))
# with(data,vioplot::vioplot(cg~baltimore,col="gray"))

# with(data,boxplot(cg~host))
# with(data,vioplot::vioplot(cg~host,col="gray"))

pca <- prcomp(data[,110:125], scale.=TRUE)
plot(pca)
# axis(1, at = 1:10, labels = paste("PC", 1:10), las=1)
# biplot(pca, col = adjustcolor(c("white", "black"),0.5), scale=0.83, xaxt="n", yaxt="n", ann=FALSE)
mat <- summary(pca)
mat$importance[3, 2]

palette(adjustcolor(baltimore,0.5))

1 -> c
# 1.5 -> c
# par.bak <- par()
# par(mfrow = c(3, 3))
#png(paste0("pca_",h,"_loadings.#png"))
biplot(pca, col = c("white", "black"), pc.biplot = TRUE, main = "", cex = c(1, c),
       xlab="", ylab="", xaxt="n", yaxt="n")
title(main="Dinucleotide frequencies\nby Baltimore")
mtext(paste0("PC1 (",format(round(mat$importance[2, 1], 2), nsmall = 2),")"),
      side = 1, line = 0.5, cex.lab = 1, las = 0, col = "black")
mtext(paste0("PC2 (",format(round(mat$importance[2, 2], 2), nsmall = 2),")"),
      side = 2, line = 0.5, cex.lab = 1, las = 0, col = "black")
abline(v=0, lty=2, col = adjustcolor("black",0.5))
abline(h=0, lty=2, col = adjustcolor("black",0.5))
#dev.off()
g = 0
# prcomp(data[,110:125]
for(i in levels(data[, "baltimore"])){
  g = g + 1
  #png(paste0("pca_",h,"_group",g,"_",i,".#png"))
  # plot(pca$x[, 1:2], type = "n", main = paste0(i,"\n", format(round(mat$importance[3, 2], 2), nsmall = 2)),
  plot(pca$x[, 1:2], type = "n", main = i,
       xlab = paste0("PC1 (",format(round(mat$importance[2, 1], 2), nsmall = 2),")"),
       ylab = paste0("PC2 (",format(round(mat$importance[2, 2], 2), nsmall = 2),")"))
  points(pca$x[data[complete.cases(data[,110:125]), "baltimore"] != i, 1:2], col = adjustcolor("gray", 0.5), pch = 20)
  points(pca$x[data[complete.cases(data[,110:125]), "baltimore"] == i, 1:2], col = adjustcolor("black", 0.5), pch = 20)
  title(paste0("n = ",format(nrow(na.omit(pca$x[data[complete.cases(data[,110:125]), "baltimore"] == i, 1:2])),big.mark=",")),adj=1,font.main=1, cex.main=1)
  # points(pca$x[data[, "baltimore"] == i, 1:2], col = adjustcolor(which(levels(data[, "baltimore"])==i), 1), pch = 20)
  # if(i == "dsDNA-RT"){
  #   points(pca$x[data[, "baltimore"] == i, 1:2], col = adjustcolor(which(levels(data[, "baltimore"])==i), 1), pch = 20)
  # }
  #dev.off()
}

palette(adjustcolor(host[-which(levels(data[, "host"])=="")],0.5))

1 -> c
# 1.5 -> c
# par.bak <- par()
# par(mfrow = c(3, 3))
#png(paste0("pca_",h,"_loadings.#png"))
biplot(pca, col = c("white", "black"), pc.biplot = TRUE, main = "", cex = c(1, c),
       xlab="", ylab="", xaxt="n", yaxt="n")
title(main="Dinucleotide frequencies\nby host")
mtext(paste0("PC1 (",format(round(mat$importance[2, 1], 2), nsmall = 2),")"),
      side = 1, line = 0.5, cex.lab = 1, las = 0, col = "black")
mtext(paste0("PC2 (",format(round(mat$importance[2, 2], 2), nsmall = 2),")"),
      side = 2, line = 0.5, cex.lab = 1, las = 0, col = "black")
abline(v=0, lty=2, col = adjustcolor("black",0.5))
abline(h=0, lty=2, col = adjustcolor("black",0.5))
#dev.off()
g = 0
# prcomp(data[,110:125]
for(i in levels(data[, "host"])[-which(levels(data[, "host"])=="")]){
  g = g + 1
  #png(paste0("pca_",h,"_group",g,"_",i,".#png"))
  # plot(pca$x[, 1:2], type = "n", main = paste0(i,"\n", format(round(mat$importance[3, 2], 2), nsmall = 2)),
  plot(pca$x[, 1:2], type = "n", main = Capitalize(i),
       xlab = paste0("PC1 (",format(round(mat$importance[2, 1], 2), nsmall = 2),")"),
       ylab = paste0("PC2 (",format(round(mat$importance[2, 2], 2), nsmall = 2),")"))
  points(pca$x[data[complete.cases(data[,110:125]), "host"] != i, 1:2], col = adjustcolor("gray", 0.5), pch = 20)
  points(pca$x[data[complete.cases(data[,110:125]), "host"] == i, 1:2], col = adjustcolor("black", 0.5), pch = 20)
  title(paste0("n = ",format(nrow(na.omit(pca$x[data[complete.cases(data[,110:125]), "host"] == i, 1:2])),big.mark=",")),adj=1,font.main=1, cex.main=1)
  # points(pca$x[data[, "host"] == i, 1:2], col = adjustcolor(which(levels(data[, "host"])==i)-1, 1), pch = 20)
  # if(i == "protist"){
  #   points(pca$x[data[, "host"] == i, 1:2], col = adjustcolor(which(levels(data[, "host"])==i)-1, 1), pch = 20)
  # }
  #dev.off()
}

#####
## Codon frequencies

names(data)[44:107]

data[44:107] <- data[44:107]/rowSums(data[44:107])

pca <- prcomp(data[complete.cases(data[,44:107]),44:107], scale.=TRUE)
plot(pca)
# biplot(pca, col = adjustcolor(c("white", "black"),0.5), scale=0.83, xaxt="n", yaxt="n", ann=FALSE)
mat <- summary(pca)
mat$importance[3, 2]

palette(adjustcolor(baltimore,0.5))

1 -> c
# 1.5 -> c
# par.bak <- par()
# par(mfrow = c(3, 3))
#png(paste0("pca_",h,"_loadings.#png"))
biplot(pca, col = c("white", "black"), pc.biplot = TRUE, main = "", cex = c(1, c),
       xlab="", ylab="", xaxt="n", yaxt="n")
title(main="Codon frequencies\nby Baltimore")
mtext(paste0("PC1 (",format(round(mat$importance[2, 1], 2), nsmall = 2),")"),
      side = 1, line = 0.5, cex.lab = 1, las = 0, col = "black")
mtext(paste0("PC2 (",format(round(mat$importance[2, 2], 2), nsmall = 2),")"),
      side = 2, line = 0.5, cex.lab = 1, las = 0, col = "black")
abline(v=0, lty=2, col = adjustcolor("black",0.5))
abline(h=0, lty=2, col = adjustcolor("black",0.5))
#dev.off()
g = 0
# prcomp(data[,110:125]
for(i in levels(data[complete.cases(data[,44:107]), "baltimore"])){
  g = g + 1
  #png(paste0("pca_",h,"_group",g,"_",i,".#png"))
  # plot(pca$x[, 1:2], type = "n", main = paste0(i,"\n", format(round(mat$importance[3, 2], 2), nsmall = 2)),
  plot(pca$x[, 1:2], type = "n", main = i, yaxt = "n",
       xlab = paste0("PC1 (",format(round(mat$importance[2, 1], 2), nsmall = 2),")"),
       ylab = paste0("PC2 (",format(round(mat$importance[2, 2], 2), nsmall = 2),")"))
  axis(2, at = seq(-10,10,5), las=2)
  points(pca$x[data[complete.cases(data[,44:107]), "baltimore"] != i, 1:2], col = adjustcolor("gray", 0.5), pch = 20)
  points(pca$x[data[complete.cases(data[,44:107]), "baltimore"] == i, 1:2], col = adjustcolor("black", 0.5), pch = 20)
  title(paste0("n = ",format(nrow(na.omit(pca$x[data[complete.cases(data[,44:107]), "baltimore"] == i, 1:2])),big.mark=",")),adj=1,font.main=1, cex.main=1)
  # points(pca$x[data[complete.cases(data[,44:107]), "baltimore"] == i, 1:2], col = adjustcolor(which(levels(data[complete.cases(data[,44:107]), "baltimore"])==i), 1), pch = 20)
  # if(i == "dsDNA-RT"){
  #   points(pca$x[data[complete.cases(data[,44:107]), "baltimore"] == i, 1:2], col = adjustcolor(which(levels(data[complete.cases(data[,44:107]), "baltimore"])==i), 1), pch = 20)
  # }
  #dev.off()
}

palette(adjustcolor(host[-which(levels(data[, "host"])=="")],0.5))

1 -> c
# 1.5 -> c
# par.bak <- par()
# par(mfrow = c(3, 3))
#png(paste0("pca_",h,"_loadings.#png"))
biplot(pca, col = c("white", "black"), pc.biplot = TRUE, main = "", cex = c(1, c),
       xlab="", ylab="", xaxt="n", yaxt="n")
title(main="Codon frequencies\nby host")
mtext(paste0("PC1 (",format(round(mat$importance[2, 1], 2), nsmall = 2),")"),
      side = 1, line = 0.5, cex.lab = 1, las = 0, col = "black")
mtext(paste0("PC2 (",format(round(mat$importance[2, 2], 2), nsmall = 2),")"),
      side = 2, line = 0.5, cex.lab = 1, las = 0, col = "black")
abline(v=0, lty=2, col = adjustcolor("black",0.5))
abline(h=0, lty=2, col = adjustcolor("black",0.5))
#dev.off()
g = 0
# prcomp(data[,110:125]
for(i in levels(data[complete.cases(data[,44:107]), "host"])[-which(levels(data[, "host"])=="")]){
  g = g + 1
  #png(paste0("pca_",h,"_group",g,"_",i,".#png"))
  # plot(pca$x[, 1:2], type = "n", main = paste0(i,"\n", format(round(mat$importance[3, 2], 2), nsmall = 2)),
  plot(pca$x[, 1:2], type = "n", main = Capitalize(i), yaxt = "n",
       xlab = paste0("PC1 (",format(round(mat$importance[2, 1], 2), nsmall = 2),")"),
       ylab = paste0("PC2 (",format(round(mat$importance[2, 2], 2), nsmall = 2),")"))
  axis(2, at = seq(-10,10,5), las=2)
  points(pca$x[data[complete.cases(data[,44:107]), "host"] != i, 1:2], col = adjustcolor("gray", 0.5), pch = 20)
  points(pca$x[data[complete.cases(data[,44:107]), "host"] == i, 1:2], col = adjustcolor("black", 0.5), pch = 20)
  title(paste0("n = ",format(nrow(na.omit(pca$x[data[complete.cases(data[,44:107]), "host"] == i, 1:2])),big.mark=",")),adj=1,font.main=1, cex.main=1)
  # points(pca$x[data[complete.cases(data[,44:107]), "host"] == i, 1:2], col = adjustcolor(which(levels(data[complete.cases(data[,44:107]), "host"])==i)-1, 1), pch = 20)
  # if(i == "protist"){
  #   points(pca$x[data[complete.cases(data[,44:107]), "host"] == i, 1:2], col = adjustcolor(which(levels(data[complete.cases(data[,44:107]), "host"])==i)-1, 1), pch = 20)
  # }
  #dev.off()
}
