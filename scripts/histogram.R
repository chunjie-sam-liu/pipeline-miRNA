com_pdf<-$outdir."expression_heatmap.pdf";
diff_pdf<-$outdir."different_expressed_miRNA.pdf";
s<- data.level1 <- read.table("$outfile", head=T, sep="\\t",row.names=1)$#files+1;

n<-(1:$s)*2
colors<- colorRampPalette(c("white", "darkblue"))(100)
data.level1 <- read.table("$outfile", head=T, sep="\\t",row.names=1)
data2<-read.table("$outfile_diff", head=T, sep="\\t",row.names=1)
data.mc <- data.level1[,c(n)]
data.mc2 <- data2[,c(n)]
matrix.mc <- as.matrix(data.mc)
matrix.mc2 <- as.matrix(data.mc2)
pdf(file="$com_pdf")
heatmap(log2(matrix.mc), col=colors, scale="none", margins=c(2,9), cexCol=1.2)
graphics.off()\n
pdf(file="$diff_pdf")
heatmap(log2(matrix.mc2), col=colors, scale="none", margins=c(2,9), cexCol=1.2)
graphics.off()\n

