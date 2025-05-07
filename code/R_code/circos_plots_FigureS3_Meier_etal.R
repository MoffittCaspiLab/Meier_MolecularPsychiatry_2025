##############################################################################################################
## Scripts to create Circos plots for Meier et al., "DNA Methylation Profiles of Long-Term Cannabis Users in Midlife: A Comprehensive Evaluation of Published Cannabis-Associated Methylation Markers in a Representative Cohort"; Figure S3
##############################################################################################################

set.seed(2024)

sessionInfo()

## Load necessary libraries

# install.packages("circlize")  # Install circlize
library(circlize)             # Load the library
library(RColorBrewer)
library(dplyr)


#### input genomic location data

methyl.df <- read.table("./ExpressionvsMeth/Data/annotated_probeMap_annotated.txt", sep="\t", header=F, stringsAsFactors=F, quote="")
exp.df <- read.table("./ExpressionvsMeth/Data/Primeview_updatedAnnotation.txt", sep="\t", header=T, stringsAsFactors=F, quote="")
Locations.df <- data.frame(
  chromosome=factor(methyl.df$V4, levels=c('1','2','3','4','5','6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')),
  position = apply(methyl.df[,5:6], 1, function(x) { x<- as.numeric(x); round(mean(c(x[1], x[2])))})
)
rownames(Locations.df) <- methyl.df$V1

expLocations.df <- data.frame(
  chromosome=factor(exp.df$Chromosome, levels=c('1','2','3','4','5','6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')),
  position=apply(exp.df[,8:9], 1, function(x) {
    if( is.numeric(x) ) { round(mean(c(x[1], x[2]))) } else { NA }
  })
)
rownames(expLocations.df) <- exp.df$ID

## remove expLocation rows where chromosome is NA
expLocations.df <- filter(expLocations.df, !is.na(chromosome), !is.na(position))

###############################################################################################################
### get the names of the dataframes containing rho and p

cgdfs <- ls(pattern = "^cg\\d+_combined") # Matches names like cgNNN_combined
geneid <- gsub("_combined", "", cgdfs[24])

fname <- paste0("./ExpressionvsMeth/Results/", geneid, "_combined.txt")

data.df <- read.table(fname, sep="\t", skip=1, header=F, stringsAsFactors=F)
data.df <- data.df[which(data.df[,2] <= 0.05/485577 & data.df[,4] <= 0),] # p-values below Bonferroni cut-off

## Identify the nearby sites

geneChr <- as.character(methyl.df$V4[which(methyl.df$V1==geneid)[1]])
geneStart <- as.numeric(methyl.df$V5[which(methyl.df$V1==geneid)[1]])
geneEnd <- as.numeric(methyl.df$V6[which(methyl.df$V1==geneid)[1]])
geneMid <- round(mean(c(geneStart, geneEnd)))
geneName <- as.character(methyl.df$v8[which(methyl.df$V1==geneid)[1]])

## link the significant relationships

## create basic circos plot

# Initialize the circos plot
color_palette <- brewer.pal(8, name = "Set3")
color_palette <- colorRampPalette(color_palette)(24) 
par(new=TRUE)
circos.par("track.height" = 0.1, "canvas.xlim" = c(-2, 2), "canvas.ylim" = c(-2, 2), "start.degree"=360)
circos.initialize(factors = expLocations.df$chromosome, x = expLocations.df$position)
circos.trackPlotRegion(
  factors = expLocations.df$chromosome,
  ylim=c(0,1),
  panel.fun = function(x, y) {
    circos.axis(labels=FALSE, major.tick=FALSE)
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + 1.5, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "black")
    # gene <- gene_names[which(expLocations.df$chromosome == sector.name)]
    # circos.text(CELL_META$xcenter, 0.5, targets, cex = 0.6, col = "black", facing = "bending")
  },
  bg.border = '#999999', 
  bg.col = color_palette, 
  bg.lwd=1
)

## add links between genes and probes

for( i in 1:nrow(data.df) ) {
  circos.link(
    geneChr, 
    geneMid, 
    as.character(expLocations.df[data.df[i,1], 1]), 
    as.numeric(expLocations.df[data.df[i,1], 2]),
    col= "black",
    directional=-1
  )
}

title(paste0("Genomic locations of gene expression correlates of ", geneid))
print(geneid)
circos.clear()

dev.off()

