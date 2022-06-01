#load the wgcna package
library(WGCNA)

## read the command line argumnts
# 1: file path, 
# 2: Soft threshold to use, 
# 3: CutHeight, 
# 4: Output base name, 
# 5: Output directory
# 6: Samples to run WGCNA (if our data set contains many samples and if we want to
# run WGCNA only on subset of samples pass sampes that we use as text file
# Where column contains samples that are need to select
# otherwise use set this argument to "all" it will use all samples)


args <- commandArgs(trailingOnly = TRUE)

softThresh=as.numeric(args[2])
cutHeight=as.numeric(args[3])
outBase=args[4]
out=args[5]

# Create the name for the output base directory
outName=paste(outBase,"_wgcna_",cutHeight,"_",softThresh,sep = "")

# Create the location to save the files
working_dir=getwd()
output_dir <- file.path(working_dir,outName)

# Check if the directory is already there and if it exists stop! 
if (!dir.exists(output_dir)){
dir.create(output_dir)
} else {
    stop("Output directory already exists!")
}

# set working directory 
setwd(output_dir)

#read the tpm file from the first command line argument (expression file)
# sample X gene matrix
# rows are samples and columns are genes

tpm=read.csv(args[1],check.names=FALSE)

## set first column of the file as the rownames 
rownames(tpm) <- tpm[,1]
tpm[,1] <- NULL

#rownames(tpm)
#colnames(tpm)

if(args[6]=="all") {
      datExpr0 = data.matrix(tpm)

} else {
   subset_sampes = read.csv(args[6], check.names = FALSE, header = TRUE)
   tpm <- subset(tpm, rownames(tpm) %in%  subset_sampes[,1])
   
}

## Filter the data using WGCNA filering function


datExpr0 = data.matrix(tpm)
gsg = goodSamplesGenes(datExpr0, verbose = 3)

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

datExpr=datExpr0

net = blockwiseModules(  datExpr, 
                         maxBlockSize = 30000, 
                         power = softThresh, 
                         networkType = "signed hybrid", 
                         TOMType = "unsigned", 
                         minModuleSize = 30,
                         corType="bicor", 
                         corOptions = "use = 'p', maxPOutliers = 0.05", 
                         reassignThreshold = 0,
                         replaceMissingAdjacencies = TRUE, 
                         mergeCutHeight = cutHeight,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "Signed_hybrid_TOM-blockwise_maxP0.05",
                         verbose = 3)

#First let's see what the bwnet dataset contains
names(net)
length(table(net$blocks))
table(net$blocks)


RobjFileName=paste(output_dir,"/","Consensus-NetworkConstruction-cutHight",cutHeight,".RDS",sep="")
GeneExpFileName=paste(output_dir,"/","Exp_with_moduleGenes",cutHeight,".csv",sep="")
modLabFileName=paste(output_dir,"/","moduleLabels",cutHeight,".csv",sep="")
netWokrCompFileName=paste(output_dir,"/","network_components",cutHeight,"_RData",sep="")

# save the R object
saveRDS(net, file=RobjFileName)

netModuleColors <- labels2colors(net$colors)

#get the module labels
netModuleLabels <- net$colors

# look at how many genes per module
table(netModuleLabels)
write.csv(table(netModuleLabels), file=modLabFileName)

netMEs <- net$MEs
netdendrograms <- net$dendrograms

save(netModuleColors, netModuleLabels, netMEs, netdendrograms, file = netWokrCompFileName)

### now want to output genes with module membership
gene_expr <- cbind(t(datExpr),netModuleLabels, netModuleColors)

#head(gene_expr)

# write genes with modules to file
write.csv(gene_expr,file=GeneExpFileName)
