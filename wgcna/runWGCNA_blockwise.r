## load the wgcna package
library(WGCNA)

## this function redetects modules obtained from blockwise module detection
##  by branch cutting of the corresponding dendrograms

## read the command line argumnts
# 1: tpm, 2: blockwise base name, 3: CutHeight, 4: blockR object


args <- commandArgs(trailingOnly = TRUE)

blockwiseResDir=args[2]
cutHeight=as.numeric(args[3])
blockNetFile=args[4]


#outBase=args[4]
#output name
outName=paste("redetected_modules_",cutHeight,sep = "")

#working_dir=getwd()
#output_dir <- file.path(working_dir,outName)

if (!dir.exists(blockwiseResDir)){
stop("Blockwise module detection has not been run!")

} else {
    setwd(blockwiseResDir)
    dir.create(outName)
}

#read the tpm file
tpm=read.csv(args[1],check.names=FALSE)

## set first column of the file as the rownames 
rownames(tpm) <- tpm[,1]
tpm[,1] <- NULL

#rownames(tpm)
#colnames(tpm)

## Filter the data 
datExpr0 = as.data.frame(t(tpm))
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

datExpr = datExpr0

net=load(blockNetFile)

net2 <- recutBlockwiseTrees(datExpr, goodSamples=net$goodSamples, goodGenes =net$goodGenes,
                    blocks = net$blocks, TOMFiles = net$TOMFiles, dendrograms = net$dendrograms,
                    corType = "bicor", corOptions = "use = 'p', maxPOutliers = 0.05", networkType = "signed hybrid",
                    minModuleSize = 30, reassignThreshold = 0, 
                    mergeCutHeight = 0.10,  
                    detectCutHeight = 0.995, # changed from the default of 0.995
                    deepSplit = 2 , # changed from default of 2
                    numericLabels = TRUE,
                    verbose = 3)


net2ModuleColors <- labels2colors(net2$colors)
#get the module labels
net2ModuleLabels <-net2$colors

# look at how many genes per module
table(net2ModuleLabels)
write.csv(table(net2ModuleLabels), file="net2_modules_mergeCutHeight_0.10.csv")


save(net2, file="net2_network_mergeCutHeight0.10.RData")
names(net2)

# also save the components separately (as the manual suggests)
netMEs <- net2$MEs
netdendrograms <- net$dendrograms
# get the modules colours
netModuleColors <- labels2colors(net2$colors)

#get the module labels
netModuleLabels <- net2$colors

# 
save(netModuleColors, netModuleLabels, netMEs, netdendrograms, file = 
       "net2_network_components_mergeCutHeight0.10.RData")


### now want to output genes with module membership
gene_expr <- cbind(t(datExpr),netModuleLabels, netModuleColors)
head(gene_expr)

# write genes with modules to file
write.csv(gene_expr,"genes_with_modules_mergeCutHeight0.10.csv")



