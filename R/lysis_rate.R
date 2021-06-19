## lysis_rate.R
## An automated R script for estimating the taxon-specific lysis and lysis rate of microbe.
## Kevin Xu ZHONG, kevinzhong2006@gmail.com; xzhong@eoas.ubc.ca
## depends on R libraries reshape2 and VennDiagram being installed
## One input file is a tabular-separated otu table (.txt) contains one column named with "OTU_ID", another column named with "taxonomy", and the remaining columns named by "sample.id"_"one of three rDNA/rRNA fraction".  For example, if a water sample with "sample.id"="May_30m", then in this otu_table the column names for each rDNA/rRNA fraction of this sample is suggested to be: May_30m_Cellular.rDNA, May_30m_Cellular.rRNA and May_30m_Extracellular.rRNA.
## 
## lysis_rate.R function will firstly calculate the relative abundance of each fraction of the samples. Thus, it's important to be noted that the taxa within the otu_table should be pretreated to include only the microbial group the you are intended. For example, if you aim for calculate the taxon-specific lysis for prokaryotes, then the "Eukaryota", "Chloroplast" and "Mitochondria" should be removed from the otu table. Another example, if you aim for calculate the taxon-specific lysis for microeukaryotes then the "Metazoa", "Bacteria" and "Archaea" should be removed from the otu table
##
## The second input file is a tabular-separated metadata file (.txt) contains at least columns "sample.id", "Abundance_Total_Cellular.rRNA", "Abundance_Total_Extracellular.rRNA" and "turnover_rate".
## For example, "sample.id" = "May_30m", which is the name of water sample.
## "Abundance_Total_Cellular.rRNA" = Abundance of total cellular rRNA measured from the water sample.
## "Abundance_Total_Extracellular.rRNA" = Abundance of total extracellular rRNA measured from the water sample.
## "turnover_rate" = The estimated turnover rate of extracellular rRNA in the water sample
##


#' A lysis_rate Function
#'
#' This function allows calculating both the taxon-specific lysis and lysis rate of microbe, based on the input of a tabular-separated otu table (.txt) and metadata (.txt)
#' The otu_table contains data from three rDNA/rRNA fractions (Cellular rDNA, Cellular rRNA and Extracellular rRNA) of each water sample.
#' The metadata contains contains at least columns "sample.id", "Abundance_Total_Cellular.rRNA", "Abundance_Total_Extracellular.rRNA" and "turnover_rate".
#'
#'
#' @param otu_table 'Path/To/Your/OTU_TABLE_FILE'; For example: otu_table="/home/kevin/Desktop/data/10reads-plus/otu_table_qiime2_otu_16S.txt".
#' 
#' @param metadata 'Path/To/Your/METADATA_FILE'; For example: metadata="/home/kevin/Desktop/data/10reads-plus/metadata.txt".
#'
#' @keywords lysis_rate
#' @export
#'
#' @examples
#' ##### If you want to estimate both the taxon-specific lysis and lysis rate for prokaryotes, then run
#' lysis_rate(otu_table="Path/To/Your/OTU_TABLE_FILE_16S.txt", metadata="Path/To/Your/METADATA.txt")
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'


lysis_rate <- function(otu_table, metadata) {
  library(reshape2)
  library(VennDiagram)

  otu_table <- as.character(otu_table)#PATH to the otu_table
  metadata <- as.character(metadata)#PATH to the metadata
  
  if(is.null(otu_table)) {
    print("ERROR: otu_table is missing!")
  }
  if(is.null(metadata)) {
    print("ERROR: metadata is missing!")
  }
  else{
    print("The R is running, please be patient!")
    
    annA <- read.table(otu_table, h=T, sep="\t", quote=NULL, comment='', fill=T, stringsAsFactors = F)
    colnames(annA)
    colnames(annA)[colnames(annA)=="ConsensusLineage"] <- "taxonomy"
    annA$taxonomy[1]
    annA$OTU_ID[1]
    
    annA$Taxon <- paste(annA$OTU_ID, annA$taxonomy, sep = "___")
    annA <- annA[, !colnames(annA) %in% c("OTU_ID", "taxonomy")]
    
    rownames(annA) <- annA$Taxon
    annA <- annA[, !colnames(annA)=="Taxon"]
    for (i in 1:ncol(annA)){
      annA[,i] <- annA[,i]/sum(annA[,i])
    }#calculate the relative abundance
    annA$Taxon <- rownames(annA)
    
    #annA$Taxon <- gsub("; ", ";", annA$Taxon)
    #annA$Taxon <- gsub("__", "_", annA$Taxon)
    annAm <- melt(annA, id="Taxon")
    sum(annA[, 3])
    annA$Taxon[1]
    
    annm1 <- annAm
    annm1d <- dcast(annm1, Taxon ~ variable, value.var='value')
    annm1dm <- melt(annm1d, id="Taxon")
    annm1dm$value[is.na(annm1dm$value)] <- 0
    annm1dm[is.na(annm1dm$value), ]
    
    annm1dm1 <- annm1dm
    annm1dm1$variable <- gsub("_Cellular.rDNA", "_Cells", annm1dm1$variable)
    annm1dm2 <- annm1dm[grepl("_Cellular.rDNA",annm1dm$variable), ]
    
    ############
    #cpn <- read.table("16S_rDNA_q2-D2-bm2-400-ne1.txt", h=T, sep="\t", quote=NULL, comment='', fill=T, stringsAsFactors = F)
    #str(cpn)
    #min(cpn$num)#1
    #max(cpn$num)#15
    #min(cpn$copynum)#1
    #max(cpn$copynum)#15
    
    #cpn <- cpn[, c(1,2)]
    #colnames(cpn) <- c("Taxon", "copynum")
    
    #library(reshape2)
    #annm1dm2$Taxon1 <- sapply(strsplit(annm1dm2$Taxon, "___"), function(x) x[2])
    #annm1dm2 <- merge(annm1dm2, cpn, by.x="Taxon1", by.y="Taxon")
    #annm1dm2$copynum <- as.numeric(annm1dm2$copynum)
    #min(annm1dm2$copynum)#1
    #max(annm1dm2$copynum)#15
    
    annm1dm2$copynum <- 1###############remove this line if cpn is known
    
    
    ##########
    annm1dm2$value <- round(1000000 * annm1dm2$value / annm1dm2$copynum, digits = 10)
    str(annm1dm2)
    #sum(annm1dm2$value[annm1dm2$variable=="mp1"])
    
    annm1dm.new <- rbind(annm1dm1, annm1dm2[, c("Taxon", "variable", "value")])
    annm1dm.new[is.na(annm1dm.new$value), ]
    str(annm1dm.new)
    #sum(annm1dm.new$value[annm1dm.new$variable=="mp1"])
    
    ######calculate the relative abundance
    ann0 <- dcast(annm1dm.new, Taxon ~ variable, value.var='value')
    ann0m <- melt(ann0, id="Taxon")
    ann0m$value[is.na(ann0m$value)] <- 0
    ann <- dcast(ann0m, Taxon ~ variable, value.var='value')
    #sum(ann$mp1)
    str(ann)
    
    rownames(ann) <- ann$Taxon
    ann <- ann[, !colnames(ann)=="Taxon"]
    for (i in 1:ncol(ann)){
      ann[,i] <- ann[,i]/sum(ann[,i])
    }#calculate the relative abundance
    ann$Taxon <- rownames(ann)
    
    #############unique name of taxa
    tax <- as.data.frame(ann$Taxon)
    tax$count <- 1
    colnames(tax) <- c("idot", "count")
    tax$Taxon <- tax$idot
    tax$Taxon <- as.character(tax$Taxon)
    tax$Taxon <- sapply(strsplit(tax$Taxon, "___"), function(x) x[2])
    tax1 <- tax[!duplicated(tax$Taxon), ]
    
    #############
    ann4 <- ann
    row.names(ann4) <- ann4$Taxon
    ann4 <- ann4[, !colnames(ann4)=="Taxon"]
    ann4 <- ann4[!rowSums(ann4)==0, ]
    
    #to calculate the growth, relative_abundance, lysis
    pc <- ann4[, grepl("_Extracellular.rRNA", colnames(ann4))]
    pc$Taxon <- row.names(pc)
    mc <- ann4[, grepl("_Cellular.rRNA", colnames(ann4))]
    mc$Taxon <- row.names(mc)
    rmp <- ann4[, grepl("_Cellular.rDNA", colnames(ann4))]#raw.mp
    rmp$Taxon <- row.names(rmp)
    mp <- ann4[, !colnames(ann4) %in% c(colnames(pc), colnames(mc), colnames(rmp)), drop=F]
    mp$Taxon <- row.names(mp)
    
    library(reshape2)
    pcm <- melt(pc, id = 'Taxon')
    pcm$variable <- as.character(pcm$variable)
    pcm$sample.id <- sapply(strsplit(pcm$variable, "_Extracellular"), function(x) x[1])
    pcm$id <- paste(pcm$sample.id, pcm$Taxon, sep="___")
    
    mcm <- melt(mc, id = 'Taxon')
    mcm$variable <- as.character(mcm$variable)
    mcm$sample.id <- sapply(strsplit(mcm$variable, "_Cellular"), function(x) x[1])
    mcm$id <- paste(mcm$sample.id, mcm$Taxon, sep="___")
    
    mpm <- melt(mp, id = 'Taxon')
    mpm$variable <- as.character(mpm$variable)
    mpm$sample.id <- sapply(strsplit(mpm$variable, "_Cells"), function(x) x[1])
    mpm$id <- paste(mpm$sample.id, mcm$Taxon, sep="___")
    
    rmpm <- melt(rmp, id = 'Taxon')
    rmpm$variable <- as.character(rmpm$variable)
    rmpm$sample.id <- sapply(strsplit(rmpm$variable, "_Cellular"), function(x) x[1])
    rmpm$id <- paste(rmpm$sample.id, rmpm$Taxon, sep="___")
    
    meg1 <- merge(pcm, mcm, by.x="id", by.y="id")
    meg1$Lysis <- meg1$value.x / meg1$value.y
    colnames(meg1)[colnames(meg1)=="value.x"] <- "Relative.Abundance_Extracellular.rRNA"
    colnames(meg1)[colnames(meg1)=="value.y"] <- "Relative.Abundance_Cellular.rRNA"
    meg1 <- meg1[, c("id", "Lysis", "Relative.Abundance_Cellular.rRNA", "Relative.Abundance_Extracellular.rRNA")]
    
    meg2 <- merge(mcm, mpm, by.x="id", by.y="id")#to calculate growth
    meg2$Growth <- meg2$value.x / meg2$value.y
    meg2 <- meg2[, c("id", "Growth")]
    
    meg <- merge(mpm, meg1, by.x="id", by.y="id")
    meg <- merge(meg, meg2, by.x="id", by.y="id")
    colnames(meg)[colnames(meg)=="value"] <- "Relative.Abundance_Cells"
    meg$mp <- meg$Rlabundance
    
    rmpm2 <- rmpm[,c("id", "value")]
    meg <- merge(meg, rmpm2, by.x="id", by.y="id")
    colnames(meg)[colnames(meg)=="value"] <- "Relative.Abundance_Cellular.rDNA"
    colnames(meg)
    meg <- meg[, c("Taxon","sample.id","Relative.Abundance_Cells", "Relative.Abundance_Cellular.rDNA", "Relative.Abundance_Cellular.rRNA","Relative.Abundance_Extracellular.rRNA", "Lysis","Growth")]
    
    length(levels(factor(meg$Taxon)))#this show there are 3796 taxa for mc,mp,pc among 15 samples
    
    meg$Lysis[meg$Lysis=="Inf" | meg$Lysis=="-Inf" | is.na(meg$Lysis)] <- 0
    meg$Growth[meg$Growth=="Inf" | meg$Growth=="-Inf" | is.na(meg$Growth)] <- 0
    meg$group <- "no"
    meg$group[meg$Relative.Abundance_Cellular.rDNA==0 & !meg$Relative.Abundance_Cellular.rRNA==0 & meg$Relative.Abundance_Extracellular.rRNA==0] <- "I"
    meg$group[(!meg$Relative.Abundance_Cellular.rDNA==0) & (!meg$Relative.Abundance_Cellular.rRNA==0) & (meg$Relative.Abundance_Extracellular.rRNA==0)] <- "II"
    meg$group[(meg$Relative.Abundance_Cellular.rDNA==0) & (!meg$Relative.Abundance_Cellular.rRNA==0) & (!meg$Relative.Abundance_Extracellular.rRNA==0)] <- "III"
    meg$group[(!meg$Relative.Abundance_Cellular.rDNA==0) & (!meg$Relative.Abundance_Cellular.rRNA==0) & (!meg$Relative.Abundance_Extracellular.rRNA==0)] <- "IV"
    meg$group[(!meg$Relative.Abundance_Cellular.rDNA==0) & (meg$Relative.Abundance_Cellular.rRNA==0) & (meg$Relative.Abundance_Extracellular.rRNA==0)] <- "V"
    meg$group[(!meg$Relative.Abundance_Cellular.rDNA==0) & (meg$Relative.Abundance_Cellular.rRNA==0) & (!meg$Relative.Abundance_Extracellular.rRNA==0)] <- "VI"
    meg$group[(meg$Relative.Abundance_Cellular.rDNA==0) & (meg$Relative.Abundance_Cellular.rRNA==0) & (!meg$Relative.Abundance_Extracellular.rRNA==0)] <- "VII"
    
    levels(factor(meg$group))
    length(levels(factor(meg$mp)))
    length(levels(factor(meg$Relative.Abundance_Cellular.rDNA)))
    
    meg$Lysis[meg$group=="III"] <- 0
    
    str(meg)
    
    
    
    meg$mp <- meg$Relative.Abundance_Cells
    meg$raw.mp <- meg$Relative.Abundance_Cellular.rDNA
    meg$mc <- meg$Relative.Abundance_Cellular.rRNA
    meg$pc <- meg$Relative.Abundance_Extracellular.rRNA
    meg$Rlabundance <- meg$raw.mp
    
    
    
    ##########subset the taxa with extracellular ribosomes (pc) detected
    pcm2 <- pcm
    pcm2$ID <- pcm2$sample.id
    pcm3 <- pcm2[pcm2$value>0,]
    pcm4 <- dcast(pcm3, Taxon ~ ID, value.var='value')#1463 out of 2757 taxa has extracelluar ribosomes detected, in the datafram pcm4, the NA indicates the pc=0
    #library(xlsx)
    pcm5 <- pcm3[pcm3$Taxon %in% pcm4$Taxon, drop=F, ]
    pcm5d <- dcast(pcm5, Taxon ~ ID, value.var='value')#1463 out of 2757 taxa has extracelluar ribosomes detected, in the datafram pcm4, the NA indicates the pc=0
    #write.csv(pcm5d, paste0("taxa_with_extracellular_rRNA_detected_", as.character(nrow(pcm5d)), ".csv"), row.names = F, col.names = TRUE, sep="\t", quote = F) 
    #write.table(pcm5d, paste0("taxa_with_extracellular_rRNA_detected_", as.character(nrow(pcm5d)), ".txt"), row.names = F, col.names = TRUE, quote=F, sep="\t") 
    
    ###########focus on the lysis############################################################
    #subset these with measurable lysis in all 12 samples
    library(reshape2)
    data2 <- meg
    data2$ID <- data2$sample.id
    data <- dcast(data2, Taxon ~ ID, value.var='Lysis')#the original data for lysis for 1454 taxa
    per <- na.omit(data)#subset taxa that with lysis detected in all 12 samples, we can see there are 164 taxa, some with Inf 
    per1 <- melt(per, id = 'Taxon')
    per1 <- subset(per1, per1$value != "Inf")
    per1 <- subset(per1, per1$value != "-Inf")
    per1 <- subset(per1, per1$value>0)
    per2 <- dcast(per1, Taxon ~ variable, value.var='value')#these Inf has been replaced by NA, and thus we can remove as follow
    per3 <- na.omit(per2)#taxa with measurable lysis in all samples
    #write.csv(per3, paste0("taxa_with_measurable_lysis_in_all_samples_", as.character(nrow(per3)), ".csv"), row.names = F, col.names = TRUE, sep="\t", quote = F) 
    #write.table(per3, paste0("taxa_with_measurable_lysis_in_all_samples_", as.character(nrow(per3)), ".txt"), row.names = F, col.names = TRUE, sep="\t", quote = F) 
    
    per4m <- melt(per, id = 'Taxon')
    per4m <- subset(per4m, !per4m$value==0)
    per4md <- dcast(per4m, Taxon ~ variable, value.var='value')#
    per4 <- per4md
    #write.csv(per4,  paste0("taxa_with_measurable_lysis_in_either_sample_", as.character(nrow(per4)), ".csv"), row.names = F, col.names = TRUE, sep="\t", quote = F) 
    #write.table(per4,  paste0("taxa_with_measurable_lysis_in_either_sample_", as.character(nrow(per4)), ".txt"), row.names = F, col.names = TRUE, sep="\t", quote = F) 
    
    
    ############################# get a table with mp,mc,pc,lysis, growth
    #meg3 <- meg[, c("Taxon","sample.id","Relative.Abundance_Cells", "Relative.Abundance_Cellular.rDNA", "Relative.Abundance_Cellular.rRNA","Relative.Abundance_Extracellular.rRNA", "Lysis","Growth")]
    meg3 <- meg[, c("Taxon","sample.id", "Relative.Abundance_Cellular.rDNA", "Relative.Abundance_Cellular.rRNA","Relative.Abundance_Extracellular.rRNA", "Lysis")]
    
    meg3m <- melt(meg3, id = c("Taxon","sample.id"))
    meg3m$ID <- paste(as.character(meg3m$sample.id), as.character(meg3m$variable), sep = "_")
    meg4 <- dcast(meg3m, Taxon ~ ID, value.var='value')
    colnames(meg4)
    #write.csv(colnames(meg4), "colnames.csv", row.names = F, col.names = TRUE, sep="\t", quote = F) 
    #clist <- read.table("colnames01.txt", h=T, sep="\t")
    #clist$x <- as.character(clist$x)
    #meg4 <- meg4[, c(clist$x)]#reorder the column like clist
    meg4[is.na(meg4)] <- 0 #if not exist or undetectable, NA replaced by 0
    meg4 <- meg4[order(meg4$Taxon), ]
    
    #write.csv(meg4, paste0("output_", as.character(nrow(meg4)), "_taxa_alltogether.csv"), row.names = F, sep="\t", quote = F) 
    #write.table(meg4, paste0("output_", as.character(nrow(meg4)), "_taxa_alltogether.txt"), row.names = F, col.names = TRUE, sep="\t", quote = F) 
    
    ##############################venndiagram summary mp, mc, pc
    reset.df <- function(df) {
      row.names(df) <- df$Taxon
      df <- df[, !colnames(df)=="Taxon"]
      df1 <- subset(df, rowSums(df)!=0)
      df1$Taxon <- row.names(df1)
      df1 <- df1[, c("Taxon", colnames(df))]#
    }#a function to reset the dataframe
    
    mp.s <- reset.df(mp)
    mc.s <- reset.df(mc)
    pc.s <- reset.df(pc)
    
    df.data.mp <- as.data.frame(levels(factor(mp.s$Taxon)))
    colnames(df.data.mp) <- "Taxon"
    df.data.mp$type <- "Cellular_16S_rDNA"
    
    df.data.mc <- as.data.frame(levels(factor(mc.s$Taxon)))
    colnames(df.data.mc) <- "Taxon"
    df.data.mc$type <- "Cellular_16S_rRNA"
    
    df.data.pc <- as.data.frame(levels(factor(pc.s$Taxon)))
    colnames(df.data.pc) <- "Taxon"
    df.data.pc$type <- "Extracellular_16S_rRNA"
    
    summary.data <- rbind(df.data.pc, df.data.mc, df.data.mp)
    summary.data$count <- 1
    str(summary.data)
    
    
    library(VennDiagram)
    SNP_pop_11 <- summary.data$Taxon[summary.data$type=="Cellular_16S_rDNA"]
    SNP_pop_12 <- summary.data$Taxon[summary.data$type=="Cellular_16S_rRNA"]
    SNP_pop_13 <- summary.data$Taxon[summary.data$type=="Extracellular_16S_rRNA"]
    
    venn.diagram(
      x = list(SNP_pop_11 , SNP_pop_12, SNP_pop_13),
      category.names = c("Cellular_16S_rDNA" , "Cellular_16S_rRNA", "Extracellular_16S_rRNA"),
      filename = 'Venndiagram_summary_for_taxa_of_3_fractions_from_all_samples.png',
      output = TRUE ,
      imagetype="png" ,
      height = 480, 
      width = 480, 
      resolution = 300,
      compression = "lzw",
      lwd = 2,
      lty = 'blank',
      fill = c('purple', 'green', "grey"),
      cex = 1,
      fontface = "bold",
      fontfamily = "sans",
      cat.cex = 0.4,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.pos = c(27, 27, -27),
      cat.dist = c(-0.4, 0, -0.4),
      cat.fontfamily = "sans"
    )
    
    ###########lysis rate calculation
    abs <- read.table(metadata, h=T, sep="\t", quote=NULL, comment='', fill=T, stringsAsFactors = F)
    meg0 <- meg
    meg0[meg0=="Inf" | meg0=="-Inf" | is.na(meg0)] <- 0
    
    meg0 <- merge(meg0, abs, by.x="sample.id", by.y="sample.id")
    str(meg0)
    meg0$Lysis_rate <- (meg0$Lysis * meg0$turnover_rate * meg0$Abundance_Total_Extracellular.rRNA) / (meg0$Abundance_Total_Cellular.rRNA)
    
    meg0.lysis_rate <- dcast(meg0[!meg0$Lysis==0, ], Taxon~sample.id, value.var = "Lysis_rate")
    
    write.csv(meg0.lysis_rate,  paste0("taxa_with_measurable_lysis-rate_in_either_sample_", as.character(nrow(meg0.lysis_rate)), ".csv"), row.names = F, col.names = TRUE, sep="\t", quote = F) 
    write.table(meg0.lysis_rate,  paste0("taxa_with_measurable_lysis-rate_in_either_sample_", as.character(nrow(meg0.lysis_rate)), ".txt"), row.names = F, col.names = TRUE, sep="\t", quote = F) 
    
    
    #meg.export <- meg0[, c("Taxon","sample.id","Relative.Abundance_Cells", "Relative.Abundance_Cellular.rDNA", "Relative.Abundance_Cellular.rRNA","Relative.Abundance_Extracellular.rRNA", "Lysis","Growth")]
    meg.export <- meg0[, c("Taxon","sample.id","Relative.Abundance_Cellular.rDNA", "Relative.Abundance_Cellular.rRNA","Relative.Abundance_Extracellular.rRNA", "Lysis", "Lysis_rate", "group")]
    
    write.table(meg.export, file = "Output_Lysis&Lysis-rate&Relative-abundance.txt", row.names = F, col.names = TRUE, sep="\t", quote = F)
    #write.csv(meg.export, file = "Output_Lysis&Lysis-rate&Relative-abundance.csv", row.names = F, col.names = TRUE, sep="\t", quote = F)
    
  }
}
