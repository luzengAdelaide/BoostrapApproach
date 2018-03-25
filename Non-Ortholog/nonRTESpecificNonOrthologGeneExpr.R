#################################################################

##########################  non-ortholog #############################


rm(list = ls())
setwd("~/Desktop/PhD manuscript plan/gene expression/NewOrthologAnalysis/specificTE/")

########
#
#     reading chicken expression data
#
########

gal_expr <- read.delim("nonortho_gal_expr",sep="\t", header=T)
length_gal <- read.delim("nonortho_length_gal_TPM", sep="\t", header=T)
gal_expr_use <- log2(gal_expr[,2:ncol(gal_expr)])
gal_expr.use <- cbind(gal_expr$gal, length_gal$effective_length, gal_expr_use)
colnames(gal_expr.use)[1] = "ID"
colnames(gal_expr.use)[2] = "length"

########
#
#     reading and sorting of TE data
#
########

# TE data
nonRTELINEGal <- read.table("noRTE_nonortho_gal_line")
colnames(nonRTELINEGal) <- "ID"
nonRTESINEGal <- read.table("noRTE_nonortho_gal_sine")
colnames(nonRTESINEGal) <- "ID"
######
#
#    add new information to the big df
#
######

# if it is in give it a 1

TE_gene <- c("nonRTELINEGal", "nonRTESINEGal", "nonRTEERVGal")

zero <- matrix(data = rep(0, dim(gal_expr.use)[1] * length(TE_gene)), nrow = dim(gal_expr.use)[1], ncol = length(TE_gene))
colnames(zero) <- TE_gene
gal.all.expr <- data.frame(gal_expr.use, zero)


for(i in TE_gene){
  te <- unique(get(i))
  gal.all.expr[as.character(gal.all.expr[,"ID"]) %in% as.character(te$ID),i] <- 1 	
}

colnames(gal.all.expr)[34:35] <- paste(colnames(gal.all.expr[,c(34:56)]), "IDs" , sep = "_")                                                       

write.table(gal.all.expr, file = "gal_2nonRTENonortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)


#######################################################################################

########
#
#     reading Anolis expression data
#
########

ano_expr <- read.delim("nonortho_ano_expr",sep="\t", header=T)
length_ano <- read.delim("nonortho_length_ano_TPM", sep="\t", header=T)
ano_expr_use <- log2(ano_expr[,2:ncol(ano_expr)])
ano_expr.use <- cbind(ano_expr$ano, length_ano$effective_length, ano_expr_use)
colnames(ano_expr.use)[1] = "ID" 
colnames(ano_expr.use)[2] = "length"
########
#
#     reading and sorting of TE data
#
########

# TE data

nonRTELINEAno <- read.table("noRTE_nonortho_ano_line")
colnames(nonRTELINEAno) <- "ID"
nonRTESINEAno <- read.table("noRTE_nonortho_ano_sine")
colnames(nonRTESINEAno) <- "ID"
nonRTEDNAAno <- read.table("noRTE_nonortho_ano_dna")
colnames(nonRTEDNAAno) <- "ID"
nonRTEERVAno <- read.table("noRTE_nonortho_ano_ervltr")
colnames(nonRTEERVAno) <- "ID"

######
#
#    add new information to the big df
#
######

# if it is in give it a 1
TE_gene <- c("nonRTELINEAno", "nonRTESINEAno", "nonRTEDNAAno", "nonRTEERVAno")

zero <- matrix(data = rep(0, dim(ano_expr.use)[1] * length(TE_gene)), nrow = dim(ano_expr.use)[1], ncol = length(TE_gene))
colnames(zero) <- TE_gene
ano.all.expr <- data.frame(ano_expr.use, zero)


for(i in TE_gene){
  te <- unique(get(i))
  ano.all.expr[as.character(ano.all.expr[,"ID"]) %in% as.character(te$ID),i] <- 1 	
}

colnames(ano.all.expr)[33:36] <- paste(colnames(ano.all.expr[,c(33:36)]), "IDs" , sep = "_")                                                       

write.table(ano.all.expr, file = "ano_4nonRTENNonortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)


#######################################################################################

########
#
#     reading Pogona expression data
#
########

bg_expr <- read.delim("nonortho_bdg_expr",sep="\t", header=T)
length_bg <- read.delim("nonortho_length_bdg_TPM", sep="\t", header=T)
bg_expr_use <- log2(bg_expr[,2:ncol(bg_expr)])
bg_expr.use <- cbind(bg_expr$bdg, length_bg$effective_length, bg_expr_use)
colnames(bg_expr.use)[1] = "ID" 
colnames(bg_expr.use)[2] = "length"
########
#
#     reading and sorting of TE data
#
########

# TE data

nonRTELINEBdg <- read.table("noRTE_nonortho_bg_line")
colnames(nonRTELINEBdg) <- "ID"
nonRTESINEBdg <- read.table("noRTE_nonortho_bg_sine")
colnames(nonRTESINEBdg) <- "ID"
nonRTEDNABdg <- read.table("noRTE_nonortho_bg_dna")
colnames(nonRTEDNABdg) <- "ID"
nonRTEERVBdg <- read.table("noRTE_nonortho_bg_ervltr")
colnames(nonRTEERVBdg) <- "ID"

######
#
#    add new information to the big df
#
######

# if it is in give it a 1
TE_gene <- c("nonRTELINEBdg", "nonRTESINEBdg", "nonRTEDNABdg", "nonRTEERVBdg")

zero <- matrix(data = rep(0, dim(bg_expr.use)[1] * length(TE_gene)), nrow = dim(bg_expr.use)[1], ncol = length(TE_gene))
colnames(zero) <- TE_gene
bg.all.expr <- data.frame(bg_expr.use, zero)


for(i in TE_gene){
  te <- unique(get(i))
  bg.all.expr[as.character(bg.all.expr[,"ID"]) %in% as.character(te$ID),i] <- 1 	
}

colnames(bg.all.expr)[17:20] <- paste(colnames(bg.all.expr[,c(17:20)]), "IDs" , sep = "_")                                                       
bg.all.expr$ID <- sub("^", "GeneID_", bg.all.expr$ID )

write.table(bg.all.expr, file = "bg_4nonRTENonortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)


#######################################################################################

########
#
#     reading Platypus expression data
#
########

oan_expr <- read.delim("nonortho_oan_expr",sep="\t", header=T)
length_oan <- read.delim("nonortho_length_oan_TPM", sep="\t", header=T)
oan_expr_use <- log2(oan_expr[,2:ncol(oan_expr)])
oan_expr.use <- cbind(oan_expr$oan, length_oan$effective_length, oan_expr_use)
colnames(oan_expr.use)[1] = "ID"
colnames(oan_expr.use)[2] = "length"
########
#
#     reading and sorting of TE data
#
########

# TE data
nonRTELINEOan <- read.table("noRTE_nonortho_oana_line")
colnames(nonRTELINEOan) <- "ID"
nonRTESINEOan <- read.table("noRTE_nonortho_oana_sine")
colnames(nonRTESINEOan) <- "ID"
nonRTEERVOan <- read.table("noRTE_nonortho_oana_ervltr")
colnames(nonRTEERVOan) <- "ID"

######
#
#    add new information to the big df
#
######

# if it is in give it a 1
TE_gene <- c("nonRTELINEOan", "nonRTESINEOan", "nonRTEERVOan")

zero <- matrix(data = rep(0, dim(oan_expr.use)[1] * length(TE_gene)), nrow = dim(oan_expr.use)[1], ncol = length(TE_gene))
colnames(zero) <- TE_gene
oan.all.expr <- data.frame(oan_expr.use, zero)


for(i in TE_gene){
  te <- unique(get(i))
  oan.all.expr[as.character(oan.all.expr[,"ID"]) %in% as.character(te$ID),i] <- 1 	
}

colnames(oan.all.expr)[35:37] <- paste(colnames(oan.all.expr[,c(35:37)]), "IDs" , sep = "_")                                                       

write.table(oan.all.expr, file = "oan_3RTENonortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)


#######################################################################################

########
#
#     reading Opossum expression data
#
########

mdo_expr <- read.delim("nonortho_mdo_expr",sep="\t", header=T)
length_mdo <- read.delim("nonortho_length_mdo_TPM", sep="\t", header=T)
mdo_expr_use <- log2(mdo_expr[,2:ncol(mdo_expr)])
mdo_expr.use <- cbind(mdo_expr$mdo, length_mdo$effective_length, mdo_expr_use)
colnames(mdo_expr.use)[1] = "ID"
colnames(mdo_expr.use)[2] = "length"
########
#
#     reading and sorting of TE data
#
########

# TE data
nonRTELINEMdo <- read.table("noRTE_nonortho_mdo_line")
colnames(nonRTELINEMdo) <- "ID"
nonRTESINEMdo <- read.table("noRTE_nonortho_mdo_sine")
colnames(nonRTESINEMdo) <- "ID"
nonRTEDNAMdo <- read.table("noRTE_nonortho_mdo_dna")
colnames(nonRTEDNAMdo) <- "ID"
nonRTEERVMdo <- read.table("noRTE_nonortho_mdo_ervltr")
colnames(nonRTEERVMdo) <- "ID"

######
#
#    add new information to the big df
#
######

# if it is in give it a 1
TE_gene <- c("nonRTELINEMdo", "nonRTESINEMdo", "nonRTEDNAMdo", "nonRTEERVMdo")

zero <- matrix(data = rep(0, dim(mdo_expr.use)[1] * length(TE_gene)), nrow = dim(mdo_expr.use)[1], ncol = length(TE_gene))
colnames(zero) <- TE_gene
mdo.all.expr <- data.frame(mdo_expr.use, zero)


for(i in TE_gene){
  te <- unique(get(i))
  mdo.all.expr[as.character(mdo.all.expr[,"ID"]) %in% as.character(te$ID),i] <- 1 	
}

colnames(mdo.all.expr)[32:35] <- paste(colnames(mdo.all.expr[,c(32:35)]), "IDs" , sep = "_")                                                       

write.table(mdo.all.expr, file = "mdo_4nonRTENonortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)

#######################################################################################

########
#
#     reading Human expression data
#
########

hg_expr <- read.delim("nonortho_hgs_expr",sep="\t", header=T)
length_hg <- read.delim("nonortho_length_hgs_TPM",sep="\t", header=T)
hg_expr_use <- log2(hg_expr[,2:ncol(hg_expr)])
hg_expr.use <- cbind(hg_expr$hgs, length_hg$effective_length, hg_expr_use)
colnames(hg_expr.use)[1] = "ID"
colnames(hg_expr.use)[2] = "length"
########
#
#     reading and sorting of TE data
#
########

# TE data
nonRTELINEHg <- read.table("noRTE_nonortho_hg_line")
colnames(nonRTELINEHg) <- "ID"
nonRTESINEHg <- read.table("noRTE_nonortho_hg_sine")
colnames(nonRTESINEHg) <- "ID"
nonRTEERVHg <- read.table("noRTE_nonortho_hg_dna")
colnames(nonRTEERVHg) <- "ID"
nonRTEDNAHg <- read.table("noRTE_nonortho_hg_ervltr")
colnames(nonRTEDNAHg) <- "ID"
######
#
#    add new information to the big df
#
######

# if it is in give it a 1
TE_gene <- c("nonRTELINEHg", "nonRTESINEHg","nonRTEERVHg", "nonRTEDNAHg")

zero <- matrix(data = rep(0, dim(hg_expr.use)[1] * length(TE_gene)), nrow = dim(hg_expr.use)[1], ncol = length(TE_gene))
colnames(zero) <- TE_gene
hg.all.expr <- data.frame(hg_expr.use, zero)


for(i in TE_gene){
  te <- unique(get(i))
  hg.all.expr[as.character(hg.all.expr[,"ID"]) %in% as.character(te$ID),i] <- 1 	
}

colnames(hg.all.expr)[31:34] <- paste(colnames(hg.all.expr[,c(31:34)]), "IDs" , sep = "_")                                                       

write.table(hg.all.expr, file = "hg_4nonRTENonortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)

