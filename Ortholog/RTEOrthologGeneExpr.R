#################################################################

##########################  ortholog #############################


rm(list = ls())
setwd("~/Desktop/gene_expression/OrthologAnalysis/")

########
#
#     reading chicken expression data
#
########

gal_expr <- read.delim("total_gal",sep="\t", header=T)
length_gal <- read.delim("length_ortho_gal", sep="\t", header=F)
gal_expr_use <- log2(gal_expr[,2:ncol(gal_expr)])
gal_expr.use <- cbind(gal_expr$gal, length_gal$V2, gal_expr_use)
colnames(gal_expr.use)[1] = "ID"
colnames(gal_expr.use)[2] = "length"

########
#
#     reading and sorting of TE data
#
########


# TE data
recentLINEGal <- read.table("recent_gal_line_ortho")
colnames(recentLINEGal) <- "ID"
######
#
#    add new information to the big df
#
######
# if it is in give it a 1
zero <- matrix(data = rep(0, dim(gal_expr.use)[1]), nrow = dim(gal_expr.use)[1], ncol = 1)
colnames(zero) <- "recentLINEGal"
gal.all.expr <- data.frame(gal_expr.use, zero)

gal.all.expr[as.character(gal.all.expr$ID) %in% as.character(recentLINEGal$ID), 34] <- 1 	

colnames(gal.all.expr)[34] <- "recentLINEGal_IDs"                                                      

write.table(gal.all.expr, file = "gal_1RTEortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)

#######################################################################################

########
#
#     reading Anolis expression data
#
########

ano_expr <- read.delim("total_ano",sep="\t", header=T)
length_ano <- read.delim("length_ortho_ano", sep="\t", header=F)
ano_expr_use <- log2(ano_expr[,2:ncol(ano_expr)])
ano_expr.use <- cbind(ano_expr$ano, length_ano$V2, ano_expr_use)
colnames(ano_expr.use)[1] = "ID" 
colnames(ano_expr.use)[2] = "length"
########
#
#     reading and sorting of TE data
#
########

# TE data

recentLINEAno <- read.table("recent_ano_line_ortho")
colnames(recentLINEAno) <- "ID"
recentSINEAno <- read.table("recent_ano_sine_ortho")
colnames(recentSINEAno) <- "ID"
recentDNAAno <- read.table("recent_ano_dna_ortho")
colnames(recentDNAAno) <- "ID"
recentERVAno <- read.table("recent_ano_ervltr_ortho")
colnames(recentERVAno) <- "ID"

######
#
#    add new information to the big df
#
######

# if it is in give it a 1
TE_gene <- c("recentLINEAno", "recentSINEAno", "recentDNAAno", "recentERVAno")

zero <- matrix(data = rep(0, dim(ano_expr.use)[1] * length(TE_gene)), nrow = dim(ano_expr.use)[1], ncol = length(TE_gene))
colnames(zero) <- TE_gene
ano.all.expr <- data.frame(ano_expr.use, zero)


for(i in TE_gene){
  te <- unique(get(i))
  ano.all.expr[as.character(ano.all.expr[,"ID"]) %in% as.character(te$ID),i] <- 1 	
}

colnames(ano.all.expr)[33:36] <- paste(colnames(ano.all.expr[,c(33:36)]), "IDs" , sep = "_")                                                       

write.table(ano.all.expr, file = "ano_4RTEortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)


#######################################################################################

########
#
#     reading Pogona expression data
#
########

bg_expr <- read.delim("total_bdg",sep="\t", header=T)
length_bg <- read.delim("length_ortho_bdg", sep="\t", header=F)
bg_expr_use <- log2(bg_expr[,2:ncol(bg_expr)])
bg_expr.use <- cbind(bg_expr$bdg, length_bg$V2, bg_expr_use)
colnames(bg_expr.use)[1] = "ID" 
colnames(bg_expr.use)[2] = "length"
########
#
#     reading and sorting of TE data
#
########

# TE data

recentLINEBdg <- read.table("recent_bg_line_ortho")
colnames(recentLINEBdg) <- "ID"
recentSINEBdg <- read.table("recent_bg_sine_ortho")
colnames(recentSINEBdg) <- "ID"
recentDNABdg <- read.table("recent_bg_dna_ortho")
colnames(recentDNABdg) <- "ID"

######
#
#    add new information to the big df
#
######

# if it is in give it a 1
TE_gene <- c("recentLINEBdg", "recentSINEBdg", "recentDNABdg")

zero <- matrix(data = rep(0, dim(bg_expr.use)[1] * length(TE_gene)), nrow = dim(bg_expr.use)[1], ncol = length(TE_gene))
colnames(zero) <- TE_gene
bg.all.expr <- data.frame(bg_expr.use, zero)


for(i in TE_gene){
  te <- unique(get(i))
  bg.all.expr[as.character(bg.all.expr[,"ID"]) %in% as.character(te$ID),i] <- 1 	
}

colnames(bg.all.expr)[17:19] <- paste(colnames(bg.all.expr[,c(17:19)]), "IDs" , sep = "_")                                                       
bg.all.expr$ID <- sub("^", "GeneID_", bg.all.expr$ID )

write.table(bg.all.expr, file = "bg_3RTEortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)


#######################################################################################

########
#
#     reading Platypus expression data
#
########

oan_expr <- read.delim("total_oan",sep="\t", header=T)
length_oan <- read.delim("length_ortho_oan", sep="\t", header=F)
oan_expr_use <- log2(oan_expr[,2:ncol(oan_expr)])
oan_expr.use <- cbind(oan_expr$oan, length_oan$V2, oan_expr_use)
colnames(oan_expr.use)[1] = "ID"
colnames(oan_expr.use)[2] = "length"
########
#
#     reading and sorting of TE data
#
########

# TE data
recentLINEOan <- read.table("recent_oana_line_ortho")
colnames(recentLINEOan) <- "ID"
recentSINEOan <- read.table("recent_oana_sine_ortho")
colnames(recentSINEOan) <- "ID"

######
#
#    add new information to the big df
#
######

# if it is in give it a 1
TE_gene <- c("recentLINEOan", "recentSINEOan")

zero <- matrix(data = rep(0, dim(oan_expr.use)[1] * length(TE_gene)), nrow = dim(oan_expr.use)[1], ncol = length(TE_gene))
colnames(zero) <- TE_gene
oan.all.expr <- data.frame(oan_expr.use, zero)


for(i in TE_gene){
  te <- unique(get(i))
  oan.all.expr[as.character(oan.all.expr[,"ID"]) %in% as.character(te$ID),i] <- 1 	
}

colnames(oan.all.expr)[35:36] <- paste(colnames(oan.all.expr[,c(35:36)]), "IDs" , sep = "_")                                                       

write.table(oan.all.expr, file = "oan_2RTEortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)


#######################################################################################

########
#
#     reading Opossum expression data
#
########

mdo_expr <- read.delim("total_mdo",sep="\t", header=T)
length_mdo <- read.delim("length_ortho_mdo", sep="\t", header=F)
mdo_expr_use <- log2(mdo_expr[,2:ncol(mdo_expr)])
mdo_expr.use <- cbind(mdo_expr$mdo, length_mdo$V2, mdo_expr_use)
colnames(mdo_expr.use)[1] = "ID"
colnames(mdo_expr.use)[2] = "length"
########
#
#     reading and sorting of TE data
#
########

# TE data
recentLINEMdo <- read.table("recent_mdo_line_ortho")
colnames(recentLINEMdo) <- "ID"
recentSINEMdo <- read.table("recent_mdo_sine_ortho")
colnames(recentSINEMdo) <- "ID"
recentDNAMdo <- read.table("recent_mdo_dna_ortho")
colnames(recentDNAMdo) <- "ID"
recentERVMdo <- read.table("recent_mdo_ervltr_ortho")
colnames(recentERVMdo) <- "ID"

######
#
#    add new information to the big df
#
######

# if it is in give it a 1
TE_gene <- c("recentLINEMdo", "recentSINEMdo", "recentDNAMdo", "recentERVMdo")

zero <- matrix(data = rep(0, dim(mdo_expr.use)[1] * length(TE_gene)), nrow = dim(mdo_expr.use)[1], ncol = length(TE_gene))
colnames(zero) <- TE_gene
mdo.all.expr <- data.frame(mdo_expr.use, zero)


for(i in TE_gene){
  te <- unique(get(i))
  mdo.all.expr[as.character(mdo.all.expr[,"ID"]) %in% as.character(te$ID),i] <- 1 	
}

colnames(mdo.all.expr)[32:35] <- paste(colnames(mdo.all.expr[,c(32:35)]), "IDs" , sep = "_")                                                       

write.table(mdo.all.expr, file = "mdo_4RTEortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)

#######################################################################################

########
#
#     reading Human expression data
#
########

hg_expr <- read.delim("total_hgs",sep="\t", header=T)
length_hg <- read.delim("length_ortho_hgs",sep="\t", header=F)
hg_expr_use <- log2(hg_expr[,2:ncol(hg_expr)])
hg_expr.use <- cbind(hg_expr$hgs, length_hg$V2, hg_expr_use)
colnames(hg_expr.use)[1] = "ID"
colnames(hg_expr.use)[2] = "length"
########
#
#     reading and sorting of TE data
#
########

# TE data
recentLINEHg <- read.table("recent_hg_line_ortho")
colnames(recentLINEHg) <- "ID"
recentSINEHg <- read.table("recent_hg_sine_ortho")
colnames(recentSINEHg) <- "ID"
recentERVHg <- read.table("recent_hg_ervltr_ortho")
colnames(recentERVHg) <- "ID"

######
#
#    add new information to the big df
#
######

# if it is in give it a 1
TE_gene <- c("recentLINEHg", "recentSINEHg","recentERVHg")

zero <- matrix(data = rep(0, dim(hg_expr.use)[1] * length(TE_gene)), nrow = dim(hg_expr.use)[1], ncol = length(TE_gene))
colnames(zero) <- TE_gene
hg.all.expr <- data.frame(hg_expr.use, zero)


for(i in TE_gene){
  te <- unique(get(i))
  hg.all.expr[as.character(hg.all.expr[,"ID"]) %in% as.character(te$ID),i] <- 1 	
}

colnames(hg.all.expr)[31:33] <- paste(colnames(hg.all.expr[,c(31:33)]), "IDs" , sep = "_")                                                       

write.table(hg.all.expr, file = "hg_3RTEortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)
