#################################################################

##########################  ortholog #############################


rm(list = ls())
setwd("~/Desktop/PhD manuscript plan/gene expression/NewOrthologAnalysis/specificTE/")

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
nonRTELINEGal <- read.table("noRTE_ortho_gal_line")
colnames(nonRTELINEGal) <- "ID"
nonRTESINEGal <- read.table("noRTE_ortho_gal_sine")
colnames(nonRTESINEGal) <- "ID"
######
#
#    add new information to the big df
#
######
# if it is in give it a 1
TE_gene <- c("nonRTELINEGal", "nonRTESINEGal")

zero <- matrix(data = rep(0, dim(gal_expr.use)[1] * length(TE_gene)), nrow = dim(gal_expr.use)[1], ncol = length(TE_gene))
colnames(zero) <- TE_gene
gal.all.expr <- data.frame(gal_expr.use, zero)

for(i in TE_gene){
  te <- unique(get(i))
  gal.all.expr[as.character(gal.all.expr[,"ID"]) %in% as.character(te$ID),i] <- 1 	
}

colnames(gal.all.expr)[34:35] <- paste(colnames(gal.all.expr[,c(34:35)]), "IDs" , sep = "_")  

write.table(gal.all.expr, file = "gal_2nonRTEortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)

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

nonRTELINEAno <- read.table("noRTE_ortho_ano_line")
colnames(nonRTELINEAno) <- "ID"
nonRTESINEAno <- read.table("noRTE_ortho_ano_sine")
colnames(nonRTESINEAno) <- "ID"
nonRTEDNAAno <- read.table("noRTE_ortho_ano_dna")
colnames(nonRTEDNAAno) <- "ID"
nonRTEERVAno <- read.table("noRTE_ortho_ano_ervltr")
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

write.table(ano.all.expr, file = "ano_4nonRTEortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)


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

nonRTELINEBdg <- read.table("noRTE_ortho_bg_line")
colnames(nonRTELINEBdg) <- "ID"
nonRTESINEBdg <- read.table("noRTE_ortho_bg_sine")
colnames(nonRTESINEBdg) <- "ID"

######
#
#    add new information to the big df
#
######

# if it is in give it a 1
TE_gene <- c("nonRTELINEBdg", "nonRTESINEBdg", "nonRTEDNABdg")

zero <- matrix(data = rep(0, dim(bg_expr.use)[1] * length(TE_gene)), nrow = dim(bg_expr.use)[1], ncol = length(TE_gene))
colnames(zero) <- TE_gene
bg.all.expr <- data.frame(bg_expr.use, zero)


for(i in TE_gene){
  te <- unique(get(i))
  bg.all.expr[as.character(bg.all.expr[,"ID"]) %in% as.character(te$ID),i] <- 1 	
}

colnames(bg.all.expr)[17:18] <- paste(colnames(bg.all.expr[,c(17:18)]), "IDs" , sep = "_")                                                       
bg.all.expr$ID <- sub("^", "GeneID_", bg.all.expr$ID )

write.table(bg.all.expr, file = "bg_2nonRTEortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)


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
nonRTELINEOan <- read.table("noRTE_ortho_oana_line")
colnames(nonRTELINEOan) <- "ID"
nonRTESINEOan <- read.table("noRTE_ortho_oana_sine")
colnames(nonRTESINEOan) <- "ID"
######
#
#    add new information to the big df
#
######

# if it is in give it a 1
TE_gene <- c("nonRTELINEOan", "nonRTESINEOan")

zero <- matrix(data = rep(0, dim(oan_expr.use)[1] * length(TE_gene)), nrow = dim(oan_expr.use)[1], ncol = length(TE_gene))
colnames(zero) <- TE_gene
oan.all.expr <- data.frame(oan_expr.use, zero)


for(i in TE_gene){
  te <- unique(get(i))
  oan.all.expr[as.character(oan.all.expr[,"ID"]) %in% as.character(te$ID),i] <- 1 	
}

colnames(oan.all.expr)[35:36] <- paste(colnames(oan.all.expr[,c(35:36)]), "IDs" , sep = "_")                                                       

write.table(oan.all.expr, file = "oan_2nonRTEortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)


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
nonRTELINEMdo <- read.table("noRTE_ortho_mdo_line")
colnames(nonRTELINEMdo) <- "ID"
nonRTEERVMdo <- read.table("noRTE_ortho_mdo_ervltr")
colnames(nonRTEERVMdo) <- "ID"

######
#
#    add new information to the big df
#
######

# if it is in give it a 1
TE_gene <- c("nonRTELINEMdo", "nonRTEERVMdo")

zero <- matrix(data = rep(0, dim(mdo_expr.use)[1] * length(TE_gene)), nrow = dim(mdo_expr.use)[1], ncol = length(TE_gene))
colnames(zero) <- TE_gene
mdo.all.expr <- data.frame(mdo_expr.use, zero)


for(i in TE_gene){
  te <- unique(get(i))
  mdo.all.expr[as.character(mdo.all.expr[,"ID"]) %in% as.character(te$ID),i] <- 1 	
}

colnames(mdo.all.expr)[32:33] <- paste(colnames(mdo.all.expr[,c(32:33)]), "IDs" , sep = "_")                                                       

write.table(mdo.all.expr, file = "mdo_2nonRTEortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)

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
nonRTELINEHg <- read.table("noRTE_ortho_hg_line")
colnames(nonRTELINEHg) <- "ID"
nonRTESINEHg <- read.table("noRTE_ortho_hg_sine")
colnames(nonRTESINEHg) <- "ID"
nonRTEERVHg <- read.table("noRTE_ortho_hg_ervltr")
colnames(nonRTEERVHg) <- "ID"

######
#
#    add new information to the big df
#
######

# if it is in give it a 1
TE_gene <- c("nonRTELINEHg", "nonRTESINEHg","nonRTEERVHg")

zero <- matrix(data = rep(0, dim(hg_expr.use)[1] * length(TE_gene)), nrow = dim(hg_expr.use)[1], ncol = length(TE_gene))
colnames(zero) <- TE_gene
hg.all.expr <- data.frame(hg_expr.use, zero)


for(i in TE_gene){
  te <- unique(get(i))
  hg.all.expr[as.character(hg.all.expr[,"ID"]) %in% as.character(te$ID),i] <- 1 	
}

colnames(hg.all.expr)[31:33] <- paste(colnames(hg.all.expr[,c(31:33)]), "IDs" , sep = "_")                                                       

write.table(hg.all.expr, file = "hg_3nonRTEortholog.txt",quote = FALSE, sep = "\t", row.names=FALSE)
