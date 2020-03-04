library(dplyr)
library(reshape2)

# BP
pop.filepath = "data/external/Enrichr/GO_Biological_Process_2018.txt"
num.col = max(count.fields(pop.filepath, sep = "\t", blank.lines.skip = TRUE), na.rm = TRUE)
GO_BP = read.table(pop.filepath, sep = "\t", header = FALSE, stringsAsFactors = FALSE, fill = TRUE, col.names = 1:num.col, blank.lines.skip = TRUE, quote = "")
rownames(GO_BP) = GO_BP[,1]
GO_BP = GO_BP[,-c(1:2)]
GO_BP.df = melt(as.matrix(GO_BP), na.rm = T) %>% select(-Var2) %>% filter(value!="")
# MF
pop.filepath = "data/external/Enrichr/GO_Molecular_Function_2018.txt"
num.col = max(count.fields(pop.filepath, sep = "\t", blank.lines.skip = TRUE), na.rm = TRUE)
GO_MF = read.table(pop.filepath, sep = "\t", header = FALSE, stringsAsFactors = FALSE, fill = TRUE, col.names = 1:num.col, blank.lines.skip = TRUE, quote = "")
rownames(GO_MF) = GO_MF[,1]
GO_MF = GO_MF[,-c(1:2)]
GO_MF.df = melt(as.matrix(GO_MF), na.rm = T) %>% select(-Var2) %>% filter(value!="")

# CC
pop.filepath = "data/external/Enrichr/GO_Cellular_Component_2018.txt"
num.col = max(count.fields(pop.filepath, sep = "\t", blank.lines.skip = TRUE), na.rm = TRUE)
GO_CC = read.table(pop.filepath, sep = "\t", header = FALSE, stringsAsFactors = FALSE, fill = TRUE, col.names = 1:num.col, blank.lines.skip = TRUE, quote = "")
rownames(GO_CC) = GO_CC[,1]
GO_CC = GO_CC[,-c(1:2)]
GO_CC.df = melt(as.matrix(GO_CC), na.rm = T) %>% select(-Var2) %>% filter(value!="")


combined.df = rbind(GO_BP.df,GO_MF.df)
combined.df = rbind(combined.df, GO_CC.df)

# find distinct pairs now
combined.df.2 = combined.df[!duplicated(t(apply(combined.df, 1, sort))),]
# report distinct genes(proteins)
# length(unique(combined.df.2[,2]))
# [1] 15714

# report distinct [genes(proteins) -- GO terms] associations
# nrow(combined.df.2)
# [1] 250734
