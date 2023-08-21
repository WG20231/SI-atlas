library(CellChat)
library(stringr)

#####fetal_firsttrim
f1_matrix<- read.csv("/gpfs/ysm/project/konnikova/wg297/final2/remove_repeat/fetal_normatrix_011223.csv", row.names=1)
f1_matrix=as.matrix(f1_matrix)

f1_metadata <- read.csv("/gpfs/ysm/project/konnikova/wg297/final2/remove_repeat/fetal_metadata_011223.tsv", row.names=1)
rownames(f1_metadata)<-str_replace(rownames(f1_metadata),"-", ".")
rownames(f1_metadata)<-str_replace(rownames(f1_metadata),"-", ".")
f1_metadata=f1_metadata %>% select('celltype2')
cellchat_f1 <- createCellChat(object = f1_matrix, meta = f1_metadata, group.by = "celltype2")
saveRDS(cellchat_f1,'cellchat_fetal_firsttrim.rds')

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = c( "Secreted Signaling",'Cell-Cell contact')) # use Secreted Signaling
cellchat_f1@DB <- CellChatDB.use
cellchat_f1 <- subsetData(cellchat_f1)
cellchat_f1 <- identifyOverExpressedGenes(cellchat_f1)
cellchat_f1 <- identifyOverExpressedInteractions(cellchat_f1)
cellchat_f1 <- computeCommunProb(cellchat_f1)
cellchat_f1 <- filterCommunication(cellchat_f1, min.cells = 10)
cellchat_f1 <- computeCommunProbPathway(cellchat_f1) 
cellchat_f1 <- aggregateNet(cellchat_f1)
cellchat_f1 <- netAnalysis_computeCentrality(cellchat_f1, slot.name = "netP")
cellchat_f1@netP$prob
saveRDS(cellchat_f1,'cellchat_fetal.rds')

#####postnatal
f2_matrix<- read.csv("/gpfs/ysm/project/konnikova/wg297/final2/remove_repeat/postnatal_normatrix_011223.csv", row.names=1)
f2_matrix=as.matrix(f2_matrix)

f2_metadata <- read.csv("/gpfs/ysm/project/konnikova/wg297/final2/remove_repeat/postnatal_metadata_011223.tsv", row.names=1)
rownames(f2_metadata)<-str_replace(rownames(f2_metadata),"-", ".")
rownames(f2_metadata)<-str_replace(rownames(f2_metadata),"-", ".")
f2_metadata=f2_metadata %>% select('celltype2')
cellchat_f2 <- createCellChat(object = f2_matrix, meta = f2_metadata, group.by = "celltype2")
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = c( "Secreted Signaling",'Cell-Cell contact')) # use Secreted Signaling
cellchat_f2@DB <- CellChatDB.use
cellchat_f2 <- subsetData(cellchat_f2)
cellchat_f2 <- identifyOverExpressedGenes(cellchat_f2)
cellchat_f2 <- identifyOverExpressedInteractions(cellchat_f2)
cellchat_f2 <- computeCommunProb(cellchat_f2)
cellchat_f2 <- filterCommunication(cellchat_f2, min.cells = 10)
cellchat_f2 <- computeCommunProbPathway(cellchat_f2) 
cellchat_f2 <- aggregateNet(cellchat_f2)
cellchat_f2 <- netAnalysis_computeCentrality(cellchat_f2, slot.name = "netP")
saveRDS(cellchat_f2,'cellchat_postnatal.rds')

object <- list(Fetal = cellchat_f1, Postnatal = cellchat_f2)
cellchat <- mergeCellChat(object, add.names = names(object))
gg1 <- compareInteractions(cellchat, show.legend = F)
gg2 <- compareInteractions(cellchat, show.legend = F, measure = "weight")
gg1 + gg2

netVisual_heatmap(cellchat, measure = "weight",comparison = c(2, 1))

####differential analysis 
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Fetal"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat2, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat2, net = net, datasets = "Fetal",ligand.logFC = 0.1, receptor.logFC = 0.1)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
#net.down <- subsetCommunication(cellchat2, net = net, datasets = "NL",ligand.logFC = -0.1, receptor.logFC = -0.1)
install.packages("wordcloud")
net.up<-drop_na(net.up)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat2)
computeEnrichmentScore(net.up, species = 'human')
netVisual_chord_gene(object[[1]], sources.use = NULL, targets.use = NULL, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object)[1]),legend.pos.x = 5,legend.pos.y = 80)


