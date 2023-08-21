
Convert('mye.h5ad', dest = "h5seurat", overwrite = F)
all <- LoadH5Seurat("mye.h5seurat",  meta.data = FALSE,misc = FALSE) ###meta.data and misc must be FALSE
all<- AddMetaData(all, mye_metadata, col.name = NULL)
Idents(all)='celltype2' ##whatever cluster name to be used.
DimPlot(mye, label = T)

library("scProportionTest")
mye<-subset(all, celltype2=='slan+ Monocytes' |  celltype2=='Migratory DCs' | celltype2=='MMP9+ MΦ' )

prop_test <- sc_utils(mye)
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype2",
  sample_1 = "Neonatal", sample_2 = "NEC_Neonatal",
  sample_identity = "Age_group"
)

permutation_plot(prop_test)
