#########################
##### Tissue Damage #####
########################
library(ggpubr)
library(scater)
library(edgeR)
library(dplyr)
library(xlsx)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Ss.eg.db)

##### DEG through Seurat pipeline ######
# PT PT_Prolif PT_VIM+ TAL_1 TAL_2 DTC IC_TypeA IC_TypeB
PTC <- subset(kidneys, subset = cell_types %in% c("PT","PT_Prolif","PT_VIM+"))
PTC_1 <- subset(PTC, Condition %in% c("kidney_1st_treated","kidney_1st_untreated_rt"))
PTC_2 <- subset(PTC, Condition %in% c("kidney_2nd_treated","kidney_2nd_untreated_rt"))

# Find Differentially Expressed genes with Seurat FindMarkers
Idents(PTC_1) <- PTC_1$Condition
deg.1.df <- FindMarkers(PTC_1,
								assay = "RNA",
								ident.1 = "kidney_1st_treated") %>% arrange(desc(avg_log2FC))

Idents(PTC_2) <- PTC_2$Condition

deg.2.df <- FindMarkers(PTC_2,
								assay = "RNA",
								ident.1 = "kidney_2nd_treated") %>% arrange(desc(avg_log2FC))
# set thresholds 
fold_cutoff = 0.5
pvalue_cutoff = 0.01

plot.deg1 <- deg.1.df %>% 
	tibble::rownames_to_column("geneIDs") %>% 
	dplyr::select(geneIDs, avg_log2FC, p_val) %>% 
	mutate(logpval = -log10(p_val))
colnames(plot.deg1) <- c("geneIDs","logFC","pvalue","logpval")
plot.deg2 <- deg.2.df %>% 
	tibble::rownames_to_column("geneIDs") %>% 
	dplyr::select(geneIDs, avg_log2FC, p_val) %>% 
	mutate(logpval = -log10(p_val))
colnames(plot.deg2) <- c("geneIDs","logFC","pvalue","logpval")

# plot representative genes
volcanoPlot <- function(results.df,
								title = "Volcano Plot of Gene Expression",
								manual_labels = NULL) {
	
	# Create factor denoting differential expression status
	stats.df <- results.df %>%
		mutate(diffexpressed = ifelse(geneIDs %in% manual_labels, "Labelled", 
												"Not labelled"))
	
	stats.df$diffexpressed <- as.factor(stats.df$diffexpressed)
	
	# Base plot
	plot <- ggplot(stats.df, aes(x = logFC, y = logpval, col = diffexpressed, text = geneIDs)) +
		geom_point() +
		ggtitle(title) +
		ylab("-log10(p-val)") +
		ylim(0, 300) +  # Set the y-axis limit
		xlim(-1.5, 5)+
		xlab("log2FC")+
		geom_vline(xintercept = 0, col = "gray") +
		geom_hline(yintercept = 0, col = "gray") +
		scale_color_manual(values = c("#F96286", "gray")) +
		theme_classic()
	
	# Add geneIDs as labels for labelled genes
	labelled_genes <- stats.df %>%
		filter(geneIDs %in% manual_labels)
	
	plot <- plot +
		geom_label_repel(data = labelled_genes, 
							  aes(label = geneIDs), 
							  box.padding = 0.4, 
							  point.padding = 0.2, 
							  max.overlaps = 50)
	
	return(plot)
}

# Figure 4 A
volcanoPlot(plot.deg1,
				manual_labels = c("SPP1",
								  		"GPNMB",
								  		"S100A6",
								  		"PCLAF", # Proliferating cells markers
								  		"HMGN2",
								  		"HMGB2",
								  		"STMN1",
								  		cc.genes$s.genes, # Unsupervised proliferation markers
								  		cc.genes$g2m.genes)) + NoLegend()
volcanoPlot(plot.deg2,
				manual_labels = c("SPP1",
								  		"GPNMB",
								  		"S100A6",
								  		"PCLAF",
								  		"HMGN2",
								  		"HMGB2",
								  		"STMN1",
								  		cc.genes$s.genes,
								  		cc.genes$g2m.genes)) + NoLegend()

# Look at the damage markers across all nephron cell types
nephrons <- subset(kidneys,
						 subset = cell_types %in% c("PT",
						 									"PT_Prolif",
						 									"PT_VIM+",
						 									"DTC",
						 									"TAL_1",
						 									"TAL_2",
						 									"IC_TypeA",
						 									"IC_TypeB"))
# Figure 4 B
DefaultAssay(nephrons) <- "SCT"
VlnPlot(nephrons,
		  features = c(
		  "SPP1",
		  "GPNMB",
		  "S100A6"),
		  group.by = "Condition",
		  split.by = "cell_types",
		  pt.size = 0.01,
		  stack = TRUE,
		  combine = TRUE,
		  flip = TRUE,
		  same.y.lims = TRUE,
		  cols = c("#68D693", "#299CA8", "#01E1B1",
					  "#F29494", "#FA55AA","#FA8BF1",
					  "#FFE96E", "#FAB32F", "#FA7623",
					  "#6C82E6",
					  "#42D6FF"))+
			theme_classic() +
			theme(legend.position = "bottom",
					axis.title = element_blank())+
			ylim(0,8) +
			ggpubr::stat_compare_means(aes(group = Condition), 
												label = "p.signif", method = "wilcox.test", 
												comparisons = list( c("kidney_1st_untreated_rt","kidney_1st_treated"), 
																		  c("kidney_2nd_untreated_rt","kidney_2nd_treated")))


###### Enrichment Analysis #####
# Gene Ontology
DEG1_GeneID <- rownames(head(deg.1.df, n = 1000))
DEG2_GeneID <- rownames(head(deg.2.df, n = 1000))

GO_DEG1 <- enrichGO(gene = DEG1_GeneID,
						  OrgDb = org.Ss.eg.db, 
						  keyType = 'SYMBOL',
						  ont = "ALL", 
						  pAdjustMethod = "BH",
						  pvalueCutoff = 0.05, 
						  qvalueCutoff = 0.99) 
GO_DEG2 <- enrichGO(gene = DEG2_GeneID,
						  OrgDb = org.Ss.eg.db, 
						  keyType = 'SYMBOL',
						  ont = "ALL", 
						  pAdjustMethod = "BH",
						  pvalueCutoff = 0.05, 
						  qvalueCutoff = 0.99) 

# Supplementary Figure 2 A
dotplot(GO_DEG1) + theme_classic() + theme(axis.text = element_text(size = 14),
														 axis.text.x = element_text(angle = 90))
# Supplementary Figure 2 B
dotplot(GO_DEG2) + theme_classic() + theme(axis.text = element_text(size = 14),
														 axis.text.x = element_text(angle = 90))

# KEGG Pathways
DEG1_ENTREZID =bitr(DEG1_GeneID,
						  fromType="SYMBOL",
						  toType="ENTREZID",
						  OrgDb="org.Ss.eg.db")
DEG1_KEGG <- enrichKEGG(gene = DEG1_ENTREZID$ENTREZID,
								organism = "ssc",
								keyType = "kegg",
								pvalueCutoff = 0.05,
								qvalueCutoff = 0.99)
DEG2_ENTREZID =bitr(DEG2_GeneID,
						  fromType="SYMBOL",
						  toType="ENTREZID",
						  OrgDb="org.Ss.eg.db")
DEG2_KEGG <- enrichKEGG(gene = DEG2_ENTREZID$ENTREZID,
								organism = "ssc",
								keyType = "kegg",
								pvalueCutoff = 0.05,
								qvalueCutoff = 0.99) 

# Supplementary Figure 2 C
barplot(DEG1_KEGG) + theme_classic() + theme(axis.text = element_text(size = 14))
# Supplementary Figure 2 D
barplot(DEG2_KEGG) + theme_classic() + theme(axis.text = element_text(size = 14))

# Visualize Damage and Repair Genes across time
plotBulk <- function(genes){
	# Filter for selected genes
	subset_data <- bulk[bulk$gene %in% genes,]
	subset_data$gene <- factor(subset_data$gene, levels = genes)
	# Plot logFC over time
	ggplot(subset_data, aes(x = as.numeric(variable), y = value, group = gene, color = gene)) +
		geom_line(size = 1.5) +
		labs(x = "Hours after xenotransplantation", 
			  y = "log2FC Gene Expression") +
		theme_classic() +
		theme(legend.position="bottom",
				legend.text = element_text(size=16),
				legend.title = element_blank(),
				plot.title = element_blank())
}

# Figure 4 J
plotBulk(c("SPP1","GPNMB","S100A6",
			  "STMN1",
			  "PCLAF",
			  "HMGN2",
			  "HMGB2")) +
	scale_color_manual(values = c(rev(colorRampPalette(brewer.pal(3, "Oranges"))(3)), 
											rev(colorRampPalette(brewer.pal(4, "Greens"))(4))))
