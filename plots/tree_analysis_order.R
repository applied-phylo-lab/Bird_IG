order<-"Waterfowl"
gene_tree<-read.tree(paste0("/local/storage/kav67/Bird_data/",order,"/genes_aligned.treefile"))
annotation<-fread(paste0("/local/storage/kav67/Bird_data/",order,"/genes_annotation.tsv"))
p <- ggtree(gene_tree, layout="fan", open.angle=15, size=0.1)
p <- p %<+% annotation



p1 <-p +
  geom_tippoint(
    mapping=aes(colour=Species),
    size=1.5,
    stroke=0,
    alpha=0.8
  ) +
  theme(
    legend.title=element_text(size=11),
    legend.text=element_text(size=10),
    legend.spacing.y = unit(0.02, "cm")
  )

p2 <-p1 +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=Productive),
    width=0.2,
    offset=0.05
  ) +
  theme(
    legend.title=element_text(size=11), 
    legend.text=element_text(size=10),
    legend.spacing.y = unit(0.1, "cm")
  )

p3 <-p2 +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=Strand),
    width=0.2,
    offset=0.05
  ) +
  theme(
    legend.title=element_text(size=11), 
    legend.text=element_text(size=10),
    legend.spacing.y = unit(0.1, "cm")
  )
p3+labs(title=order)
