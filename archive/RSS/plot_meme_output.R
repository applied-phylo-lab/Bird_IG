library(readr)
library(ggplot2)
fimo <- read_tsv("/local/storage/kav67/Bird_data/fimo_out/fimo.tsv")
fimo$`p-value`<-as.numeric(fimo$`p-value`)
fimo$start<-as.numeric(fimo$start)
fimo$logp<--log10(fimo$`p-value`)

fimo<-fimo[fimo$start<100,]
  
ggplot(fimo, aes(x = start, y = logp, color = motif_id)) +
  geom_point() 

ggplot(fimo[fimo$start<5,], aes(x = start, y = logp, color = motif_id)) +
  geom_point() 

ggplot(fimo, aes(x = start)) +
  geom_histogram(bins = 50, alpha = 0.8) +
  facet_wrap(~ motif_id, ncol = 1, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Distribution of Motif Hits Position (Local 1–100 bp Window)",
    x = "Position in 100 bp downstream window",
    y = "Count"
  )
