# figure 1: QGIS  ----------------------

# figure 2:  ----------------------

### histogram most abundant genera, families and phyla

Red_fam <- read.csv("outputs/Most_abundant_fam_input.csv", row.names=1) 

Red_fam <- Red_fam %>% 
  group_by(Family) %>% 
  summarise(across(everything(), ~ sum(.x)))

write.xlsx(Red_fam, file = "outputs/Most_abundant_fam.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

Red_phyla <- read.csv("outputs/Most_abundant_phylum_input.csv", row.names=1) 

Red_phyla <- Red_phyla %>% 
  group_by(Phylum) %>% 
  summarise(across(everything(), ~ sum(.x)))

write.xlsx(Red_phyla, file = "outputs/Most_abundant_phyla.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

### genera

abund_gen <- read.csv("outputs/Most_abundant_genera.csv")
colnames(abund_gen) <- c("genus","Percentages", "family", "phylum")

gen <- ggplot(abund_gen, aes(x=reorder(genus, -Percentages), y=Percentages, fill = phylum)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x="Genera accounting for >1% reads", y="Percentage of reads") +
  theme_bw() + 
  theme(legend.position = "none")

### family

abund_fam <- read.csv("outputs/Most_abundant_fam.csv")
colnames(abund_fam) <- c("family","Percentages", "phylum")

fam <- ggplot(abund_fam, aes(x=reorder(family, -Percentages), y=Percentages, fill = phylum)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x="Families accounting for >1% reads", y="Percentage of reads") +
  theme_bw() +
  theme(legend.position = "none")

### phyla

abund_phyla <- read.csv("outputs/Most_abundant_phyla.csv")
colnames(abund_phyla) <- c("phylum","Percentages")

phy <- ggplot(abund_phyla, aes(x=reorder(phylum, -Percentages), y=Percentages, fill = phylum)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x="Phyla accounting for >1% reads", y="Percentage of reads") +
  theme_bw() + 
  theme(legend.position = "none")

composite_barplots <- ggarrange(phy, fam, gen, labels = c("A", "B", "C"), ncol = 3, nrow = 1)
composite_barplots

ggsave("Figure_2_tiff.tiff", device='tiff', dpi=300)



# Figure 3: ----------------------
divdat <- read.csv("outputs/Diversity_indices_groups.csv")

my_comparisons4 = list( c("Agro_High", "Conv_High"), c("Agro_Low", "Conv_Low") )
p4 <- ggboxplot(
  divdat,
  x = "Conditions",
  y = "Faith_Phylogenetic_Diversity",
  fill = "Conditions",
  ylab = "Phylogenetic Diversity",
  xlab = "High altitude                   Low altitude                       ",
  legend = "none", #"none", "right",
  palette = c("darkolivegreen4", "darkgoldenrod3", "darkolivegreen2", "darkgoldenrod1")) + # palette = c("#00AFBB", "#E7B800", "#00AFBB", "#E7B800"))
  scale_x_discrete(labels = c("Agroecological", "Conventional", "Agroecological", "Conventional")) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), text = element_text(size = 10)) + 
  stat_compare_means(comparisons = my_comparisons4, label.y = c(185, 185),tip.length=0.0,  #label.y = c(0.01, 0.04)
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "***", "**", "*", "ns")))


#my_comparisons1 = list( c("Agro_High", "Conv_High"), c("Conv_High", "Agro_Low"), c("Conv_High", "Conv_Low") )
my_comparisons1 = list( c("Agro_High", "Conv_High"), c("Conv_Low", "Agro_Low"))
p1 <- ggboxplot(
  divdat,
  x = "Conditions",
  y = "S.ACE",
  # y = "new_x_exact",
  fill = "Conditions",
  ylab = "Species richness ACE ",
  #xlab = "High altitude                                       Low altitude",
  xlab = "",
  legend = "none", #"right",
  palette = c("darkolivegreen4", "darkgoldenrod3", "darkolivegreen2", "darkgoldenrod1")) + # palette = c("#00AFBB", "#E7B800", "#00AFBB", "#E7B800"))
  theme(text = element_text(size = 10)) + 
  scale_x_discrete(labels = c("", "", "", "")) +
  stat_compare_means(comparisons = my_comparisons1, label.y = c(580, 580),tip.length=0.0,  #label.y = c(0.01, 0.04)
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("*", "*", "*", "*", "ns")))

my_comparisons2 = list( c("Agro_High", "Conv_High"), c("Conv_Low", "Agro_Low"))
p2 <- ggboxplot(
  divdat,
  x = "Conditions",
  y = "Shannon",
  fill = "Conditions",
  ylab = "Shannon",
  xlab = "High altitude                   Low altitude                       ",
  legend = "none", #"none", "right",
  palette = c("darkolivegreen4", "darkgoldenrod3", "darkolivegreen2", "darkgoldenrod1")) + # palette = c("#00AFBB", "#E7B800", "#00AFBB", "#E7B800"))
  scale_x_discrete(labels = c("Agroecological", "Conventional", "Agroecological", "Conventional")) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), text = element_text(size = 10)) +
  stat_compare_means(comparisons = my_comparisons1, label.y = c(5.2, 5.2),tip.length=0.0,  #label.y = c(0.01, 0.04)
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "***", "**", "*", "ns")))


my_comparisons3 = list( c("Agro_High", "Conv_High"), c("Agro_Low", "Conv_Low"))
p3 <- ggboxplot(
  divdat,
  x = "Conditions",
  # y = "FourRTtransf_InvS",
  y = "Inverse_Simpson",
  fill = "Conditions",
  ylab = "Inverse Simpson",
  xlab =  "", #"High altitude                                       Low altitude",
  legend = "none", #"left"
  palette = c("darkolivegreen4", "darkgoldenrod3", "darkolivegreen2", "darkgoldenrod1")) + # palette = c("#00AFBB", "#E7B800", "#00AFBB", "#E7B800"))
  scale_x_discrete(labels = c("", "", "", "")) +
  theme(text = element_text(size = 10)) +
  stat_compare_means(comparisons = my_comparisons3, label.y = c(50, 50),tip.length=0.0,  #label.y = c(0.01, 0.04)
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "***", "**", "ns", "ns")))


composite <- ggarrange(p3, p1, p2, p4,labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
composite

ggsave("Figure_3_tiff.tiff", device='tiff', dpi=300)

# Figure 4:  ----------------------

## High altitude clr

PCA_df <- read.csv("outputs/PCA_input_coordinates_HIGH.csv", row.names = 1)
colnames(PCA_df) <- c("pcoa1", "pcoa2", "Management")
PCA_df %>% as_tibble(rownames="samples")

pca_High <- PCA_df %>% 
  ggplot(aes(x=pcoa1, y=pcoa2, color = Management)) +
  geom_point(aes(col = Management)) +
  labs(x="PCo 1 (16.7%)", y="PCo 2 (12.6%)",  title="Reads central log transformed") +
  stat_ellipse(level = 0.95) +
  theme_bw() + 
  theme(legend.position = "none") +
  scale_color_manual(values = c("darkolivegreen4", 
                                "darkgoldenrod3"))  


## Low altitude clr

PCA_df_low <- read.csv("outputs/PCA_input_coordinates_LOW.csv", row.names = 1)
colnames(PCA_df_low) <- c("pcoa1", "pcoa2", "Management")
PCA_df_low %>% as_tibble(rownames="samples")

pca_Low <- PCA_df_low %>% 
  ggplot(aes(x=pcoa1, y=pcoa2, color = Management)) +
  geom_point(aes(col = Management)) +
  labs(x="PCo 1 (13.3%)", y="PCo 2 (11.1%)", title="") +
  stat_ellipse(level = 0.95) +
  theme_bw() + 
  theme(legend.position = "none") + 
  scale_color_manual(values = c("darkolivegreen2", 
                                "darkgoldenrod1"))

## High altitude frequencies

PCA_df_freq <- read.csv("outputs/PCA_Freq_HIGH.csv", row.names = 1)
colnames(PCA_df_freq) <- c("pcoa1", "pcoa2", "Management")
PCA_df_freq %>% as_tibble(rownames="samples")

pca_High_freq <- PCA_df_freq %>% 
  ggplot(aes(x=pcoa1, y=pcoa2, color = Management)) +
  geom_point(aes(col = Management)) +
  labs(x="PCo 1 (19.7%)", y="PCo 2 (12.7%)", title="Reads in frequencies") +
  stat_ellipse(level = 0.95) +
  theme_bw() + 
  scale_color_manual(values = c("darkolivegreen4", 
                                "darkgoldenrod3"))

## Low altitude frequencies
PCA_df_low_freq <- read.csv("outputs/PCA_Freq_LOW.csv", row.names = 1)
colnames(PCA_df_low_freq) <- c("pcoa1", "pcoa2", "Management")
PCA_df_low_freq %>% as_tibble(rownames="samples")

pca_Low_freq <- PCA_df_low_freq %>% 
  ggplot(aes(x=pcoa1, y=pcoa2, color = Management)) +
  geom_point(aes(col = Management)) +
  labs(x="PCo 1 (16.9%)", y="PCo 2 (10.8%)", title="") +
  stat_ellipse(level = 0.95) +
  theme_bw() + 
  scale_color_manual(values = c("darkolivegreen2", 
                                "darkgoldenrod1"))


pca_Low
pca_High
pca_High_freq
pca_Low_freq

composite_pca <- ggarrange(pca_High, pca_High_freq, pca_Low, pca_Low_freq, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2, widths = c(0.45, 0.6, 0.45, 0.6))
composite_pca

ggsave("Figure_4_tiff.tiff", device='tiff', dpi=300)

# Figure 5:  ----------------------

# boxplots abundances per condition 
# Aldex_boxplot <- read.csv("outputs/Aldex_abundant_genera_boxplot_input.csv", row.names = 1)
Aldex_boxplot <- read.csv("outputs/Aldex_abundant_genera_boxplot_input_freq.csv", row.names = 1)

box_Romboutsia <-
  boxplot(Aldex_boxplot$Romboutsia_agro,
          Aldex_boxplot$Romboutsia_conv,
          main = "Differential abundance of the genus Romboutsia",
          #xlab="Agroecological                                 Conventional", 
          ylab="Abundance in percentages (%)", 
          col = c("darkolivegreen4", "darkgoldenrod3"),
          names=c("Agroecological", "Conventional"))

legend("topleft", inset=.02, title="Management", c("Agroecological", "Conventional"),
       fill = c("darkolivegreen4", "darkgoldenrod3"), horiz=FALSE, cex=1.2)


library(Cairo)
Cairo(file="Figure_5.png", width=2000, height=2000, type="png", dpi=300, bg = "white")
Cairo(file="Figure_5.tiff", width=2000, height=2000, type="tiff", dpi=300, bg = "white", fallback_resolution=300)
boxplot(Aldex_boxplot$Romboutsia_agro,
          Aldex_boxplot$Romboutsia_conv,
          main = "Differential abundance of the genus Romboutsia",
          ylab="Abundance in percentages (%)", 
          col = c("darkolivegreen4", "darkgoldenrod3"),
          names=c("Agroecological", "Conventional"))

  legend("topleft", inset=.02, title="Management", c("Agroecological", "Conventional"),
       fill = c("darkolivegreen4", "darkgoldenrod3"), horiz=FALSE, cex=1.2)

dev.off()

# ggsave("Figure_5_tiff.tiff", device='tiff', dpi=300)


# Figure 6:  ----------------------
box_others <- 
  boxplot(Aldex_boxplot$Lysinibacillus_agro,
          Aldex_boxplot$Lysinibacillus_conv,
          Aldex_boxplot$Empedobacter_agro,
          Aldex_boxplot$Empedobacter_conv,
          Aldex_boxplot$Propionispira_agro,
          Aldex_boxplot$Propionispira_conv,
          Aldex_boxplot$Erysipelothrix_agro,
          Aldex_boxplot$Erysipelothrix_conv,
          main = "Genera differing significant, but precent in low abundances (<0.5% reads)",
          ylab="Abundance in percentages (%)", 
          names=c("                         Lysinibacillus","", "                         Empedobacter","","                        Propionispira","","                       Erysipelothrix", ""),
          col = c("darkolivegreen4", "darkgoldenrod3","darkolivegreen4", "darkgoldenrod3","darkolivegreen4", "darkgoldenrod3","darkolivegreen4", "darkgoldenrod3"))

legend("topleft", inset=.02, title="Management", c("Agroecological", "Conventional"),
       fill = c("darkolivegreen4", "darkgoldenrod3"), horiz=FALSE, cex=1.2)

Cairo(file="Figure_6.tiff", width=3500, height=2000, type="tiff", dpi=300, bg = "white", fallback_resolution=300)
boxplot(Aldex_boxplot$Lysinibacillus_agro,
                Aldex_boxplot$Lysinibacillus_conv,
                Aldex_boxplot$Empedobacter_agro,
                Aldex_boxplot$Empedobacter_conv,
                Aldex_boxplot$Propionispira_agro,
                Aldex_boxplot$Propionispira_conv,
                Aldex_boxplot$Erysipelothrix_agro,
                Aldex_boxplot$Erysipelothrix_conv,
                #main = "Genera differing significant, but precent in low abundances (<0.5% reads)",
                main= "Low abundance genera",
                ylab="Abundance in percentages (%)", 
                #names=c("                         Lysinibacillus","","                         Empedobacter","","                        Propionispira","","                       Erysipelothrix"),
                xlab=expression(italic("Lysinibacillus                               Empedobacter                              Propionispira                             Erysipelothrix")),
                col = c("darkolivegreen4", "darkgoldenrod3","darkolivegreen4", "darkgoldenrod3","darkolivegreen4", "darkgoldenrod3","darkolivegreen4", "darkgoldenrod3"))
legend("topleft", inset=.02, c("Agroecological", "Conventional"),
       fill = c("darkolivegreen4", "darkgoldenrod3"), horiz=FALSE, cex=1.2)

dev.off()
