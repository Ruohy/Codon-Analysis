# Analysing Metagenomics Data
# Installing packages
install.packages("ggrepel")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyr")
install.packages("RColorBrewer")
install.packages("psych")
install.packages("MuMIn")
install.packages("lme4")

#### Link GitHub ####
install.packages("usethis")
library(usethis)
use_git_config(user.name="Ruohy", user.email="rhys.tuohy@hdr.mq.edu.au")
use_git()

# create token
create_github_token()
token <- Sys.getenv("GITHUB_PAT")
install.packages("gitcreds")
library(gitcreds)
gitcreds_set()

# connect specific project
use_github()

#### Data exploration (hclust + errors/fixes) ####

# load data
Codon_usage<-read.csv("Data/codon_usage.csv") #load data
View(Codon_usage) #inspect data manually

# heirarchial clustering showing codon frequency relatedness
cdn<-Codon_usage[,6:ncol(Codon_usage)] #select only codon frequency columns
plot(hclust(dist(t(cdn)))) # cluster to look at related/unrelated data
# we see an outlier datapoint (CUA)

#perform a principal component analysis (PCA)
p<-princomp(cdn) #error with data (contains infinite values?)

#so we investigate problematic values by summing
sum(cdn) #says there are non-numeric values in the matrix

#data was causing problems, need to define non-numeric values
table(unlist(cdn)) # this makes it easy to manually spot non-numbers like "-"
cdn[cdn=="-"]=0 # change dash to 0 to make it numeric
table(unlist(cdn)) #it worked
sum(cdn) # but sum still wont work 

sum(is.na(cdn)) # check for existance of NA values
bad<-which(is.na(cdn)) #find position of NA values
unlist(cdn)[bad] # prints the list so we can see the inputs of problem values

#change these values to 0 too!
cdn[cdn=="non-B hepatitis virus"]=0 
cdn[cdn=="12;I"]=0

#convert entire matrix to numeric
cdn=matrix(as.numeric(unlist(cdn)),nrow(cdn),ncol(cdn))




#

#### PCA + initial visual exploration of $loadings [RUN START] ####

#load data
Codon_usage<-read.csv("Data/codon_usage.csv") #load data
View(Codon_usage) #inspect data manually

cdn=matrix(as.numeric(unlist(cdn)),nrow(cdn),ncol(cdn))
p<-princomp(cdn) #run PCA worked!
View(p) #shows princomp data

p$loadings[1:30,]
plot(p$loadings) #made a scatter plot of loadings
colnames(cdn)=colnames(Codon_usage[,6:ncol(Codon_usage)]) 
p$loadings[7,] #shows all the rows
shade=c(rep("coral",16),rep("green2",16),rep("deepskyblue",16),rep("gold",16))
#create a list of each colour in repeates of 16
plot(p$loadings,col=shade,pch=19) 
#plotting different colours based on codon loadings to examine patterns??
#seems a littel scattered so far
regexpr("^U",colnames(cdn)) #an attempt to find a pattern in codon 1st letter U
#colour the first letter of each differently
shade[regexpr("^U",colnames(cdn))==1]="coral"
shade[regexpr("^C",colnames(cdn))==1]="green2"
shade[regexpr("^A",colnames(cdn))==1]="deepskyblue"
shade[regexpr("^G",colnames(cdn))==1]="gold"
plot(p$loadings,col=shade,pch=19) #plot again to see how it looks
#a little more clustering, maybe something?

# p$loadings = which codons contribute most to patterns

#### Plot codon frequencies against kingdoms ($scores = patterns!) ####

#try plotting relating to kingdoms
plot(p$scores,cex=0.1) # a definite pattern!
#lets bring out some colors relating to where the kingdoms lay on this plot
points(p$scores[Codon_usage$Kingdom=="vrl",],col="coral",cex=0.2)
points(p$scores[Codon_usage$Kingdom %in% 
                  c("mam","rod","pri"),],col="deepskyblue",cex=0.2)

#grouping some kingdoms together to reduce unnecessary clutter
Codon_usage$Group <- NA
Codon_usage$Group[Codon_usage$Kingdom %in% c("mam", "pri", "rod")] <- "Mammals"
Codon_usage$Group[Codon_usage$Kingdom == "vrt"] <- "Vertebrates"
Codon_usage$Group[Codon_usage$Kingdom == "inv"] <- "Invertebrates"
Codon_usage$Group[Codon_usage$Kingdom == "bct"] <- "Bacteria"
Codon_usage$Group[Codon_usage$Kingdom == "arc"] <- "Archaea"
Codon_usage$Group[Codon_usage$Kingdom %in% c("vrl", "phg")] <- "Viruses"
Codon_usage$Group[Codon_usage$Kingdom == "pln"] <- "Plants"
Codon_usage$Group[Codon_usage$Kingdom %in% "plm"] <- "Protists"

pca_df <- data.frame(p$scores)
pca_df$Group <- Codon_usage$Group
pca_df$Kingdom <- Codon_usage$Kingdom

#assign colors to these groups
group_colors <- c(
  "Mammals" = "deepskyblue",
  "Vertebrates" = "dodgerblue",
  "Invertebrates" = "turquoise4",
  "Bacteria" = "firebrick",
  "Archaea" = "tomato",
  "Viruses" = "darkorange2",
  "Plants" = "forestgreen",
  "Protists" = "darkviolet"
)
# Create plot using scores from PCA
plot(p$scores, type = "n", main = "Codon Usage PCA by Biological Group")
points(p$scores, 
       col = group_colors[Codon_usage$Group], 
       pch = 19, cex = 0.4)
# Add a legend
legend("bottomright", legend = names(group_colors), 
       col = group_colors, pch = 19, cex = 0.8)

#Maybe make a nicer looking plot using ggplot2
  

#### Codon PCA plot with ggplot2 ####

library(ggplot2) # load package

# create data frame with PCA data and metadata
pca_df <- data.frame(p$scores)
pca_df$Kingdom <- Codon_usage$Kingdom
pca_df$Group <- Codon_usage$Group
df_pca <- data.frame(p$scores, Kingdom = Codon_usage$Kingdom)

kingdom_group_plot<-ggplot(pca_df, aes(x = Comp.1, y = Comp.2, color = Group)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = group_colors) +
  theme_minimal() +
  labs(
    title = "PCA of Codon Usage by Biological Group",
    x = "PC1",
    y = "PC2",
    color = "Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )
print(kingdom_group_plot)

# add labelled ellipses to show clusters

library(dplyr)
library(ggrepel)

# Group centroids only (no manual nudging yet)
centroids <- pca_df %>%
  group_by(Group) %>%
  summarise(x = mean(Comp.1), y = mean(Comp.2))

# Plot with auto-positioned repel labels
kingdom_group_ellipses<-ggplot(pca_df, 
                               aes(x = Comp.1, y = Comp.2, color = Group)) +
  geom_point(size = 1.5, alpha = 0.75) +
  stat_ellipse(type = "norm", linetype = 2) +
  geom_text_repel(
    data = centroids,
    aes(x = x, y = y, label = Group),
    size = 4,
    fontface = "bold",
    segment.size = 0.3,
    segment.color = "grey40",
    color = "black",
    force = 1,           # slightly repels away from centroid
    max.overlaps = Inf,  # ensure all are shown
    show.legend = FALSE
  ) +
  scale_color_manual(values = group_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA with Clean Ellipse Labels",
    x = "PC1", y = "PC2", color = "Group"
  )
print(kingdom_group_ellipses)
# zoom in on a specific group extra function

# zooming in on viruses
subset_df <- subset(pca_df, Group == "Viruses")
ggplot(subset_df, aes(x = Comp.1, y = Comp.2)) +
  geom_point(color = "darkorange2") +
  theme_minimal() +
  labs(title = "Codon Usage Within Viruses")

# zooming in on bacteria
subset_df <- subset(pca_df, Group == "Bacteria")
ggplot(subset_df, aes(x = Comp.1, y = Comp.2)) +
  geom_point(color = "firebrick") +
  theme_minimal() +
  labs(title = "Codon Usage Within Bacteria")

# colours and clustering make it hard to understand what is happening,
# best to try to seperate PCA distributions to better compare them

# Change the colours?
custom_palette <- c(
  "Mammals" = "#1f78b4",         
  "Vertebrates" = "#559ccf",     
  "Invertebrates" = "#76b76a",   
  "Bacteria" = "#e31a1c",        
  "Archaea" = "#d85f71",         
  "Viruses" = "#ff7f00",         
  "Plants" = "#33a02c",          
  "Protists" = "#6a3d9a"         
)

# Re-plot with fixed colours
library(ggplot2)
KG_PCA_fix<-ggplot(pca_df, aes(x = Comp.1, y = Comp.2, color = Group)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = custom_palette) +
  theme_minimal() +
  labs(
    title = "Codon Usage Patterns by Biological Group",
    x = "PC1",
    y = "PC2",
    color = "Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  )

print(KG_PCA_fix)
ggsave(
  filename = "Output/KG_PCA_fix.png",
  plot = KG_PCA_fix,   # optional if it's the last plot created
  width = 8,                 # in inches
  height = 6,
  dpi = 300,
  bg = "white"
)

# Split the view between groups
library(ggplot2)
library(dplyr)

pca_df$Group <- factor(pca_df$Group)

gray_df <- expand.grid(FacetGroup = levels(pca_df$Group),
                       Row = 1:nrow(pca_df)) %>%
  mutate(Comp.1 = pca_df$Comp.1[Row],
         Comp.2 = pca_df$Comp.2[Row])

color_df <- pca_df %>%
  mutate(FacetGroup = Group)

KG_PCA_split <- ggplot() +
  geom_point(data = gray_df,
             aes(x = Comp.1, y = Comp.2),
             color = "gray80", size = 1.2, alpha = 0.5) +
  geom_point(data = color_df,
             aes(x = Comp.1, y = Comp.2, color = Group),
             size = 1.6, alpha = 0.85) +
  facet_wrap(~ FacetGroup, scales = "fixed") +
  scale_color_manual(values = custom_palette) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none") +
  labs(
    title = "Codon Usage PCA — Group Highlighted per Panel",
    x = "PC1", y = "PC2"
  )

print(KG_PCA_split)

ggsave(
  filename = "Output/KG_PCA_split.png",
  plot = KG_PCA_split,   # optional if it's the last plot created
  width = 8,                 # in inches
  height = 6,
  dpi = 300,
  bg = "white"
)


#### Codon usage PCA by domain ####
Codon_usage$Domain <- NA
Codon_usage$Domain[Codon_usage$Kingdom %in%
                    c("pri", "rod", "mam", "pln",
                      "vrt", "inv", "fun", "plm", "alg" )] <- "Eukarya"
Codon_usage$Domain[Codon_usage$Kingdom %in% c("bct")] <- "Bacteria"
Codon_usage$Domain[Codon_usage$Kingdom %in% c("arc")] <- "Archaea"
Codon_usage$Domain[Codon_usage$Kingdom %in% c("vrl", "phg")] <- "Virus"

# Contrasting colours
domain_colors <- c(
  "Bacteria" = "#E41A1C",
  "Archaea"  = "#984EA3",
  "Eukarya"  = "#4DAF4A",
  "Virus"    = "#FF7F00"
)

#plot with ggplot

df_pca <- data.frame(p$scores)

library(ggplot2)
library(dplyr)

df_pca$Domain <- Codon_usage$Domain

domain_plot<-ggplot(df_pca, aes(x = Comp.1, y = Comp.2, color = Domain)) +
  geom_point(size = 1, alpha = 0.85) +
  scale_color_manual(values = domain_colors) +
  theme_minimal(base_size = 15)

print(domain_plot)
ggsave(
  filename = "Output/Domain_plot.png",
  plot = domain_plot, 
  width = 8,                 
  height = 6,
  dpi = 300,
  bg = "white"
)

# plot 4 domains separately and panel because there's too much overlap

# Ensure Domain is a factor
df_pca$Domain <- factor(df_pca$Domain)

gray_df <- expand.grid(FacetDomain = levels(df_pca$Domain),
                       Row = 1:nrow(df_pca)) %>%
  mutate(Comp.1 = df_pca$Comp.1[Row],
         Comp.2 = df_pca$Comp.2[Row])

color_df <- df_pca %>%
  mutate(FacetDomain = Domain)


domain_split_panel <- ggplot() +
  geom_point(data = gray_df,
             aes(x = Comp.1, y = Comp.2),
             color = "gray80", size = 1.2, alpha = 0.5) +
  geom_point(data = color_df,
             aes(x = Comp.1, y = Comp.2, color = Domain),
             size = 1.6, alpha = 0.85) +
  facet_wrap(~ FacetDomain) +
  scale_color_manual(values = domain_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Codon Usage PCA — Domain Highlighted per Panel",
    x = "PC1", y = "PC2"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

print(domain_split_panel)
ggsave(
  filename = "Output/Domain_split_panel.png",
  plot = domain_split_panel,   
  width = 8,                 
  height = 6,
  dpi = 300,
  bg = "white"
)

#
#### Amino acid <- Codon plot by Kingdom ####

# Define codon-to-amino-acid map
cdn_AA <- c("UUU"="Phe","UUC"="Phe","UUA"="Leu","UUG"="Leu(S)","CUU"="Leu",
            "CUC"="Leu","CUA"="Leu","CUG"="Leu","AUU"="Ile","AUC"="Ile",
            "AUA"="Ile","AUG"="Met(S)","GUU"="Val","GUC"="Val","GUA"="Val",
            "GUG"="Val(S)","UCU"="Ser","UCC"="Ser","UCA"="Ser","UCG"="Ser",
            "CCU"="Pro","CCC"="Pro","CCA"="Pro","CCG"="Pro","ACU"="Thr",
            "ACC"="Thr","ACA"="Thr","ACG"="Thr","GCU"="Ala","GCC"="Ala",
            "GCA"="Ala","GCG"="Ala","UAU"="Tyr","UAC"="Tyr","UAA"="STOP",
            "UAG"="STOP","CAU"="His","CAC"="His","CAA"="Gln","CAG"="Gln",
            "AAU"="Asn","AAC"="Asn","AAA"="Lys","AAG"="Lys","GAU"="Asp",
            "GAC"="Asp","GAA"="Glu","GAG"="Glu","UGU"="Cys","UGC"="Cys",
            "UGA"="STOP","UGG"="Trp","CGU"="Arg","CGC"="Arg","CGA"="Arg",
            "CGG"="Arg","AGU"="Ser","AGC"="Ser","AGA"="Arg","AGG"="Arg",
            "GGU"="Gly","GGC"="Gly","GGA"="Gly","GGG"="Gly")

# Create long format codon data
library(tidyr)
cdn_aa_df <- as.data.frame(cdn)
cdn_aa_df$Sample <- 1:nrow(cdn_aa_df)

cdn_long <- pivot_longer(cdn_aa_df, -Sample, names_to = "Codon", values_to = "Freq")
cdn_long$AminoAcid <- cdn_AA[cdn_long$Codon]

# Average codon frequency per amino acid per sample
library(dplyr)
amino_summary <- cdn_long %>%
  group_by(Sample, AminoAcid) %>%
  summarise(MeanFreq = mean(as.numeric(Freq)), .groups = "drop")

# Convert to wide format: samples × amino acids
cdn_amino <- pivot_wider(amino_summary, names_from = AminoAcid, values_from = MeanFreq)
cdn_amino <- as.data.frame(cdn_amino[,-1])  # drop sample column

# PCA on amino acid usage
amino_pca <- prcomp(cdn_amino, scale. = TRUE)
amino_pca_df <- data.frame(amino_pca$x)
amino_pca_df$Group <- Codon_usage$Group
amino_pca_df$Kingdom <- Codon_usage$Kingdom

# Plot PCA of amino acid usage
library(ggplot2)
Amino_PCA<-ggplot(amino_pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size =2, alpha = 0.8) +
  scale_color_manual(values = custom_palette) +
  theme_minimal(base_size = 15) +
  labs(
    title = "PCA of Amino Acid Usage by Group",
    x = paste0("PC1 (", round(100 * summary(amino_pca)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(amino_pca)$importance[2, 2], 1), "%)"),
    color = "Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

print(Amino_PCA)
ggsave(
  filename = "Output/Amino_PCA_Kingdom.png",
  plot = Amino_PCA, 
  width = 8,                 
  height = 6,
  dpi = 300,
  bg = "white"
)


# Amino acid usage by Kingdom split plot
library(dplyr)
library(ggplot2)

amino_pca_df$Group <- factor(amino_pca_df$Group)

gray_amino_df <- expand.grid(FacetGroup = levels(amino_pca_df$Group),
                             Row = 1:nrow(amino_pca_df)) %>%
  mutate(PC1 = amino_pca_df$PC1[Row],
         PC2 = amino_pca_df$PC2[Row])

color_amino_df <- amino_pca_df %>%
  mutate(FacetGroup = Group)

Amino_PCA_split <- ggplot() +
  geom_point(data = gray_amino_df,
             aes(x = PC1, y = PC2),
             color = "gray80", size = 1.2, alpha = 0.5) +
  geom_point(data = color_amino_df,
             aes(x = PC1, y = PC2, color = Group),
             size = 1.8, alpha = 0.85) +
  facet_wrap(~ FacetGroup, scales = "fixed") +
  scale_color_manual(values = custom_palette) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none") +
  labs(
    title = "Amino Acid Usage PCA — Group Highlighted per Panel",
    x = "PC1", y = "PC2"
  )

print(Amino_PCA_split)
ggsave(
  filename = "Output/Amino_PCA_split.png",
  plot = Amino_PCA_split, 
  width = 8,                 
  height = 6,
  dpi = 300,
  bg = "white"
)

#### Amino acid <- Codon plot by Domain ####

amino_pca <- prcomp(cdn_amino, scale. = TRUE)
amino_pca_df <- data.frame(amino_pca$x)
amino_pca_df$Group <- Codon_usage$Group
amino_pca_df$Kingdom <- Codon_usage$Kingdom

amino_pca_df$Domain <- NA
amino_pca_df$Domain[Codon_usage$Kingdom %in% 
                      c("pri", "rod", "mam", "pln", "vrt", "inv", "fun", "plm", "alg")] <- "Eukarya"
amino_pca_df$Domain[Codon_usage$Kingdom == "bct"] <- "Bacteria"
amino_pca_df$Domain[Codon_usage$Kingdom == "arc"] <- "Archaea"
amino_pca_df$Domain[Codon_usage$Kingdom %in% c("vrl", "phg")] <- "Virus"

Amino_PCA_Domain <- ggplot(amino_pca_df, aes(x = PC1, y = PC2, color = Domain)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = domain_colors) +
  theme_minimal(base_size = 15) +
  labs(
    title = "Amino Acid Usage PCA by Domain",
    x = paste0("PC1 (", round(100 * summary(amino_pca)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(amino_pca)$importance[2, 2], 1), "%)"),
    color = "Domain"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(Amino_PCA_Domain)
ggsave(
  filename = "Output/Amino_PCA_Domain.png",
  plot = Amino_PCA_Domain, 
  width = 8,                 
  height = 6,
  dpi = 300,
  bg = "white"
)

# Panel plots separately
library(dplyr)
library(ggplot2)

amino_pca_df$Domain <- factor(amino_pca_df$Domain)

gray_domain_df <- expand.grid(FacetDomain = levels(amino_pca_df$Domain),
                              Row = 1:nrow(amino_pca_df)) %>%
  mutate(PC1 = amino_pca_df$PC1[Row],
         PC2 = amino_pca_df$PC2[Row])

color_domain_df <- amino_pca_df %>%
  mutate(FacetDomain = Domain)

Amino_PCA_Domain_Split <- ggplot() +
  geom_point(data = gray_domain_df,
             aes(x = PC1, y = PC2),
             color = "gray80", size = 1.2, alpha = 0.5) +
  geom_point(data = color_domain_df,
             aes(x = PC1, y = PC2, color = Domain),
             size = 1.8, alpha = 0.85) +
  facet_wrap(~ FacetDomain, scales = "fixed") +
  scale_color_manual(values = domain_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Amino Acid Usage PCA — Domain Highlighted per Panel",
    x = "PC1", y = "PC2"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

print(Amino_PCA_Domain_Split)
ggsave(
  filename = "Output/Amino_PCA_Domain_split.png",
  plot = Amino_PCA_Domain_Split, 
  width = 8,                 
  height = 6,
  dpi = 300,
  bg = "white"
)

#### Linear Mixed Effects Regression of Amino Acid w/ PCA ####
library(lme4)
library(MuMIn)  # for r.squaredGLMM()

amino_pca_df$Group <- as.factor(amino_pca_df$Group)
amino_pca_df$Kingdom <- as.factor(amino_pca_df$Kingdom)
# Make sure Group and Kingdom are factors
m <- lmer(PC1 ~ Group + (1 | Kingdom), data = amino_pca_df)
summary(m)
ranef(m)
rand(m)
r.squaredGLMM(m)

# Fixed effect = differences between groups
# Random effect = variance within the groups
# Kingdom adds meaningful structure but is secondary.

# Violin plot to show variance
# Wider spread = greater intra-group variability in amino acid usage
# Tighter spread = more conserved patterns of usage across species in that group
# Higher PC1 score = species use certain amino acids more in different proportions

library(ggplot2)
Amino_acid_pca<-ggplot(amino_pca_df, aes(x = Group, y = PC1, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  scale_fill_manual(values = custom_palette) +
  theme_minimal(base_size = 14) +
  labs(title = "Amino Acid Usage Patterns (PC1) by Biological Group",
       x = "Group", y = "PC1 Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") + 
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "white") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2)

print(Amino_acid_pca)
ggsave(
  filename = "Output/Amino_Acid_PCA_Plot.png",
  plot = Amino_acid_pca, 
  width = 8,                 
  height = 6,
  dpi = 300,
  bg = "white"
)

#### GC + GC3 content ####

# calculate proportion of codons that contain G or C or both in any position
gc_codons <- colnames(cdn)[grepl("[GC]",colnames(cdn))] #isolate only GC codons
GC_content <- rowSums(cdn[, gc_codons], na.rm = TRUE)/1.00008
# scaled down by 1.00008 because of a miscalculation that caused a range of over 1.0

# Filter only codons that have G or C in the third position
gc3_codons <- colnames(cdn)[substr(colnames(cdn), 3, 3) %in% c("G", "C")]
GC3_content <- rowSums(cdn[, gc3_codons], na.rm = TRUE)

# Add to dataframe
amino_pca_df$GC_content <- GC_content
amino_pca_df$GC3_content <- GC3_content

library(ggplot2)

#GC content
GC_PC1_plot<-ggplot(amino_pca_df, aes(x = GC_content, y = PC1, color = Group)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = custom_palette) +
  theme_minimal(base_size = 14) +
  labs(title = "PC1 vs GC Content by Group",
       x = "GC Content", y = "PC1 (Amino Acid Usage)") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(GC_PC1_plot)
ggsave(
  filename = "Output/GC_PC1_plot.png",
  plot = GC_PC1_plot, 
  width = 8,                 
  height = 6,
  dpi = 300,
  bg = "white"
)

#GC content split plot
library(dplyr)
library(ggplot2)

amino_pca_df$Group <- factor(amino_pca_df$Group)

gray_gc_df <- expand.grid(FacetGroup = levels(amino_pca_df$Group),
                          Row = 1:nrow(amino_pca_df)) %>%
  mutate(GC_content = amino_pca_df$GC_content[Row],
         PC1 = amino_pca_df$PC1[Row])

color_gc_df <- amino_pca_df %>%
  mutate(FacetGroup = Group)

GC_PC1_facet_plot <- ggplot() +
  geom_point(data = gray_gc_df,
             aes(x = GC_content, y = PC1),
             color = "gray80", size = 1, alpha = 0.4) +
  geom_point(data = color_gc_df,
             aes(x = GC_content, y = PC1, color = Group),
             size = 1.5, alpha = 0.85) +
  geom_smooth(data = color_gc_df,
              aes(x = GC_content, y = PC1),
              method = "lm", se = TRUE, linewidth = 0.8, color = 'black') +
  facet_wrap(~ FacetGroup, scales = "fixed") +
  scale_color_manual(values = custom_palette) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PC1 vs GC Content — Group Highlighted per Panel",
    x = "GC Content", y = "PC1 (Amino Acid Usage)"
  ) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(GC_PC1_facet_plot)
ggsave(
  filename = "Output/GC_PC1_facet_plot.png",
  plot = GC_PC1_facet_plot, 
  width = 8,                 
  height = 6,
  dpi = 300,
  bg = "white"
)

#GC3 content
GC3_PC1_plot<-ggplot(amino_pca_df, aes(x = GC3_content, y = PC1, color = Group)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = custom_palette) +
  theme_minimal(base_size = 14) +
  labs(title = "PC1 vs GC3 Content by Group",
       x = "GC3 Content", y = "PC1 (Amino Acid Usage)") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(GC3_PC1_plot)
ggsave(
  filename = "Output/GC3_PC1_plot.png",
  plot = GC3_PC1_plot, 
  width = 8,                 
  height = 6,
  dpi = 300,
  bg = "white"
)

#GC3 content split plot
library(dplyr)
library(ggplot2)

amino_pca_df$Group <- factor(amino_pca_df$Group)

gray_gc3_df <- expand.grid(FacetGroup = levels(amino_pca_df$Group),
                           Row = 1:nrow(amino_pca_df)) %>%
  mutate(GC3_content = amino_pca_df$GC3_content[Row],
         PC1 = amino_pca_df$PC1[Row])

color_gc3_df <- amino_pca_df %>%
  mutate(FacetGroup = Group)

GC3_PC1_facet_plot <- ggplot() +
  geom_point(data = gray_gc3_df,
             aes(x = GC3_content, y = PC1),
             color = "gray80", size = 1, alpha = 0.4) +
  geom_point(data = color_gc3_df,
             aes(x = GC3_content, y = PC1, color = Group),
             size = 1.5, alpha = 0.85) +
  geom_smooth(data = color_gc3_df,
              aes(x = GC3_content, y = PC1),
              method = "lm", se = TRUE, size = 0.8, color = "black") +
  facet_wrap(~ FacetGroup, scales = "fixed") +
  scale_color_manual(values = custom_palette) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PC1 vs GC3 Content — Group Highlighted per Panel",
    x = "GC3 Content", y = "PC1 (Amino Acid Usage)"
  ) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(GC3_PC1_facet_plot)
ggsave(
  filename = "Output/GC3_PC1_facet_plot.png",
  plot = GC3_PC1_facet_plot, 
  width = 8,                 
  height = 6,
  dpi = 300,
  bg = "white"
)

#### Meaningful LMER with PCA ####

q <- qnorm(GC_content)

m1 <- lmer(PC1 ~ Codon_usage$Ncodons + qnorm(GC_content) + qnorm(GC3_content) + 
            (1 | Group) + (1| Codon_usage$DNAtype) + 
             (1 | Kingdom), data = amino_pca_df)
summary(m1)
ranef(m1)
rand(m1)
r.squaredGLMM(m1)

m2 <- lmer(PC2 ~ Codon_usage$Ncodons + qnorm(GC_content) + qnorm(GC3_content) + 
            (1 | Group) + (1| Codon_usage$DNAtype) + 
             (1 | Kingdom), data = amino_pca_df)
summary(m2)
ranef(m2)
rand(m2)
r.squaredGLMM(m2)

#### Quick tests done with John ####
table(Codon_usage$DNAtype)

m <- lmer(amino_pca_df$PC1 ~ log(Ncodons) + (1 | Group) + (1 | DNAtype) + (1 | Kingdom), data = Codon_usage)
summary(m)
ranef(m)
rand(m)
r.squaredGLMM(m)


m <- lmer(amino_pca_df$PC2 ~ log(Ncodons) + (1 | Group) + (1 | DNAtype) + (1 | Kingdom), data = Codon_usage)
summary(m)
ranef(m)
rand(m)
r.squaredGLMM(m)

hist(Codon_usage$Ncodons)

m <- lmer(log(Ncodons) ~ amino_pca_df$PC2 + amino_pca_df$PC1 + (1 | Group) + (1 | DNAtype), data = Codon_usage)
summary(m)
ranef(m)
rand(m)
r.squaredGLMM(m)

hist(qnorm(GC_content))

plot(hclust(dist(t(cdn_amino))))

plot(p$scores,cex=0.1), 
points(p$scores[Codon_usage$DNAtype==0,],col="coral",cex=0.5),
points(p$scores[Codon_usage$DNAtype==1,],col="blue",cex=0.5),
points(p$scores[Codon_usage$DNAtype==2,],col="green",cex=0.5)

plot(p$scores,cex=log(Codon_usage$Ncodons)/10, col=hsv(h=pnorm(scale(log(Codon_usage$Ncodons)))))

     
#### Rearrange figures for report ####
