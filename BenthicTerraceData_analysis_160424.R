##Terrace macrobenthic dataset
##PI Amanda Poste
#Script containing macrobenthic abundance and traits data for the 30 

#All dataset
##Scripts

#Packages
library (readxl)
library(gplots) #for heatwave plots
library(ggplot2) 
library(ggforce)
library(stringr)
library (reshape2)
library(sf)
library(vegan)
library(tidyverse)
library(dplyr)
library(ape) # for cluster and nMDS
library(dendextend) #dendogram
library(pvclust) #species fidelity
library(mFD)
library(labdsv) #indval
library(FactoMineR) #for PCA
library(factoextra) #for PCA, merge plots
library(missMDA)
library(ggcorrplot)
library(RColorBrewer)
library(ggrepel) #labels smaller
library(fmsb)
library(ade4) # for diversity analysis
library(stats) # for distance calculations
library(fundiversity) 

##Macrobenthic abundance, not transformed 
df <- read_excel("C:/Users/JPA/OneDrive - NIVA/PhD thesis work/Backup/Scripts/Terrace_Rfile_df_190324.xlsx") #Colummn name, first row data
df_hel <- decostand(df[c(22:210)], method = "hellinger")
dft <- data.frame(df)                 # Create a  copy of df
dft[c(22:210)] <- df_hel    # Replace columns 20 to 208 with the transformed data
dftrait <- read_excel("C:/Users/JPA/OneDrive - NIVA/PhD thesis work/Backup/Scripts/Terrace_Rfile_dftrait_110424.xlsx") #traitsdata
dftab <- read_excel("C:/Users/JPA/OneDrive - NIVA/PhD thesis work/Backup/Scripts/Terrace_Rfile_dftab_110424.xlsx") #hellinger transformed, for cwm calculation matching
dfperma 

## Environmental data
##General characterisation
environmental_data <- df[, c(8:17, 19:21)]
#by fjord branch
mean_per_fjordbranch <- aggregate(environmental_data, by = list(df$Fjord), FUN = mean)
mean_sd_per_fjordbranch <- aggregate(environmental_data, by = list(df$Fjord), FUN = function(x) c(mean = mean(x), sd = sd(x)))
mean_sd_per_fjordbranch[, -1] <- round(mean_sd_per_fjord[, -1], 2)
print(mean_sd_per_fjordbranch)
#by station type
mean_sd_per_stationtype <- aggregate(environmental_data, by = list(df$Stationtype), FUN = function(x) c(mean = mean(x), sd = sd(x)))
mean_sd_per_stationtype[, -1] <- round(mean_sd_per_stationtype[, -1], 2)
print(mean_sd_per_stationtype)

# Compute the Pearson correlation matrix
correlation_matrix <- cor(environmental_data, use = "complete.obs", method = "pearson")
# Print the correlation matrix with highlighted values greater than or equal to 0.8
print_correlation <- function(cor_matrix, threshold) {
  col_names <- colnames(cor_matrix)  # Get the column names
  n <- ncol(cor_matrix)  # Number of columns and rows
  
  cat("\t")  # Start with a tab for better alignment
  cat(paste(col_names, collapse = "\t"), "\n")  # Print header row
  
  # Loop through the matrix and print each value with row names
  for (i in 1:n) {
    cat(col_names[i], "\t")  # Print row name
    for (j in 1:n) {
      if (abs(cor_matrix[i, j]) >= threshold && i != j) {
        cat(sprintf("*%.2f*\t", cor_matrix[i, j]))  # Highlight strong correlations
      } else {
        cat(sprintf("%.2f\t", cor_matrix[i, j]))  # Normal print for others
      }
    }
    cat("\n")  # New line for the next row
  }
}
# Use the function to print the correlation matrix with highlighted values
print_correlation(correlation_matrix, 0.8)  # Updated threshold to 0.8

## Community data, PCA
## Abundance, PCA, env as passive
species_data <- dft[, 22:210]  # Selecting species community columns
env_variables <- dft[, 8:21]  # Selecting environmental variables
env_variables$Texture.group <- NULL  # Exclude the 'Texture.group' column

# Perform PCA on species data
pca_res <- PCA(species_data, scale.unit = TRUE, graph = FALSE)

# Standardize environmental variables (mean = 0, sd = 1)
env_scaled <- scale(env_variables)

# Calculate loadings of standardized environmental variables on PCA
loadings <- cor(env_scaled, pca_res$ind$coord)

# Create a data frame for the PCA scores with Station labels and types
pca_scores <- data.frame(pca_res$ind$coord, Station = dft$Stations, StationType = dft$Stationtype, Position = dft$Position, Fjord = dft$Fjord)

# Calculate the scaling factor as described earlier
pca_range_dim1 <- max(pca_scores$Dim.1) - min(pca_scores$Dim.1)
pca_range_dim2 <- max(pca_scores$Dim.2) - min(pca_scores$Dim.2)
max_loading_length <- max(sqrt(loadings[, "Dim.1"]^2 + loadings[, "Dim.2"]^2))
max_arrow_proportion <- 0.3
scaling_factor <- min(pca_range_dim1, pca_range_dim2) * max_arrow_proportion / max_loading_length
loadings_scaled <- loadings * scaling_factor

# Calculate covariance matrix of standardized variables
cov_matrix <- cov(species_data)
# Perform eigen decomposition
eigen_decomp <- eigen(cov_matrix)
# Extract eigenvalues
eigenvalues <- eigen_decomp$values
# Calculate the proportion of variance explained by each principal component
var_explained_by_axis <- eigenvalues / sum(eigenvalues) * 100
# Print the proportion of variance explained by each axis
cat("Proportion of Variance Explained by each PCA Axis:\n")
cat("PCA Axis 1: ", round(var_explained_by_axis[1], 2), "%\n")
cat("PCA Axis 2: ", round(var_explained_by_axis[2], 2), "%\n")



# Define station colors
station_colors <- c("Marine endpoint" = "#000080",    # Navy Blue
                    "Fjord mouth" = "#1E90FF",        # Dodger Blue
                    "Mid fjord transect" = "#87CEEB",  # Sky Blue
                    "Nearshore control" = "#32CD32",   # Lime Green
                    "River Estuary" = "#008000",       # Green
                    "Glacier Influenced" = "#FFD700"  # Gold
)

# Update the plot code to include this information in the axis labels and change the theme
p <- ggplot(data = pca_scores, aes(x = Dim.1, y = Dim.2, color = StationType)) +
  geom_point(size = 3) +  # Plot points colored by 'StationType'
  geom_segment(data = as.data.frame(loadings_scaled), aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2), 
               arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = "blue", alpha = 0.1) +
  geom_text(data = as.data.frame(loadings_scaled), aes(x = Dim.1, y = Dim.2, label = rownames(loadings_scaled)), 
            vjust = -1, color = "blue", alpha = 0.5, size = 3) +
  geom_hline(yintercept = mean(pca_scores$Dim.2), linetype = "dashed", color = "black") +  # Horizontal line at the mean of Dim.2
  geom_vline(xintercept = mean(pca_scores$Dim.1), linetype = "dashed", color = "black") +  # Vertical line at the mean of Dim.1
  xlab(paste("PC1 -", round(var_explained_by_axis[1], 2), "% variance")) +  # X-axis label with variance explained
  ylab(paste("PC2 -", round(var_explained_by_axis[2], 2), "% variance")) +  # Y-axis label with variance explained
  scale_color_manual(values = station_colors, 
                     breaks = names(station_colors)) +  # Apply custom color scale and ensure legend order
  theme_classic() +  # Use classic theme for a clean look
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))  # Add black border

# Display the plot
print(p)

# Calculate centroids for each Fjord
fjord_centroids <- aggregate(cbind(Dim.1, Dim.2) ~ Fjord, data = pca_scores, mean)

# Calculate centroids for each Stationtype
stationtype_centroids <- aggregate(cbind(Dim.1, Dim.2) ~ StationType, data = pca_scores, mean)

# Calculate centroids for each Position
position_centroids <- aggregate(cbind(Dim.1, Dim.2) ~ Position, data = pca_scores, mean)

# Now create the plot with centroids for Fjords, StationTypes, and Position
pfs <- ggplot() +
  geom_point(data = fjord_centroids, aes(x = Dim.1, y = Dim.2), size = 5, shape = 15, color = "blue") +  # Centroids for Fjords
  geom_text(data = fjord_centroids, aes(x = Dim.1, y = Dim.2, label = Fjord), vjust = -1, color = "blue") +  # Labels for Fjord centroids
  geom_point(data = stationtype_centroids, aes(x = Dim.1, y = Dim.2), size = 5, shape = 15, color = "#FFD700") +  # Centroids for StationTypes in gold
  geom_text(data = stationtype_centroids, aes(x = Dim.1, y = Dim.2, label = StationType), vjust = 2, color = "#FFD700") +  # Labels for StationType centroids in gold
  geom_point(data = position_centroids, aes(x = Dim.1, y = Dim.2), size = 5, shape = 15, fill = "white", color = "black") +  # Centroids for Position in white with black border
  geom_text(data = position_centroids, aes(x = Dim.1, y = Dim.2, label = Position), vjust = 2, color = "black") +  # Labels for Position centroids
  geom_segment(data = as.data.frame(loadings_scaled), aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2, label = rownames(loadings_scaled)), 
               arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = "blue", alpha = 0.1) +
  geom_text(data = as.data.frame(loadings_scaled), aes(x = Dim.1, y = Dim.2, label = rownames(loadings_scaled)), 
            vjust = -1, color = "blue") +  # Labels for environmental drivers
  geom_hline(yintercept = mean(pca_scores$Dim.2), linetype = "dashed", color = "black") +
  geom_vline(xintercept = mean(pca_scores$Dim.1), linetype = "dashed", color = "black") +
  xlab(paste("PC1 -", round(var_explained_by_axis[1], 2), "% variance")) +
  ylab(paste("PC2 -", round(var_explained_by_axis[2], 2), "% variance")) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")  # Remove the legend

# Display the plot
print(pfs)


###

# Macrobenthic community general info
# Abundance, top 10, all dataset
abundance_data <- df[, 22:210]
species_abundance <- colSums(abundance_data)
sorted_species <- sort(species_abundance, decreasing = TRUE) # Sorting species by abundance in descending order
top_10_species <- head(sorted_species, 10)
total_abundance <- sum(species_abundance)
percentage_representation <- (top_10_species / total_abundance) * 100
top_10_species_df <- data.frame(
  Species = names(top_10_species),
  Abundance = top_10_species,
  Percentage = percentage_representation
)
print(top_10_species_df)
print(total_abundance)

# Calculate the total abundance for each Stationtype
total_abundance_by_stationtype <- tapply(rowSums(abundance_data), df$Stationtype, sum)
mean_abundance_by_stationtype <- tapply(rowSums(abundance_data), df$Stationtype, mean) * 10
sd_abundance_by_stationtype <- tapply(rowSums(abundance_data), df$Stationtype, sd) * 10

# Calculate number of taxa per station
taxa_per_station <- apply(abundance_data, 1, function(x) sum(x > 0))

# Combine with Stationtype
taxa_per_station_type <- aggregate(taxa_per_station, by = list(df$Stationtype), FUN = sum)

# Get the count of stations per Stationtype
station_counts <- table(df$Stationtype)

# Divide total number of taxa by the count of stations for each Stationtype
mean_taxa_per_station_type <- taxa_per_station_type$x / station_counts[taxa_per_station_type$Group.1]

# Calculate standard deviation of taxa per station type
standard_deviation_taxa_per_station_type <- aggregate(taxa_per_station, by = list(df$Stationtype), FUN = sd)$x

# Combine mean and standard deviation into a data frame
mean_std_taxa_per_station_type <- data.frame(
  Stationtypeabundance = total_abundance_by_stationtype,
  StationType = taxa_per_station_type$Group.1,
  MeanTaxa = mean_taxa_per_station_type,
  StdDevTaxa = standard_deviation_taxa_per_station_type,
  MeanAbundance = mean_abundance_by_stationtype,
  SdAbundance = sd_abundance_by_stationtype
)

print(mean_std_taxa_per_station_type)


## RDA
#Normalize and Standardization of environmental factors
#plot the relationships
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
pairs(dft[,c(7:15, 17, 19)], panel=panel.smooth, diag.panel=panel.hist)

# Set "Stations" column as row names
rownames(dft) <- dft$Stations

#different scales, standardised
t.env <- decostand(dft[, c(7:11, 13, 19)] - min(dft[, c(7:11, 13, 19)]) + 1, method = "log")   #the ones that are not standardised, 9,13,18(SI)- negative values
env.stand <-cbind(t.env,dft[,c(12,14:15,17)])

#RDA
spec.rda <- rda(df_hel ~ ., env.stand)
summary(spec.rda)

#plotRDA
model <- ordiplot(spec.rda, type = "none", scaling = 2, cex=10, xlab = "RDA1 (X%)", ylab = "RDA2 (X%)", cex.lab=1.25)
points(spec.rda, col="black", cex=1)
points(spec.rda, dis="sp", col="darkgrey",  cex=1)
text(spec.rda, dis="sp", col="darkgrey", cex=1)
text(spec.rda, dis="bp", col="black")




# Select relevant columns
env_data <- dft[, c(7:15, 17, 19)]
species_data <- dft[, 20:208]

# RDA
rda <- rda(species_data ~ ., data = env_data)

# Perform forward selection
# Define scope variables
scope_vars <- colnames(env_data)

# Perform forward selection with double-stopping criterion
forward_sel <- ordiR2step(rda, scope = scope_vars, perm.max = 9999)

forward_sel <- ordiR2step(rda, perm.max = 9999)

scope <- formula(~ .)  # Specify the model space
forward_sel <- ordiR2step(rda, scope = scope, data = env_data, perm.max = 9999)

# Variation partitioning
var_part <- varpart(rda)

# Output results
print("RDA:")
print(summary(forward_sel))
print("Variation partitioning:")
print(var_part)
##

## PCA
## Clustering 
## UPMGA, incl mult bootstrap (9999) and species fidelity 

macrodata <- dft[, (22:210)] # Subset the data for species abundance
rownames(macrodata) <- dft$Stations # Assign station names as row names
dist_matrix <- dist(macrodata, method = "euclidean") # Compute distance matrix
upgma_cluster <- hclust(dist_matrix, method = "average")  # Perform hierarchical clustering
dend <- as.dendrogram(upgma_cluster) # Convert hclust object to dendrogram
dend <- color_branches(dend, k = 9) # Colour clustering

par(mar = c(5, 4, 4, 2) + 0.1)  # Adjust margin to make room for station names
plot(dend, main = "UPGMA Clustering Dendrogram", xlab = "Distance", ylab = "", 
     sub = "", cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.5)  # Plot dendrogram
labels_colors(dend) <- 1  # Change colour of labels

result <- pvclust(t(macrodata), method.dist="euclidean", method.hclust="average", nboot=9999) # Perform multiscale Bootstrap Resampling
plot(result) 
pvrect(result)  # Add rectangles to identify clusters
summary(result) # Summary of clustering result

# Clustering description#

cluster_labels <- cutree(upgma_cluster, k = 9)
station_clusters <- data.frame(Station = rownames(macrodata), Cluster = cluster_labels)
for (i in 1:9) {
  cat("Cluster", i, ":", paste(station_clusters$Station[station_clusters$Cluster == i], collapse = ", "), "\n")
}

## NMDS plot
bc_dissimilarity <- vegdist(macrodata, method = "bray")
nmds_result <- metaMDS(bc_dissimilarity)
cluster_assignments <- cutree(as.hclust(result$hclust), k = 9)  # Adjusted k, number of clusters from the boots
plot(nmds_result$points, type = "n", main = "NMDSPlot", xlab = "NMDS1", ylab = "NMDS2") # colouring plot by cluster
text(nmds_result$points, labels = rownames(nmds_result$points), cex = 0.7, col = cluster_assignments)
mtext(paste("Stress Level:", round(nmds_result$stress, 4)), side = 1, line = -22, adj = 1, cex = 0.7)

# IndVal
# Make sure your Stationtype column is a factor
dft$Stationtype <- as.factor(dft$Stationtype)

# Subset your data for the species abundance
species_data <- dft[, 22:210]

# Run the indval function
result <- indval(species_data, dft$Stationtype, permutation =9999)

# Print the result
print(result)

# Convert the result to a data frame
result_df <- as.data.frame(result$indval)

# Add a column for the species names
result_df$Species <- rownames(result_df)

# Reshape the data from wide to long format
long_df <- melt(result_df, id.vars = "Species", variable.name = "Stationtype", value.name = "IndVal.g")

# Find the top 5 species with the highest IndVal.g for each Stationtype
top_species <- long_df %>% group_by(Stationtype) %>% top_n(5, IndVal.g) %>% arrange(Stationtype, desc(IndVal.g))

# Print the result
print(top_species)



###

##Taxonomic indices

# Calculate diversity indices separated by Stations
stations <- unique(df$Stations)
for (station in stations) {
  cat("Station:", station, "\n")
  
  # Subset data for each Station
  subset_data <- df[df$Stations == station, ]
  
  # Check if there are any rows to calculate diversity indices
  if (nrow(subset_data) == 0) {
    cat("No data for this Station.\n\n")
    next  # Skip to the next Station
  }
  
  # Calculate Shannon diversity index (H???(log e))
  shannon_index_loge <- diversity(subset_data[, 22:210], index = "shannon")
  
  # Calculate the number of species (non-zero counts)
  num_species <- sum(subset_data[, 22:210] > 0)
  
  # Calculate Pielou's evenness index
  pielou_index <- shannon_index_loge / log(num_species)
  
  # Calculate Shannon diversity index (H???(log 2))
  shannon_index_log2 <- diversity(subset_data[, 22:210], index = "shannon", base = 2)
  
  # Calculate mean for each index
  pielou_mean <- mean(pielou_index)
  
  shannon_loge_mean <- mean(shannon_index_loge)
  
  shannon_log2_mean <- mean(shannon_index_log2)
  
  # Format mean for each index to two decimal places
  pielou_mean <- round(pielou_mean, 2)
  
  shannon_loge_mean <- round(shannon_loge_mean, 2)
  
  shannon_log2_mean <- round(shannon_log2_mean, 2)
  
  # Print mean for each index
  cat("Pielou's Evenness Index (Mean):", pielou_mean, "\n")
  cat("Shannon Diversity Index (H???(log e)) (Mean):", shannon_loge_mean, "\n")
  cat("Shannon Diversity Index (H???(log 2)) (Mean):", shannon_log2_mean, "\n")
  
  cat("\n")
}

# Calculate diversity indices separated by Stationtype
stationtypes <- unique(dft$Stationtype)
for (st in stationtypes) {
  cat("Stationtype:", st, "\n")
  
  # Subset data for each Stationtype
  subset_data <- dft[dft$Stationtype == st, ]
  
  # Check if there are any rows to calculate diversity indices
  if (nrow(subset_data) == 0) {
    cat("No data for this Stationtype.\n\n")
    next  # Skip to the next Stationtype
  }
  
  # Calculate Shannon diversity index (H???(log e))
  shannon_index_loge <- diversity(subset_data[, 22:210], index = "shannon")
  
  # Calculate the total number of individuals
  total_individuals <- sum(subset_data[, 22:210])
  
  # Calculate the number of species (non-zero counts)
  num_species <- sum(subset_data[, 22:210] > 0)
  
  # Calculate Pielou's evenness index
  pielou_index <- shannon_index_loge / log(num_species)
  
  # Calculate Shannon diversity index (H???(log 2))
  shannon_index_log2 <- diversity(subset_data[, 22:210], index = "shannon", base = 2)
  
  # Calculate mean and standard deviation for each index
  pielou_mean <- mean(pielou_index)
  pielou_sd <- sd(pielou_index)
  
  shannon_loge_mean <- mean(shannon_index_loge)
  shannon_loge_sd <- sd(shannon_index_loge)
  
  shannon_log2_mean <- mean(shannon_index_log2)
  shannon_log2_sd <- sd(shannon_index_log2)
  
  # Format mean and standard deviation for each index to two decimal places
  pielou_mean <- round(pielou_mean, 2)
  pielou_sd <- round(pielou_sd, 2)
  
  shannon_loge_mean <- round(shannon_loge_mean, 2)
  shannon_loge_sd <- round(shannon_loge_sd, 2)
  
  shannon_log2_mean <- round(shannon_log2_mean, 2)
  shannon_log2_sd <- round(shannon_log2_sd, 2)
  
  # Print mean and standard deviation for each index
  cat("Pielou's Evenness Index (Mean ?? SD):", pielou_mean, "??", pielou_sd, "\n")
  cat("Shannon Diversity Index (H???(log e)) (Mean ?? SD):", shannon_loge_mean, "??", shannon_loge_sd, "\n")
  cat("Shannon Diversity Index (H???(log 2)) (Mean ?? SD):", shannon_log2_mean, "??", shannon_log2_sd, "\n")
  
  cat("\n")
}

##CWM
# Extract the relevant columns for taxa abundance from 'dft'
abundance_data <- dftab[, 22:155]
abundance_data <- as.data.frame(abundance_data)

# Convert raw abundance data to relative abundances
relative_abundance_data <- sweep(abundance_data, 1, rowSums(abundance_data), FUN="/")

# Convert dftrait from a tibble to a standard data frame
dftrait <- as.data.frame(dftrait)

# Set the first column of dftrait as row names
rownames(dftrait) <- dftrait[, 1]

# Remove the first column from dftrait as it's now set as row names
dftrait <- dftrait[, -1]

# Ensure that the taxa order matches in both abundance and trait dataframes
if(!all(colnames(abundance_data) %in% rownames(dftrait))) {
  stop("Some taxa in abundance data are not present in trait data. Please align them.")
}

# Align dftrait to match the order of taxa in abundance_data
dftrait_aligned <- dftrait[colnames(abundance_data), , drop = FALSE]
dftrait_aligned[is.na(dftrait_aligned)] <- 0

# Calculate Community Weighted Means (CWM)
abundance_matrix <- as.matrix(relative_abundance_data)
trait_matrix <- as.matrix(dftrait_aligned)
cwm <- abundance_matrix %*% trait_matrix

# Convert the resulting matrix back to a dataframe for easier handling
cwm_df <- as.data.frame(cwm)
rownames(cwm_df) <- rownames(abundance_data)  # Assigning row names from the abundance data (station names)
colnames(cwm_df) <- colnames(trait_matrix)  # Assigning column names from the trait data (trait names)

# Add Factor Information to CWM Dataframe
cwm_df$Stationtype <- dft$Stationtype
cwm_df$Position <- dft$Position
cwm_df$Fjord <- dft$Fjord

# Calculate Mean CWM Values for Each Factor
# For Stationtype
mean_cwm_by_stationtype <- cwm_df %>%
  group_by(Stationtype) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# For Position
mean_cwm_by_position <- cwm_df %>%
  group_by(Position) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# For Fjord
mean_cwm_by_fjord <- cwm_df %>%
  group_by(Fjord) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

#CWM plot, radar, graveyard data
#CWM Weighted Averaging, single value per trait, graveyard data
#Correlation among traits to check
# Combine all trait modality columns into one vector for simplicity
trait_columns <- c("AH1", "AH2", "AH3", "AH4", "AH5", "AH6", 
                   "FH1", "FH3", "FH5", "FH6", "FH7", "FH8", "FH9", "FH10", "FH11",
                   "NS1", "NS2", "NS3", "NS4", "NS5", "NS6",
                   "LD1", "LD2", "LD3",
                   "RT1", "RT2", "RT3", "RT4",
                   "LT1", "LT2", "LT3")


# Calculate Spearman's correlation matrix for the traits
cor_matrix <- cor(cwm_df[, trait_columns], method = "spearman")

# Find correlations that are greater than or equal to 0.55
threshold <- 0.55
high_cor_mask <- abs(cor_matrix) >= threshold

# Applying the mask and setting non-significant values to NA for clarity
significant_cor_matrix <- cor_matrix * high_cor_mask
significant_cor_matrix[!high_cor_mask] <- NA

# View the significant correlation matrix
print(significant_cor_matrix)

# Extract and list significant correlations
significant_pairs <- which(abs(cor_matrix) >= threshold, arr.ind = TRUE)
cor_list <- setNames(object = cor_matrix[significant_pairs], nm = apply(significant_pairs, 1, function(idx) {
  paste(colnames(cor_matrix)[idx[1]], colnames(cor_matrix)[idx[2]], sep = " and ")
}))

# Remove duplicates (since matrix is symmetric)
cor_list <- cor_list[!duplicated(t(apply(significant_pairs, 1, sort)))]

# Remove pairs with the same initial acronyms
cor_list <- cor_list[sapply(names(cor_list), function(name) {
  # Split the trait names and get the first two characters of each
  traits <- strsplit(name, " and ")[[1]]
  substr(traits[1], 1, 2) != substr(traits[2], 1, 2)
})]

# Print the list of significant correlations
print(cor_list)


#Functional analysis
## Rao Index Calculation
# Exclude columns 7 to 16 from trait_matrix (mobility and body form)
trait_matrix_filtered <- trait_matrix[, -c(7:16)]

# Calculate the trait distance matrix
trait_distance_matrix <- as.dist(dist(trait_matrix, method = "euclidean"))
trait_distance_matrix <- as.matrix(trait_distance_matrix)

# Normalize the trait distance matrix to have a maximum of 1
trait_distance_matrix <- trait_distance_matrix / max(trait_distance_matrix)

# Calculate Rao's index for each station using the normalized trait distance matrix
rao_indices <- numeric(nrow(abundance_matrix))

# Extract station identifiers and types from the dftab data frame
station_names <- dftab$Stations  # Adjust this if the column name or position differs
station_types <- dftab$Stationtype  # Ensure this matches your data frame column name for station types

# Calculate Rao's index for each station
for (i in 1:nrow(abundance_matrix)) {
  current_abundance <- abundance_matrix[i, , drop = FALSE]
  weighted_trait_dist <- current_abundance %*% trait_distance_matrix %*% t(current_abundance)
  rao_indices[i] <- sum(weighted_trait_dist) / (sum(current_abundance) ^ 2)  # Ensure normalization
}

# Combine station names, types, and their corresponding Rao's indices into one data frame
rao_results <- data.frame(Station = station_names, Station_Type = station_types, Rao_Index = rao_indices)

# Calculate mean and standard deviation of Rao's indices per Station Type
summary_stats <- rao_results %>%
  group_by(Station_Type) %>%
  summarise(
    Mean_Rao_Index = mean(Rao_Index, na.rm = TRUE),
    SD_Rao_Index = sd(Rao_Index, na.rm = TRUE)
  )

# Print the Rao's indices with station identifiers and types
print(rao_results)
# Print summary statistics for each station type
print(summary_stats)

##Functional Redundancy
# Calculate Shannon-Wiener index for each station
shannon_indices <- apply(abundance_matrix, 1, function(x) {
  # Replace NAs with zeros for calculation (if any NAs present)
  x[is.na(x)] <- 0
  diversity(x, index = "shannon")
})

# Convert natural log base to log base e
shannon_indices <- shannon_indices / log(exp(1))

# Add Shannon indices to the existing Rao results data frame
rao_results$Shannon_Index <- shannon_indices

# Calculate Functional Redundancy as the ratio of Rao to Shannon-Wiener index
rao_results$Functional_Redundancy <- rao_results$Rao_Index / rao_results$Shannon_Index

# Optionally, summarize the results by Station_Type if desired
summary_stats <- rao_results %>%
  dplyr::group_by(Station_Type) %>%
  dplyr::summarise(
    Mean_Rao_Index = mean(Rao_Index, na.rm = TRUE),
    SD_Rao_Index = sd(Rao_Index, na.rm = TRUE),
    Mean_Shannon_Index = mean(Shannon_Index, na.rm = TRUE),
    Mean_Functional_Redundancy = mean(Functional_Redundancy, na.rm = TRUE)
  )

# Print the extended Rao's indices with station identifiers, types, and functional redundancy
print(rao_results)
# Print updated summary statistics for each station type
print(summary_stats)

##Functional Richness
#functional richness at each location using the fd_fric() function

fric = fd_fric(trait_matrix, abundance_matrix, stand = FALSE)

head(fric)

##PERMANOVA


help(fd_fric)
data(traits_birds)
fd_fric(traits_birds)

