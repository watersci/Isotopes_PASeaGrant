# Load required packages for trophic analysis, plotting, data reshaping, and data manipulation.
library(tRophicPosition) # For trophic position estimation using isotope data.
library(ggplot2)         # For data visualization
library(reshape2)        # For reshaping data.
library(plyr)            # For data manipulation and summarization
library(dplyr)           # to summarize model output
library(patchwork)       # to plot in more than one frame
library(coda)            # For MCMC diagnostics
library(tidyr)           # to summarize model output


# Load the stable isotope dataset from a CSV file. Citation for Data below
# Ashley, J. T. F., Vasquez, M. A., Zelanko, P., McKinley, E., Schafer, M.,
# Zaoudeh, L., Horwitz, R., Stapleton, H. M., & Velinsky, D. J. (2012). Trophic transfer of 
# polybrominated diphenyl ethers and polychlorinated biphenyls in a tidal freshwater marsh.
# Chemistry and Ecology, 28(4), 305–325.

dat_full <- read.csv("archive_data_2007_2008.csv")

# trim out the water, and sediment data, we might need this late, but not for now
dat = dat_full[-c(which(dat_full$Species == "water"),which(dat_full$Species == "sediment")),]

# Filter 'dat' to include only consumer species with more than 3 samples.
dat <- dat %>%
  group_by(Species) %>%
  filter(
    (all(Class %in% c("algae", "plant"))) | (n() > 3)
  ) %>%
  ungroup() %>%
  as.data.frame()

# Summarize the data by River and Species, calculating means and standard deviations 
# for δ13C and δ15N values.
data.sum<-ddply(dat, c("River", "Species"), summarise,
                d13Cmn=mean(d13C), 
                d13Csd=sd(d13C),
                d15Nmn=mean(d15N),
                d15Nsd=sd(d15N))

#### Define trophic position ####

# Define the baseline species used for the mixing model.
baseline_species <- c("algae", "plant")

# Identify consumer species by excluding baseline species.
consumers <- unique(dat$Species[!(dat$Class %in% baseline_species)])

# Retrieve trophic discrimination factors (TDF) for both δ13C and δ15N based on 
# McCutchan Jr, J. H., Lewis Jr, W. M., Kendall, C., & McGrath, C. C. (2003). 
# Variation in trophic shift for stable isotope ratios of carbon, nitrogen, 
# and sulfur. Oikos, 102(2), 378–390.
TDF_values <- TDF(author = "McCutchan", element = "both", type= "all")

# Extract isotope data using species as the consumer grouping, this is the 
# needed format for the tRophicPosition package 
TinicumList <- extractIsotopeData(dat,
                               consumersColumn = "Species",
                               consumer = consumers,
                               b1 = "algae",
                               b2 = "plant",
                               baselineColumn = "Class",
                               groupsColumn = NULL, 
                               deltaC = TDF_values$deltaC,
                               deltaN = TDF_values$deltaN)

# For each River-specific dataset in TinicumList, generate diagnostic plots and summaries.
# for (dataset in TinicumList) {
#   plot(dataset)            # Generate a diagnostic plot for each dataset.
#   print(summary(dataset))  # Print summary statistics to the console.
# }

speciesTP <- multiSpeciesTP(TinicumList,
                            lambda = 1,           # Source trophic position.
                            adapt = 50000,      # Adaptation iterations.
                            n.iter = 50000,     # Total iterations.
                            burnin = 10000,       # Burnin period.
                            n.chains = 3,         # Number of MCMC chains.
                            thin = 10,
                            model = "twoBaselinesFull")  # Model type.

# look for well mixed chains and histograms to QC model
# for (consumer in names(speciesTP$multiSpeciesTP)) {
#   mcmc_samples <- speciesTP$multiSpeciesTP[[consumer]]$samples[[1]]
#   plot(mcmc_samples,  xlab = consumer)
# }
# All look good, let's examine the results!

# Extract the posterior TP draws from the speciesTP object.
TP_list <- speciesTP$TPs  # Each element is a vector of TP draws for a consumer

# Combine the draws into a single data frame with a "consumer" column.
TP_df <- do.call(rbind, lapply(names(TP_list), function(nm) {
  data.frame(consumer = nm, TP = TP_list[[nm]])
}))

# Compute the median TP for each consumer and arrange in descending order
# We do this so the plots can be ordered in a way that looks mice
order_df <- TP_df %>%
  group_by(consumer) %>%
  summarize(medTP = median(TP)) %>%
  arrange(desc(medTP))

# Reorder the consumer factor in the combined data frame.
TP_df$consumer <- factor(TP_df$consumer, levels = order_df$consumer)

# Define a named vector for the new labels so the plots look good
consumer_labels <- c("channel.catfish.2bf" = "Channel Catfish",
                     "mummichug.2bf" = "Mummichog",
                     "banded.killfish.2bf" = "Banded Killifish",
                     "white.sucker.2bf" = "White Sucker",
                     "amphipods.2bf" = "Amphipoda spp.",
                     "corbiculae.2bf" = "Corbicula fluminea")

# Create the boxplot of Trophic Positions
p = ggplot(TP_df, aes(x = consumer, y = TP, fill = consumer)) +
  geom_boxplot(outlier.shape = NA) +         # Boxplot without outliers
  theme_bw() +
  labs(title = "Posterior Distributions of Trophic Positions by Consumer",
       x = "",
       y = "Trophic Position") +
  scale_x_discrete(labels = consumer_labels) +  # Apply custom labels
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  ) +
  annotate("text", x = Inf, y = Inf, label = "Data from Ashley et. al., 2012",
           hjust = 1.1, vjust = 1.5, size = 3, fontface = "italic") +
  scale_y_continuous(breaks = 1:8, limits = c(1, 8))

# ggsave("TrophicPositions_Ashley_2012.png", plot = p, width = 6, height = 4, units = "in", dpi = 300)

# Extract the posterior alphas draws from the speciesTP object.
Alphas_list <- speciesTP$Alphas  # Each element is a vector of TP draws for a consumer

# Combine the draws into a single data frame with a "consumer" column just like above
Alphas_list_df <- do.call(rbind, lapply(names(Alphas_list), function(nm) {
  data.frame(consumer = nm, alpha = Alphas_list[[nm]])
}))

# again this is just to order the 
Alphas_list_df$consumer <- factor(Alphas_list_df$consumer, levels = order_df$consumer)

# boxplot of the Alphas
w <- ggplot(Alphas_list_df, aes(x = consumer, y = alpha, fill = consumer)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  labs(title = "Posterior Distributions Proportion Diet Derives from Algae vs. Plants",
       x = "Consumer",
       y = "Proportion Derived from Algae") +
  scale_x_discrete(labels = consumer_labels) +
  annotate("text", x = -Inf, y = -Inf, 
           label = "Data from Ashley et. al., 2012",
           hjust = -0.1, vjust = -0.5, size = 3, fontface = "italic") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )


# Combine plots side by side and add tags 'a' and 'b'
combined_plot <- p / w + plot_annotation(tag_levels = "a")

# Save the combined plot as a PNG
# ggsave("TrophicPositions_Combined.png", plot = combined_plot, width = 6, height = 6.5, units = "in", dpi = 300)

median_TP <- TP_df %>%
  group_by(consumer) %>%
  summarize(med_TP = median(TP, na.rm = TRUE))

median_alpha <- Alphas_list_df %>%
  group_by(consumer) %>%
  summarize(med_alpha = median(alpha, na.rm = TRUE))

median_summary <- left_join(median_TP, median_alpha, by = "consumer")
print(median_summary)

median_summary_sentences <- median_summary %>%
  # Map the internal consumer code to the full label.
  mutate(consumer_label = consumer_labels[as.character(consumer)]) %>%
  rowwise() %>%
  # Create a sentence using the rounded values.
  mutate(sentence = paste0("For ", consumer_label, 
                           ", the median trophic position is ", round(med_TP, 2), 
                           " and the median alpha is ", round(med_alpha, 3), ".")) %>%
  ungroup()

# Print each sentence on a new line.
median_summary_sentences$sentence %>% 
  paste(collapse = "\n") %>%
  cat()


#### Use of older-simpiler techniques ####
# Based off Vander Zanden et al. (1997) citation
# Vander Zanden, M. J., Cabana, G., & Rasmussen, J. B. (1997). Comparing trophic 
# position of freshwater fish calculated using stable nitrogen isotope ratios 
# (δ15N) and literature dietary data. Canadian Journal of Fisheries and 
# Aquatic Sciences, 54(5), 1142–1158.

baseline_data <- dat_full %>% 
  filter(Species == "corbiculae")

# Calculate the mean δ15N for the baseline for each River.
baseline_means <- baseline_data %>%
  group_by(River) %>%
  summarize(d15N_baseline = mean(d15N, na.rm = TRUE))

# Subset consumer data by removing water, sediment, and the baseline organism.
consumer_data <- dat_full %>% 
  filter(!Species %in% c("corbiculae")) %>%
  filter(!Class %in% c("algae","plant","filter","sediment"))

# Merge consumer data with baseline_means by River
merged_data <- merge(consumer_data, baseline_means, by = "River")

# fix the names, the original data have rough names, here we clean them up
merged_data <- merged_data %>%
  mutate(Species = recode(Species,
                          "anros"                  = "American eel",
                          "darter"                 = "Darter",
                          "banded killfish"        = "Banded Killifish",
                          "eastern silvery minnow" = "Eastern Silvery Minnow",
                          "white sucker"           = "White Sucker",
                          "isopoda"                = "Isopoda spp.",
                          "amphipods"              = "Amphipoda spp.",
                          "channel catfish"        = "Channel Catfish",
                          "crayfish"               = "Crayfish",
                          "mummichug"              = "Mummichog",
                          "worms"                  = "Worms (various)",
                          "Gizzard Shad"           = "Gizzard Shad",
                          "brown bullhead catfish" = "Brown Bullhead Catfish",
                          "Pumpkinseed"            = "Pumpkinseed",
                          "white catfish"          = "White Catfish",
                          "White Perch"            = "White Perch",
                          "red breast"             = "Red Breast",
                          "small prey fish"        = "Small Prey Fish"
  ))

# same here, just capitialize these to make the plot better
merged_data <- merged_data %>%
  mutate(River = recode(River,
                          "creek"                  = "Creek",
                          "marsh"                 = "Marsh"))

# Calculate trophic position (TP) using the Vander Zanden equation.
# Here, 3.4 is used as the per-trophic-level enrichment in δ15N and the baseline trophic level is set to 2.
merged_data <- merged_data %>%
  mutate(TP = ((d15N - d15N_baseline) / 3.4) + 2) %>%
  ungroup()

ordered_species <- merged_data %>%
  group_by(Species) %>%
  summarize(mean_TP = mean(TP, na.rm = TRUE),
            sd_TP = sd(TP, na.rm = TRUE),
            n = n()) %>%
  arrange(desc(mean_TP)) %>%
  select(Species) %>%
  as.vector() %>%
  unlist() %>%
  as.character()

TP_summary <- merged_data %>%
  group_by(Species, River) %>%
  summarize(mean_TP = mean(TP, na.rm = TRUE),
            sd_TP = sd(TP, na.rm = TRUE),
            n = n()) %>%
  ungroup() %>%
  mutate(Species = factor(Species, levels = ordered_species))

# need this to add the data reference in the plot
annotation_df <- data.frame(
  River = levels(factor(TP_summary$River))[1],
  x = Inf,
  y = Inf,
  label = "Data from Ashley et. al., 2012"
)

# box plots by river of observations
q = ggplot(TP_summary, aes(x = mean_TP, y = Species)) +
  geom_col(fill = "grey", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(xmin = mean_TP - sd_TP, xmax = mean_TP + sd_TP),
                width = 0.2,
                size = 1,
                position = position_dodge(width = 0.8)) +
  facet_wrap(~ River, scales = "fixed") +
  scale_y_discrete(drop = FALSE) +
  coord_cartesian(xlim = c(1, NA)) +  # Set lower limit to 1
  theme_bw() +
  labs(title = "Mean Trophic Position by Species per River",
       x = "Mean Trophic Position",
       y = "Species") +
  geom_text(data = annotation_df, 
            aes(x = x, y = y, label = label),
            hjust = 1.1, vjust = 1.5, size = 3, fontface = "italic") +
  theme(axis.text.y = element_text(angle = 30, hjust = 1),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

# ggsave("TP_Nonly_Ashley_2012.png", plot = q, width = 6.5, height = 6, units = "in", dpi = 300)

TP_summary_wide <- data.frame(Species = unique(TP_summary$Species), Creek_n = NA, Marsh_n = NA)
for(spp in TP_summary_wide$Species){
  if(length(intersect(which(TP_summary$Species == spp),which(TP_summary$River == "Creek"))) > 0){
  TP_summary_wide$Creek_n[which(TP_summary_wide$Species == spp)] = TP_summary$n[intersect(which(TP_summary$Species == spp),which(TP_summary$River == "Creek"))]
  } else {
  TP_summary_wide$Creek_n[which(TP_summary_wide$Species == spp)] = 0
  }
  
  if(length(intersect(which(TP_summary$Species == spp),which(TP_summary$River == "Marsh"))) > 0){
    TP_summary_wide$Marsh_n[which(TP_summary_wide$Species == spp)] = TP_summary$n[intersect(which(TP_summary$Species == spp),which(TP_summary$River == "Marsh"))]
  } else {
    TP_summary_wide$Marsh_n[which(TP_summary_wide$Species == spp)] = 0
  }
}

# write.csv(TP_summary_wide, "Sample_n_by_River.csv", row.names = FALSE)
