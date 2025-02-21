# Load required packages for trophic analysis, plotting, data reshaping, and data manipulation.

library(tRophicPosition) # For trophic position estimation using isotope data.
library(ggplot2)         # For data visualization.
library(reshape2)        # For reshaping data.
library(plyr)            # For data manipulation and summarization.
library(dplyr)           # to summarize model output

# Load the stable isotope dataset from a CSV file.
dat_full <- read.csv("archive_data_2007_2008.csv")
# trim out the water, and sediment data
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

# Retrieve trophic discrimination factors (TDF) for both δ13C and δ15N based on McCutchan.
TDF_values <- TDF(author = "McCutchan", element = "both", type= "all")

# Extract isotope data using species as the consumer grouping.
TinicumList <- extractIsotopeData(dat,
                               consumersColumn = "Species",
                               consumer = consumers,
                               b1 = "algae",
                               b2 = "plant",
                               baselineColumn = "Class",
                               groupsColumn = NULL, # this was run previously with groups, the post per species were not different
                               deltaC = TDF_values$deltaC,
                               deltaN = TDF_values$deltaN)

# For each River-specific dataset in TinicumList, generate diagnostic plots and summaries.
for (dataset in TinicumList) {
  plot(dataset)            # Generate a diagnostic plot for each dataset.
  print(summary(dataset))  # Print summary statistics to the console.
}

speciesTP <- multiSpeciesTP(TinicumList,
                            lambda = 1,           # Source trophic position.
                            adapt = 50000,      # Adaptation iterations.
                            n.iter = 500000,     # Total iterations.
                            burnin = 100000,       # Burnin period.
                            n.chains = 3,         # Number of MCMC chains.
                            thin = 100,
                            model = "twoBaselines")  # Model type.


# Extract the posterior TP draws from the speciesTP object.
TP_list <- speciesTP$TPs  # Each element is a vector of TP draws for a consumer

# Combine the draws into a single data frame with a "consumer" column.
TP_df <- do.call(rbind, lapply(names(TP_list), function(nm) {
  data.frame(consumer = nm, TP = TP_list[[nm]])
}))

# Compute the median TP for each consumer and arrange in descending order.
order_df <- TP_df %>%
  group_by(consumer) %>%
  summarize(medTP = median(TP)) %>%
  arrange(desc(medTP))

# Reorder the consumer factor in the combined data frame.
TP_df$consumer <- factor(TP_df$consumer, levels = order_df$consumer)

# Define a named vector for the new labels.
consumer_labels <- c("channel.catfish.2b" = "Channel Catfish",
                     "mummichug.2b" = "Mummichog",
                     "banded.killfish.2b" = "Banded Killifish",
                     "white.sucker.2b" = "White Sucker",
                     "amphipods.2b" = "Amphipoda spp.",
                     "corbiculae.2b" = "Corbicula fluminea")

# Create the boxplot with custom x-axis labels and annotation.
p = ggplot(TP_df, aes(x = consumer, y = TP, fill = consumer)) +
  geom_boxplot(outlier.shape = NA) +         # Boxplot without outliers
  theme_bw() +
  labs(title = "Posterior Distributions of Trophic Positions by Consumer",
       x = "Consumer",
       y = "Trophic Position") +
  scale_x_discrete(labels = consumer_labels) +  # Apply custom labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  annotate("text", x = Inf, y = Inf, label = "Data from Ashley et. al., 2012",
           hjust = 1.1, vjust = 1.5, size = 3, fontface = "italic") +
  scale_y_continuous(breaks = 1:8, limits = c(1, 8))

ggsave("TrophicPositions_Ashley_2012.png", plot = p, width = 6, height = 4, units = "in", dpi = 300)


# These pairwise trophic position comparisons were used previously when the isotope data was
# grouped by site. None of them indicated any statistical difference between sites
# the data setup was altered to remove the grouping variable at Line 49 
# I would guess this is because the number of samples is so low.

# pairwiseTP <- pairwiseComparisons(speciesTP$TPs[c(1,14)], print = TRUE) # this indexing is no longer logicial as the model is differnet
