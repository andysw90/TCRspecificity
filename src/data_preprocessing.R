# Load libraries
library(ggplot2)

# Load data
fname <- "data/vdjdb_CDR3b_HLAA02.tsv"
data_raw <- readLines(fname)
data_raw_split <- strsplit(data_raw,split = "\t")
data_raw_split_lengths <- sapply(data_raw_split,length)
data_in_raw <- data.frame(do.call(rbind,data_raw_split),stringsAsFactors = FALSE)
colnames(data_in_raw) <- data_in_raw[1,]
data_in <- data_in_raw[-1,]

# Remove duplicate CDR-epitope pairings
duplicate_rows <- which(duplicated(paste0(data_in$Epitope,"_",data_in$CDR3)))
data_in <- data_in[-duplicate_rows,]

# Remove very long/short CDR3b's
cdr3_lengths <- nchar(data_in$CDR3)
proportion_of_points_per_length <- table(cdr3_lengths) / length(data_in$CDR3)
top_lengths <- as.numeric(names(which(proportion_of_points_per_length>0.01)))
good_length_indices <- which(cdr3_lengths %in% top_lengths)
data_in <- data_in[good_length_indices,]

barplot(proportion_of_points_per_length,xlab = "CDR3b length",ylab = "Proportion of dataset")

# Take only the CDR3b entries
data_in <- data_in[which(data_in$Gene=="TRB"),]



# Count and plot data points per epitope
epitope_counts <- sort(table(data_in$Epitope))
epitope_counts <- epitope_counts[which(epitope_counts>10)]
epitope_counts_df <- data.frame(Epitope = names(epitope_counts),
                                Count = as.numeric(epitope_counts),stringsAsFactors = FALSE)
epitope_counts_df$Epitope <- factor(epitope_counts_df$Epitope,levels = as.character(epitope_counts_df$Epitope))
ggplot(epitope_counts_df,aes(x = Epitope, y = Count)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  ylab("Number of unique TCRs")

epitope_counts_df_top <- epitope_counts_df[which(epitope_counts_df$Count>100),]
epitope_counts_df_top <- epitope_counts_df_top[sort.int(epitope_counts_df_top$Count,
                                                        index.return = TRUE,
                                                        decreasing = TRUE)$ix,]
epitope_counts_df_top$Epitope <- factor(epitope_counts_df_top$Epitope,levels = unique(epitope_counts_df_top$Epitope))
ggplot(epitope_counts_df_top,aes(x = Epitope, y = Count)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  ylab("Number of unique TCRs")

# Subset on just top epitopes
top_epitopes <- as.character(epitope_counts_df$Epitope[which(epitope_counts_df$Count>500)])
top_epitopes_dataIndices <- unlist(sapply(top_epitopes,function(this_epitope) which(data_in$Epitope==this_epitope)))
data_in <- data_in[top_epitopes_dataIndices,]

# Extract just columns of interest
data_out <- data.frame(CDR3 = data_in$CDR3,
                       Epitope = data_in$Epitope)
write.csv(data_out,"data/vdjdb_CDR3b_HLAA02_processed.csv",row.names = FALSE)

