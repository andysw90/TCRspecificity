# Load packages
library(mltools)
library(data.table)
library(reshape2)
library(xgboost)
library(cowplot)
library(Biostrings)
library(Rtsne)
library(ggplot2)
library(dbscan)
library(cowplot)
library(ggseqlogo)
data(BLOSUM62)

# Load data and some parameters
xgboost_rounds <- 5000
training_size <- 0.8

fname <- "data/vdjdb_CDR3b_HLAA02_processed.csv"
data_in <- read.csv(fname,stringsAsFactors = FALSE)

AAs <- c(sort(unique(unlist(strsplit(data_in$CDR3,split="")))),"X")
epitopes <- unique(data_in$Epitope)

# Extend the CDR3b sequences
extend_seq <- function(input_seq,max_length){
  # input_seq <- data_in$CDR3[1]
  # max_length <- max(nchar(data_in$CDR3))
  input_vector <- strsplit(input_seq,split="")[[1]]
  nchar_input <- nchar(input_seq)
  
  if(nchar_input < max_length){
    middle_point <- round(nchar_input / 2)
    insert_length <- max_length - nchar_input + 1
    
    pep_start  <- paste0(input_vector[1:middle_point],collapse = "")
    pep_middle <- paste0(rep("X",insert_length),collapse = "")
    pep_end    <- paste0(input_vector[(middle_point+1):nchar_input],collapse = "")
    return(paste0(pep_start,pep_middle,pep_end,collapse = ""))
  }else{
    return(input_seq)
  }
}
data_in$CDR3_extended <- sapply(data_in$CDR3,extend_seq,max_length = max(nchar(data_in$CDR3))+1)

# One-hot the CDR3b
one_hot_peps <- function(string_in){
  vector_in <- strsplit(string_in,split = "")[[1]]
  cdr3_dt <- data.table(cdr3 = factor(vector_in,levels = AAs))
  one_hot_matrix <- as.matrix(one_hot(cdr3_dt))
  one_hot_matrix_l <- melt(one_hot_matrix)
  binary_out <- one_hot_matrix_l$value
  names(binary_out) <- paste0(one_hot_matrix_l$Var2,one_hot_matrix_l$Var1)
  return(binary_out)
}
CDR3_oh <- t(sapply(data_in$CDR3_extended,one_hot_peps))

# One-hot the epitopes
epitope_dt <- data.table(factor(data_in$Epitope,levels = unique(data_in$Epitope)))
epitope_oh <- as.matrix(one_hot(epitope_dt))

# Combine data and one hot
data_with_oh <- data.frame(cbind(CDR3_oh,epitope_oh),stringsAsFactors = FALSE)

# Split training and test data
indices_per_epitope <- lapply(unique(data_in$Epitope),function(this_epitope) which(data_in$Epitope==this_epitope))
training_indices_per_epitope <- lapply(indices_per_epitope,function(this_indices) sample(x = this_indices,
                                                                                         size = length(this_indices)*0.8,
                                                                                         replace = FALSE))
training_indices <- unlist(training_indices_per_epitope)

training_data <- data_with_oh[training_indices,]
testing_data  <- data_with_oh[-training_indices,]

shuffled_successes <- sapply(1:1000,function(x) sum(data_in$Epitope == sample(x = data_in$Epitope,size = length(data_in$Epitope),replace = FALSE))/nrow(data_in))
mean(shuffled_successes)

# Run xgboost model on each and build CMs
getCM <- function(preds,truth,threshold){
  rounded_preds_vs_truth <- data.frame(pred = as.numeric(preds>threshold),
                                       truth = truth)
  overall_CM <- 
    sapply(c(0,1),function(i)
      sapply(c(0,1),function(j)
        length(intersect(which(rounded_preds_vs_truth$pred==i),
                         which(rounded_preds_vs_truth$truth==j)))))
  rownames(overall_CM) <- c("predn=0","predn=1")
  colnames(overall_CM) <- c("truth=0","truth=1")
  return(overall_CM)
}

CDR3_oh_indices <- grep("cdr3",colnames(training_data))
training_prediction <- matrix(0,nrow(training_data),length(epitopes))
testing_prediction  <- matrix(0,nrow(testing_data),length(epitopes))
overall_CM_training <- overall_CM_testing <- bst_model <- vector(mode = "list",length = length(epitopes))
for(i in 1:length(epitopes)){
  this_epitope             <- epitopes[i]
  epitope_col_indx         <- grep(this_epitope,colnames(training_data))
  bst_model[[i]]           <- xgboost(data = as.matrix(training_data[,CDR3_oh_indices]),
                                      label = training_data[,epitope_col_indx],
                                      nrounds = xgboost_rounds,early_stopping_rounds = 50,
                                      params = list(objective = "binary:logistic"))
  training_prediction[,i]  <- predict(bst_model[[i]],as.matrix(training_data[,CDR3_oh_indices]))
  testing_prediction[,i]   <- predict(bst_model[[i]],as.matrix(testing_data[,CDR3_oh_indices]))
  
}

# Plot confusion matrices
acceptance_threshold <- 0.5
overall_CM_training <- lapply(1:length(epitopes),
                              function(i) getCM(preds = training_prediction[,i],
                                                truth = training_data[,grep(epitopes[i],colnames(training_data))],
                                                threshold = acceptance_threshold))
overall_CM_testing <- lapply(1:length(epitopes),
                             function(i) getCM(preds = testing_prediction[,i],
                                               truth = testing_data[,grep(epitopes[i],colnames(testing_data))],
                                               threshold = acceptance_threshold))


names(overall_CM_testing) <- names(overall_CM_training) <- epitopes
# overall_CM_testing_propn <- lapply(overall_CM_testing,function(x) x/sum(x))
# overall_CM_training_propn <- lapply(overall_CM_training,function(x) x/sum(x))

testing_CM_long <- do.call(rbind,lapply(1:length(overall_CM_testing),
                                        function(i) cbind(melt(overall_CM_testing[[i]]),
                                                          target = epitopes[i])))
training_CM_long <- do.call(rbind,lapply(1:length(overall_CM_training),
                                         function(i) cbind(melt(overall_CM_training[[i]]),
                                                           target = epitopes[i])))
testing_cm_plot <- 
  ggplot(testing_CM_long,aes(x = Var1,y=Var2,fill = value,label = value)) + 
  geom_tile() + 
  geom_text() + 
  theme_bw() +
  xlab("") + ylab("") +
  theme(legend.position = "none") +
  facet_grid(.~target,scales = "free") + 
  scale_fill_continuous()


training_cm_plot <- 
  ggplot(training_CM_long,aes(x = Var1,y=Var2,fill = value,label = value)) + 
  geom_tile() + 
  geom_text() + 
  theme_bw() +
  xlab("") + ylab("") +
  theme(legend.position = "none") +
  facet_grid(.~target)
cm_plots_combined <- plot_grid(testing_cm_plot,
                               training_cm_plot,nrow = 2)
cm_plots_combined
ggsave(filename = "plots/cm_heatmaps.png",
       plot = cm_plots_combined,height = 5,width = 7)

# For each output, pick the highest confidence one, return epitope and confidence
make_evaluation_df <- function(prediction_data,input_data){
  prediction_epitope_indx <- apply(prediction_data,1,which.max)
  prediction_epitope      <- epitopes[prediction_epitope_indx]
  prediction_confidence   <- apply(prediction_data,1,function(x) max(x)/sum(x))
  
  true_epitope_indx       <- apply(input_data[,grep("V1_",colnames(input_data))],1,which.max)
  true_epitope            <- epitopes[true_epitope_indx]
  
  data_evaluation         <- data.frame(Epitope = true_epitope,
                                        Predicted_epitope = prediction_epitope,stringsAsFactors = FALSE)
  data_evaluation$Match   <- data_evaluation$Epitope == data_evaluation$Predicted_epitope
  data_evaluation$Confidence <- prediction_confidence
  return(data_evaluation)
}
training_data_evaluation  <- make_evaluation_df(prediction_data = training_prediction,input_data = training_data)
testing_data_evaluation   <- make_evaluation_df(prediction_data = testing_prediction,input_data = testing_data)

training_percent_accuracy <- round(100 * (sum(training_data_evaluation$Match) / nrow(training_data_evaluation)),1)
testing_percent_accuracy  <- round(100 * (sum(testing_data_evaluation$Match) / nrow(testing_data_evaluation)),2)

# Analyse confidence
assess_confidence <- function(data_evaluation_in,plot_window_size = 50){
  
  data_evaluation_sorted <- data_evaluation_in[sort.int(data_evaluation_in$Confidence,index.return = TRUE)$ix,]
  data_evaluation_sorted$Index <- 1:nrow(data_evaluation_sorted)
  waterfall_plot <- 
    ggplot(data_evaluation_sorted,aes(x = Index,
                                      y = Confidence,
                                      fill = Match)) +
    geom_bar(stat = "identity") + theme_bw() +
    theme(legend.position = "none")
  
  
  window_open <- 1:(nrow(data_evaluation_sorted)-plot_window_size+1)
  sliding_inds <- lapply(window_open,function(i) c(i:(i+plot_window_size-1)))
  sliding_values <- sapply(sliding_inds,function(inds) sum(data_evaluation_sorted$Match[inds])/length(inds))
  sliding_df <- data.frame(Index = window_open,Proportion_correct = sliding_values)
  window_plot <- 
    ggplot(sliding_df,aes(x = Index,y = Proportion_correct)) + 
    geom_line() + 
    theme_bw() + 
    theme(legend.position = "none")
  
  histogram_plot <-
    ggplot(data_evaluation_sorted,aes(x = Confidence)) + 
    geom_histogram(bins = 50) + 
      xlim(c(0,1))
  
  return(plot_grid(window_plot,
                   waterfall_plot,
                   histogram_plot,nrow = 3))
}

testing_conf_plot <- assess_confidence(testing_data_evaluation)
training_conf_plot <- assess_confidence(training_data_evaluation)

confidence_plots <- plot_grid(testing_conf_plot,training_conf_plot)
ggsave(plot = confidence_plots,filename = "plots/confidence_plots.png",height = 12,width = 8)

# Make non_binary confusion matrix
make_long_cm <- function(input_data){
  data_evaluation_by_epitope <- split(input_data,f = input_data$Epitope)
  
  # Build CM per epitope
  confusion_matrix <- 
    sapply(data_evaluation_by_epitope,function(x) 
      sapply(epitopes,function(this_epitope) 
        length(which(x$Predicted_epitope==this_epitope))))
  colnames(confusion_matrix) <- paste0("pred_",colnames(confusion_matrix))
  rownames(confusion_matrix) <- paste0("true_",rownames(confusion_matrix))
  
  prop_preds_correct <- apply(confusion_matrix,2,function(x) x/sum(x))
  prop_truth_correct <- t(apply(confusion_matrix,1,function(x) x/sum(x)))
  
  # Make long and fix epitope strings
  confusion_matrix_l <- melt(confusion_matrix,
                             varnames = c("Truth","Prediction"),
                             value.name = "Count")
  confusion_matrix_l$Truth <- gsub(pattern = "true_",
                                   replacement = "",
                                   x = confusion_matrix_l$Truth)
  confusion_matrix_l$Truth <- factor(confusion_matrix_l$Truth,levels = sort(unique(confusion_matrix_l$Truth)))
  
  confusion_matrix_l$Prediction <- gsub(pattern = "pred_",
                                        replacement = "",
                                        x = confusion_matrix_l$Prediction)
  confusion_matrix_l$Prediction <- factor(confusion_matrix_l$Prediction,levels = sort(unique(confusion_matrix_l$Prediction)))
  
  confusion_matrix_l$PropPredsCorrect <- melt(prop_preds_correct)[,3]
  confusion_matrix_l$PropTruthCorrect <- melt(prop_truth_correct)[,3]
  return(confusion_matrix_l)
}
training_cm_long <- make_long_cm(training_data_evaluation)
testing_cm_long  <- make_long_cm(testing_data_evaluation)

training_cm_long_highConfidence <- make_long_cm(training_data_evaluation[which(training_data_evaluation$Confidence>0.8),])
testing_cm_long_highConfidence  <- make_long_cm(testing_data_evaluation[which(testing_data_evaluation$Confidence>0.8),])

sum(testing_cm_long_highConfidence$Count[testing_cm_long_highConfidence$Truth==testing_cm_long_highConfidence$Prediction]) / 
  sum(testing_cm_long_highConfidence$Count)

sum(training_cm_long_highConfidence$Count[training_cm_long_highConfidence$Truth==training_cm_long_highConfidence$Prediction]) / 
  sum(training_cm_long_highConfidence$Count)

allhc <- rbind(testing_cm_long_highConfidence,training_cm_long_highConfidence)
truth <- unlist(sapply(1:nrow(allhc),function(i) rep(as.character(allhc$Truth)[i],allhc$Count[i])))
sum(truth==sample(x = truth,size = length(truth),replace = FALSE))/length(truth)

# Plot confusion matrices
cm_plot <- function(data_in,fill_column,title_text){
  ggplot(data_in,aes_string(x = "Truth",
                            y = "Prediction",
                            fill = fill_column,
                            label = "Count")) + 
    geom_tile() + 
    geom_text() + 
    # scale_fill_viridis(limits=c(0, 1),option = "B") + 
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(title_text) 
}
g1_training <- cm_plot(data_in = training_cm_long,
                       fill_column = "PropPredsCorrect",
                       title_text = "Proportion of predictions that are correct\n Training data")
g2_training <- cm_plot(data_in = training_cm_long,
                       fill_column = "PropTruthCorrect",
                       title_text = "Proportion of true matches that are correct\n Training data")
g1_testing  <- cm_plot(data_in = testing_cm_long,
                       fill_column = "PropPredsCorrect",
                       title_text = "Proportion of predictions that are correct\n Testing data")
g2_testing  <- cm_plot(data_in = testing_cm_long,
                       fill_column = "PropTruthCorrect",
                       title_text = "Proportion of true matches that are correct\n Testing data")
g_grid <- plot_grid(g1_training,g2_training,
                    g1_testing,g2_testing)
ggsave(filename = "plots/non_binary_confusion_matrices.png",
       plot = g_grid,height = 9,width = 8)

g1_training_highConfidence <- cm_plot(data_in = training_cm_long_highConfidence,
                                      fill_column = "PropPredsCorrect",
                                      title_text = "Proportion of HC predictions that are correct\n Training data")
g2_training_highConfidence <- cm_plot(data_in = training_cm_long_highConfidence,
                                      fill_column = "PropTruthCorrect",
                                      title_text = "Proportion of HC true matches that are correct\n Training data")
g1_testing_highConfidence  <- cm_plot(data_in = testing_cm_long_highConfidence,
                                      fill_column = "PropPredsCorrect",
                                      title_text = "Proportion of HC predictions that are correct\n Testing data")
g2_testing_highConfidence  <- cm_plot(data_in = testing_cm_long_highConfidence,
                                      fill_column = "PropTruthCorrect",
                                      title_text = "Proportion of HC true matches that are correct\n Testing data")
g_grid_highConfidence <- plot_grid(g1_training_highConfidence,g2_training_highConfidence,
                                   g1_testing_highConfidence,g2_testing_highConfidence)
ggsave(filename = "plots/non_binary_confusion_matrices_highConfidence.png",
       plot = g_grid_highConfidence,height = 9,width = 8)


# Pull out top features that are predictive, plot predictivity along the peptide


# Do tSNE clustering of CDR3, colour by epitope
peptide_distance <- function(pepA,pepB){
  pepA_vector <- unlist(strsplit(pepA,split=""),use.names = FALSE)
  pepB_vector <- unlist(strsplit(pepB,split=""),use.names = FALSE)
}
CDR3_vector_list <- strsplit(data_in$CDR3_extended,split="")
CDR3_vector_list_inds <- lapply(CDR3_vector_list,function(x) match(x,colnames(BLOSUM62)))
comb_ind <- t(combn(1:length(CDR3_vector_list_inds),m = 2))
# blosum_scores <- 
#   sapply(1:nrow(comb_ind),function(i) {
#     ind1 <- CDR3_vector_list_inds[[comb_ind[i,1]]]
#     ind2 <- CDR3_vector_list_inds[[comb_ind[i,2]]]
#     mean(BLOSUM62[((ind2-1) * nrow(BLOSUM62)) + ind1])
#   })
blosum_scores <- 
  sapply(1:nrow(comb_ind),function(i) {
    ind1 <- CDR3_vector_list_inds[[comb_ind[i,1]]][5:15]
    ind2 <- CDR3_vector_list_inds[[comb_ind[i,2]]][5:15]
    mean(BLOSUM62[((ind2-1) * nrow(BLOSUM62)) + ind1])
  })
blosum_matrix <- matrix(0,length(CDR3_vector_list),length(CDR3_vector_list))
blosum_matrix[comb_ind] <- blosum_scores
blosum_matrix[comb_ind[,c(2,1)]] <- blosum_scores
duplicated_rows <- which(duplicated(blosum_matrix))
colnames(blosum_matrix) <- data_in$CDR3_extended
rownames(blosum_matrix) <- data_in$CDR3_extended
if(length(duplicated_rows)>0){
  blosum_matrix <- blosum_matrix[-duplicated_rows,-duplicated_rows]
}

sample_indx <- sample(1:nrow(blosum_matrix),
                      size = 1000,
                      replace = FALSE)
blosum_matrix_samples <- blosum_matrix[sample_indx,sample_indx]


tsne_output <- Rtsne(blosum_matrix_samples,max_iter = 1000,num_threads = 3)
tsne_output_coords <- data.frame(tsne_output$Y)
colnames(tsne_output_coords) <- c("x","y")
tsne_output_coords$cdr3 <- colnames(blosum_matrix_samples)
tsne_output_coords$target <- data_in$Epitope[match(tsne_output_coords$cdr3,data_in$CDR3_extended)]

hbdscan_out <- hdbscan(data.frame(x = tsne_output_coords$x,
                                  y = tsne_output_coords$y),minPts = 50)
tsne_output_coords$cluster <- hbdscan_out$cluster + 1


g_points <- 
  ggplot(tsne_output_coords,aes(x = x,y = y)) + 
  geom_point(size = 3,fill = "grey20",alpha = 0.3) + 
  theme_bw() + 
  theme(legend.position = "none")

g_numbers <- 
  ggplot(tsne_output_coords,aes(x = x,y = y,label = cluster)) + 
  geom_text(aes(col = as.factor(cluster))) +
  theme_bw() + 
  theme(legend.position = "none")
g_grid <- plot_grid(g_points, g_numbers)
ggsave(filename = "plots/cdr3_tSNE_out.png",plot = g_grid,width = 14,height = 7)
# geom_point(size = 2,aes(fill = as.factor(cluster)))
# facet_wrap(.~target)

g_byEpitope <- ggplot(tsne_output_coords,aes(x = x,y = y,label = cdr3)) + 
  geom_point(size = 2,shape = 21,aes(fill = target),alpha = 0.5) +
  facet_wrap(.~target)
ggsave(filename = "plots/cdr3_tSNE_byepitope.png",plot = g_byEpitope,width = 10,height = 8)

cdr3_by_cluster <- split(tsne_output_coords,
                         f = tsne_output_coords$cluster)

g_seqlogo <- ggseqlogo(lapply(cdr3_by_cluster,"[[","cdr3"),facet = "grid")
ggsave(filename = "plots/cdr3_cluster_seqlogos.png",plot = g_seqlogo,width = 30,height = 4)


save.image("2020-09-07b.RData")
# Do RoC analysis
# library(ROCR)
# 
# # Use ROCR package to plot ROC Curve
# xgb.pred <- prediction(pred, test.label)
# xgb.perf <- performance(xgb.pred, "tpr", "fpr")
# 
# plot(xgb.perf,
#      avg="threshold",
#      colorize=TRUE,
#      lwd=1,
#      main="ROC Curve w/ Thresholds",
#      print.cutoffs.at=seq(0, 1, by=0.05),
#      text.adj=c(-0.5, 0.5),
#      text.cex=0.5)
# grid(col="lightgray")
# axis(1, at=seq(0, 1, by=0.1))
# axis(2, at=seq(0, 1, by=0.1))
# abline(v=c(0.1, 0.3, 0.5, 0.7, 0.9), col="lightgray", lty="dotted")
# abline(h=c(0.1, 0.3, 0.5, 0.7, 0.9), col="lightgray", lty="dotted")
# lines(x=c(0, 1), y=c(0, 1), col="black", lty="dotted")

