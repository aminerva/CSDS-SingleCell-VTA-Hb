## Functions for 1D SDA (training on two groups only)

run_1D_sda_projection <- function(data_t, labels, shuffle = FALSE, genes = NULL, 
                                           train_conditions = c("Control", "Resilient"), 
                                           project_condition = NULL, 
                                           test_data = NULL, test_labels = NULL) {
    
    # Shuffle the labels if requested
    if (shuffle == TRUE) {
        message("Shuffling training labels before creating folds")
        set.seed(123)
        labels <- sample(labels)
    } 
   
    # If genes is NULL, use all genes (otherwise will taken variable genes as input)
    if (is.null(genes)) {
        genes <- colnames(data_t) 
    }

    message("Subsetting for training conditions")
    train_idx <- rownames(data_t)[labels %in% train_conditions]
    X_train <- as.matrix(data_t[train_idx, genes, drop = FALSE])
    
    y_train <- as.factor(labels[train_idx])
  
    message("Removing genes with constant expression in at least one group")
    constant_genes <- colnames(X_train)[apply(X_train, 2, function(gene_values) {
        tapply(gene_values, y_train, function(x) length(unique(x)) == 1) %>% any()})]
    X_train <- X_train[, !colnames(X_train) %in% constant_genes, drop = FALSE]
  
    message("Training SDA model")
    sda_model <- sda(X_train, y_train)

    message("Extracting vector of LD1 projections for the training conditions (", 
        train_conditions[1], " and ", train_conditions[2], ")")
    sda_train_scores <- as.matrix(X_train) %*% sda_model$beta[1, ] 
    
    message("Evaluating training accuracy")
    preds_train <- predict(sda_model, X_train)$class # Extract predicted class labels
    acc_train <- mean(preds_train == y_train) # Determine overall accuracy
    
    sda_train_df <- data.frame(
        LD1 = sda_train_scores, 
        TrueLabel = y_train, 
        PredictedLabel = preds_train,
        Source = "Train", 
        Cell = rownames(X_train))

    # Project data if project_condition = TRUE
    sda_project_df <- data.frame(
        LD1 = numeric(0), 
        TrueLabel = character(0), 
        PredictedLabel = character(0), 
        Source = character(0), 
        Cell = character(0))
    if (!is.null(project_condition)) {
        message("Projecting ", project_condition, " onto LDA space (using the original labels for this condition)")
        project_idx <- rownames(data_t)[labels == project_condition]
        X_project <- data_t[project_idx, genes, drop = FALSE] # Grab only projected condition and genes of interest
        X_project <- as.matrix(X_project[, colnames(X_train), drop = FALSE]) # Ensure matching columns with training data
        sda_project_scores <- X_project %*% sda_model$beta[1, ]  # Project the new data onto LD1
        sda_project_df <- data.frame(
            LD1 = sda_project_scores, 
            TrueLabel = project_condition, 
            PredictedLabel = NA,
            Source = "Project",
            Cell = rownames(X_project))
    }

    # Test held out data 
    sda_test_df <- data.frame(
        LD1 = numeric(0), 
        TrueLabel = character(0), 
        PredictedLabel = character(0), 
        Source = character(0), 
        Cell = character(0))
    if (!is.null(test_data) && !is.null(test_labels)) {
        message("Projecting held out test data onto LDA space")
        X_test <- as.matrix(test_data[, colnames(X_train), drop = FALSE])
        sda_test_scores <- X_test %*% sda_model$beta[1, ]
        preds_test <- predict(sda_model, X_test)$class
        acc_test <- mean(preds_test == test_labels)
        sda_test_df <- data.frame(
            LD1 = sda_test_scores, 
            TrueLabel = test_labels, 
            PredictedLabel = preds_test,
            Source = "Test",
            Cell = rownames(X_test))
    } else {
        acc_test <- NA
    }

    sda_combined_df <- rbind(sda_train_df, sda_test_df, sda_project_df)
    
    message("All done!")

    return(list(
        sda_model = sda_model,
        sda_df = sda_combined_df, 
        train_accuracy = acc_train, 
        test_accuracy = acc_test))
}


cross_validate_1D_sda <- function(data_t, labels, genes, sample_ids, shuffle_bool=FALSE, 
                                train_conditions = c("Control", "Resilient")) { 

    # Only keep data from training conditions
    is_train <- labels %in% train_conditions
    data_train <- data_t[is_train, ]
    labels_train <- labels[is_train]
    sample_ids <- sample_ids[is_train]

    if (shuffle_bool) {
        message("Shuffling training labels before creating folds")
        set.seed(123)
        labels_train <- sample(labels_train)
    }

    # Leave-one-sample-out: each sample ID gets held out once
    sample_info <- data.frame(sample = unique(sample_ids), stringsAsFactors = FALSE)
    sample_info$group <- sapply(sample_info$sample, function(s) { unique(labels[sample_ids == s]) })
    folds <- as.list(sample_info$sample)
    
    n_folds <- length(folds) 
    fold_accuracies <- numeric(n_folds)
    fold_projections <- vector("list", n_folds)
    held_out_sample_list <- vector("list", n_folds)
    fold_lambdas <- numeric(n_folds)

    message("Running ", n_folds, "-fold cross validation with ", 
             ifelse(shuffle_bool, "shuffled", "unshuffled"), " labels")

    for (i in seq_along(folds)) {
        message("Starting fold ", i)

        # Get test and train indices
        held_out <- folds[[i]]
        held_out_sample_list[[i]] <- held_out

        test_idx <- sample_ids %in% held_out
        train_idx <- !test_idx

        # Subset data
        data_fold_train <- data_train[train_idx, ]
        data_fold_test <- data_train[test_idx, ]
        labels_fold_train <- labels_train[train_idx]
        labels_fold_test <- labels_train[test_idx]

        message("Fitting SDA model")
        sda_result <- run_1D_sda_projection(
            data_t = rbind(data_fold_train, data_fold_test),  # full data, but model only sees training portion
            labels = c(labels_fold_train, labels_fold_test),  # full labels
            genes = genes,
            shuffle = FALSE, #### Shuffle has already been handled above
            train_conditions = train_conditions,
            project_condition = NULL, 
            test_data = data_fold_test,
            test_labels = labels_fold_test)

        sda_result$sda_df$Fold <- i
        fold_projections[[i]] <- subset(sda_result$sda_df, Source == "Test")
        fold_accuracies[i] <- sda_result$test_accuracy
        fold_lambdas[i]<- sda_result$sda_model$regularization[['lambda']]
        
        message("Lambda used in fold ", i, ": ", sda_result$sda_model$regularization[['lambda']])

        message("Done with fold ", i)
        message("")
    }

    combined_df <- do.call(rbind, fold_projections)

    message("All done!")
    
    return(list(
        accuracies = fold_accuracies,
        lambdas = fold_lambdas, 
        sda_projections = combined_df,
        held_out_sample = held_out_sample_list))
}


## Functions for 2D SDA (training on all three groups)

run_2D_sda_projection <- function(data_t, labels, shuffle = FALSE, genes = NULL, 
                                 train_conditions = c("Control", "Susceptible", "Resilient"), 
                                 test_data = NULL, test_labels = NULL) {

  if (shuffle) {
    labels <- sample(labels)
  }

  if (is.null(genes)) {
    genes <- colnames(data_t)
  }

  train_idx <- rownames(data_t)[labels %in% train_conditions]
  X_train <- as.matrix(data_t[train_idx, genes, drop = FALSE])
  y_train <- as.factor(labels[train_idx])

  constant_genes <- colnames(X_train)[apply(X_train, 2, function(g) {
    tapply(g, y_train, function(x) length(unique(x)) == 1) %>% any()
  })]
  X_train <- X_train[, !colnames(X_train) %in% constant_genes, drop = FALSE]

  sda_model <- sda(X_train, y_train)

  proj_train <- as.matrix(X_train) %*% t(sda_model$beta[1:2, , drop = FALSE])
  colnames(proj_train) <- c("LD1", "LD2")
  
  preds_train <- predict(sda_model, X_train)$class

  df_train <- data.frame(
    LD1 = proj_train[, 1], 
    LD2 = proj_train[, 2], 
    TrueLabel = y_train, 
    PredictedLabel = preds_train,
    Source = "Train")
  df_train$Cell <- rownames(X_train)

  acc_train <- mean(predict(sda_model, X_train)$class == y_train)

  df_test <- data.frame()
  acc_test <- NA
  if (!is.null(test_data) && !is.null(test_labels)) {
    X_test <- as.matrix(test_data[, colnames(X_train), drop = FALSE])
    proj_test <- X_test %*% t(sda_model$beta[1:2, , drop = FALSE])
    colnames(proj_test) <- c("LD1", "LD2")

    preds_test <- predict(sda_model, X_test)$class # Get predicted labels

    df_test <- data.frame(
        LD1 = proj_test[, 1],
        LD2 = proj_test[, 2],
        TrueLabel = test_labels, 
        PredictedLabel = preds_test, 
        Source = "Test")
    df_test$Cell <- rownames(X_test)

    acc_test <- mean(predict(sda_model, X_test)$class == test_labels)
  }

  df_combined <- rbind(df_train, df_test)

  return(list(
    sda_model = sda_model,
    sda_df = df_combined,
    train_accuracy = acc_train,
    test_accuracy = acc_test
  ))
}

cross_validate_2D_sda <- function(data_t, labels, genes, sample_ids, shuffle_bool = FALSE, 
                                  train_conditions = c("Control", "Susceptible", "Resilient")) {

    is_train <- labels %in% train_conditions
    data_train <- data_t[is_train, ]
    labels_train <- labels[is_train]
    sample_ids <- sample_ids[is_train]

    if (shuffle_bool) {
        set.seed(123)
        labels_train <- sample(labels_train)
    }

    # Generate folds for cv by holding out one sample until all samples are held out once
    sample_info <- data.frame(sample = unique(sample_ids), stringsAsFactors = FALSE)
    sample_info$group <- sapply(sample_info$sample, function(s) {unique(labels[sample_ids == s]) })
    folds <- as.list(sample_info$sample) 

    n_folds <- length(folds)
    fold_accuracies <- numeric(n_folds)
    fold_projections <- vector("list", n_folds)
    fold_lambdas <- numeric(n_folds)
    held_out_sample_list <- vector("list", n_folds)

    for (i in seq_along(folds)) {
        message("Starting fold ", i)

        held_out_sample <- folds[[i]]
        held_out_sample_list[[i]] <- held_out_sample
        message("Holding out ", held_out_sample, " for training")

        test_idx <- sample_ids %in% held_out_sample
        train_idx <- !test_idx

        data_fold_train <- data_train[train_idx, ]
        data_fold_test <- data_train[test_idx, ]
        labels_fold_train <- labels_train[train_idx]
        labels_fold_test <- labels_train[test_idx]

        sda_result <- run_2D_sda_projection(
            data_t = rbind(data_fold_train, data_fold_test),
            labels = c(labels_fold_train, labels_fold_test),
            genes = genes,
            shuffle = FALSE,
            train_conditions = train_conditions,
            test_data = data_fold_test,
            test_labels = labels_fold_test)

        sda_result$sda_df$Fold <- i
        fold_projections[[i]] <- subset(sda_result$sda_df, Source == "Test")

        fold_accuracies[i] <- sda_result$test_accuracy
        fold_lambdas[i] <- sda_result$sda_model$regularization[['lambda']]
        message("Lambda used in fold ", i, ": ", sda_result$sda_model$regularization[['lambda']])

        message("Done with fold ", i)
        message("")
    }

    combined_df <- do.call(rbind, fold_projections)

    message("All done!")

    return(list(
        accuracies = fold_accuracies,
        lambdas = fold_lambdas, 
        sda_projections = combined_df, 
        held_out_samples = held_out_sample_list))
}

