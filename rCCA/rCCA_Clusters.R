---
title: "CCA_Clusters_R"
author: "Moemen Hussein"
output: html_document
---

## 1. Package setup

```{r package-setup, echo=FALSE}

library(mixOmics)
library(vegan)
library(readxl)
library(permute)
library(lattice)
library(MASS)
library(ggplot2)
library(ggExtra)
library(stats)
library(base)
library(dplyr)
library(cowplot)
library(igraph)
library(openxlsx)
library(plotly)
library(readxl)
library(writexl)

set.seed(3189)

```

## 2. Load metadata

```{r load-metadata, echo=FALSE}
## 

# Load the data
meta_cleaned_imputed <- read.csv("HERE") |> as.data.frame() 
```

## 3. Load cluster data

```{r load-align-cluster-data, echo=FALSE}

clusters_data <- read_excel(path = "HERE", sheet = 4)

```

## 4. Fit ridge regularized CCA

```{r fit-rcca, echo=FALSE}

# Perform rCCA using ridge regularization
Y <-  meta_cleaned_imputed
X_all <- clusters_data  

ridge_rcc_all_meta <- mixOmics::rcc(X_all, Y, ncomp = ncol(clusters_data), method = "ridge", lambda1 = 0.35, lambda2 = 0.2) 
# ridge_rcc_all_meta <- mixOmics::rcc(X_all, Y, ncomp = ncol(X_all) ) 

# Display the optimal lambda values after shrinkage
print(ridge_rcc_all_meta$lambda)

# Display the correlations between the X & Y components
print(ridge_rcc_all_meta$cor)

# Bar plot of canonical correlations
plot(ridge_rcc_all_meta, type = "barplot", main = "Ridge")
```

## 5. Plot variables across canonical components

```{r plot-canonical-variables, echo=FALSE}

# Variable plots
plotVar(ridge_rcc_all_meta, var.names = c(TRUE, FALSE), cex = c(2, 2), cutoff = 0, comp = c(1,2), title = 'rCCA Ridge comp 1 - 2')

```

## 6. Extract and export canonical variates

```{r export-variates, echo=FALSE}



# Objects:
variates_X <- ridge_rcc_all_meta$variates$X
variates_Y <- ridge_rcc_all_meta$variates$Y
variates_XY <- ridge_rcc_all_meta$variates$X + ridge_rcc_all_meta$variates$Y
variates_XY_df <- as.data.frame(variates_XY)

# Export each to Excel

```

## 7. Calculate and export loadings

```{r export-loadings, echo=FALSE}

loadings_meta <- cor(variates_XY_df,X_all)
loadings_meta_df <- as.data.frame(loadings_meta)


loadings_mic <- cor(variates_XY_df,Y)
loadings_mic_df <- as.data.frame(loadings_mic)


loadings_X <- cor(as.data.frame(variates_X), X_all)
loadings_Y <- cor(as.data.frame(variates_Y), Y)

```

## 8. Generate clustered image maps

```{r clustered-image-maps, echo=FALSE}


# Generate the CIM plot
cim_result <- cim(ridge_rcc_all_meta, comp=1:i,col.cex = 0.2, row.cex = 0.5, center = FALSE, clust.method = c("ward", "ward"))
  
# Convert CIM matrix to data frame

# Save to Excel


```


## 9. Grid search

```{r grid-search, echo=FALSE}

# Perform grid search for optimal sCCA penalties

Y <- meta_cleaned_imputed
X_all <- clusters_data 

X <- X_all  # If you want to switch to X_selected, change this line accordingly

lambda_x_values <- seq(0.0001, 1, length.out = 15)
lambda_y_values <- seq(0.0001, 1, length.out = 15)
k <- 5

# For reproducibility
set.seed(123)
folds <- sample(1:k, nrow(X), replace = TRUE)

# Perform baseline CCA
baseline_cca <- mixOmics::rcc(X, Y, ncomp = 1, method = 'ridge', lambda1 = 0.0001, lambda2 = 0.0001)
baseline_score <- baseline_cca$cor[1]

# Initialize data frame to store CV results
CV_results <- data.frame(lambda_x = numeric(), lambda_y = numeric(), CV_scores = numeric())

for (lambda_x in lambda_x_values) {
   for (lambda_y in lambda_y_values) {
     CV_scores <- numeric(k)
     for (fold in 1:k) {
       train_index <- which(folds != fold)
      test_index <- which(folds == fold)

       X_train <- X[train_index, ]
       Y_train <- Y[train_index, ]
       X_test <- X[test_index, ]
       Y_test <- Y[test_index, ]

       # Fit the model on the training data
       model <- mixOmics::rcc(X_train, Y_train, ncomp = 1, method = 'ridge', lambda1 = lambda_x, lambda2 = lambda_y)

       # Extract the canonical variates
       X_variate_train <- model$variates$X
       Y_variate_train <- model$variates$Y

       # Predict the test set variates using the model coefficients
       X_variate_test <- as.matrix(X_test) %*% model$loadings$X
       Y_variate_test <- as.matrix(Y_test) %*% model$loadings$Y

       # Compute the correlation between the canonical variates for the test set
       CV_scores[fold] <- cor(X_variate_test[, 1], Y_variate_test[, 1])
     }
     CV_mean_score <- mean(CV_scores, na.rm = TRUE)
     CV_results <- rbind(CV_results, data.frame(lambda_x = lambda_x, lambda_y = lambda_y, CV_scores = CV_mean_score))
   }
 }

 # The results
 print(CV_results)
 
 # Plot with text labels
 ggplot(CV_results, aes(x = lambda_x, y = lambda_y, fill = CV_scores)) +
   geom_tile() +
   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = baseline_score) +
   geom_text(aes(label = round(CV_scores, 2)), color = "black", size = 3) +
   labs(title = "Grid Search Results for sCCA Penalties", x = "Lambda X", y = "Lambda Y", fill = "Score") +
   theme_minimal()

 # Summarize the data to calculate the mean score for each lambda_y
 mean_scores_y <- CV_results %>%
   group_by(lambda_y) %>%
   summarize(CV_mean_score_y = mean(CV_scores))


 # Create the plot with ggplot and add a smooth curve
 ggplot(mean_scores_y, aes(x = lambda_y, y = CV_mean_score_y)) +
   geom_point() +
   geom_smooth(method = "loess", span = 0.5) +
   labs(title = "Mean Score Across Lambda X for Each Lambda Y", x = "Lambda Y", y = "Mean Score") +
   theme_minimal()
 
```

## 10. Calculate and plot mean scores for Lambda X

```{r lambda-x-mean-scores, echo=FALSE}
# Summarize the data to calculate the mean score for each lambda_x
mean_scores_x <- CV_results %>%
  group_by(lambda_x) %>%
  summarize(CV_mean_score = mean(CV_scores))

# Create the plot with ggplot and add a smooth curve
ggplot(mean_scores_x, aes(x = lambda_x, y = CV_mean_score)) +
  geom_point() +
  geom_smooth(method = "loess", span = 0.5) +
  labs(title = "Mean Score Across Lambda Y for Each Lambda X", x = "Lambda X", y = "Mean Score") +
  theme_minimal()
```

## 11. visualization of grid scores

```{r optional-grid-visualizations, echo=FALSE}

# Summarize mean scores for lambda_x and lambda_y
 mean_scores <- CV_results %>%
   group_by(lambda_x, lambda_y) %>%
   summarize(CV_mean_score = mean(CV_scores), .groups = 'drop')
  


# Smoothed curves for lambda_x and lambda_y
 mean_scores_x <- CV_results %>%
   group_by(lambda_x) %>%
   summarize(CV_mean_score = mean(CV_scores))

 mean_scores_y <- CV_results %>%
   group_by(lambda_y) %>%
   summarize(CV_mean_score_y = mean(CV_scores))


# 

 loess_x <- loess(CV_mean_score ~ lambda_x, data = mean_scores_x, span = 0.5)
 smoothed_x <- data.frame(lambda_x = mean_scores_x$lambda_x, score = predict(loess_x, mean_scores_x$lambda_x))

 loess_y <- loess(CV_mean_score_y ~ lambda_y, data = mean_scores_y, span = 0.5)
 smoothed_y <- data.frame(lambda_y = mean_scores_y$lambda_y, score = predict(loess_y, mean_scores_y$lambda_y))



 
 # Plot
 ggplot() +
   geom_line(data = smoothed_x, aes(x = lambda_x, y = score), color = "blue") +
   geom_line(data = smoothed_y, aes(x = lambda_y, y = score), color = "red") +
   labs(title= "Grid Scores" , x = "Lambda", y = "Mean Score") +
   theme_minimal()


```
