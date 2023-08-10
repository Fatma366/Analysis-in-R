library(dslabs)
mnist <- read_mnist()
ncol(mnist$train$images)

library(dslabs)
library(dplyr)
library(lubridate)
data(reported_heights)
data("reported_heights")
dat <- mutate(reported_heights, date_time = ymd_hms(time_stamp)) |>
  filter(date_time >= make_date(2016, 01, 25) & 
           date_time < make_date(2016, 02, 1)) |>
  mutate(type = ifelse(day(date_time) == 25 & hour(date_time) == 8 & 
                         between(minute(date_time), 15, 30),
                       "inclass", "online")) |> select(sex, type)
x <- dat$type
y <- factor(dat$sex, c("Female", "Male"))

# Count the number of females in the "inclass" and "online" group
dat %>% group_by(type) %>% summarize(prop_female = mean(sex == "Female"))

# To predict sex based on the "type" variable and report the accuracy of your prediction
y_hat <- ifelse(x == "online", "Male", "Female") %>% 
  factor(levels = levels(y))
mean(y_hat==y)
table(y_hat, y)
library(caret)
sensitivity(y_hat, y)
specificity(y_hat, y)
prevalence <- (table(dat$sex)["Female"] / nrow(dat)) * 100
prevalence_formatted <- sprintf("%.2f%%", prevalence)
prevalence_formatted
library(caret)
data(iris)
iris <- iris[-which(iris$Species=='setosa'),]
y <- iris$Species
set.seed(76)
test_index <- createDataPartition(y, times=1, p=0.5,list=FALSE)
test <- iris[test_index,]
train <- iris[-test_index,]

library(caret)
data(iris)
iris <- iris[-which(iris$Species == 'setosa'), ]
y <- iris$Species

# Create the train and test partitions
set.seed(76)
test_index <- createDataPartition(y, times = 1, p = 0.5, list = FALSE)
test <- iris[test_index, ]
train <- iris[-test_index, ]
foo <- function(x){
  rangedValues <- seq(range(x)[1], range(x)[2], by=0.1)
  sapply(rangedValues, function(i){
    y_hat <- ifelse(x>i, 'virginica', 'versicolor')
    mean(y_hat==train$Species)
  })
}
predictions <- apply(train[,-5], 2, foo)
sapply(predictions, max)

predictions <- foo(train[,4])
rangedValues <- seq(range(train[,4])[1], range(train[,4])[2], by=0.1)
cutoffs <-rangedValues[which(predictions==max(predictions))]

y_hat <- ifelse(test[,4]>cutoffs[1], 'virginica', 'versicolor')
mean(y_hat==test$Species)

# Assuming you have already defined the train and test datasets

# Find optimal cutoffs for Petal.Length and Petal.Width in the train dataset
length_seq <- seq(min(train$Petal.Length), max(train$Petal.Length), by = 0.1)
width_seq <- seq(min(train$Petal.Width), max(train$Petal.Width), by = 0.1)

best_length_cutoff <- 0
best_width_cutoff <- 0
best_accuracy <- 0

for (length_cutoff in length_seq) {
  for (width_cutoff in width_seq) {
    # Create a rule to predict "virginica" if Petal.Length is greater than the length cutoff AND Petal.Width is greater than the width cutoff
    predictions <- ifelse(train$Petal.Length > length_cutoff & train$Petal.Width > width_cutoff, "virginica", "versicolor")
    
    # Calculate accuracy for this combination of cutoffs
    accuracy <- sum(predictions == train$Species) / length(predictions)
    
    # Update best cutoffs and accuracy if the current combination has higher accuracy
    if (accuracy > best_accuracy) {
      best_length_cutoff <- length_cutoff
      best_width_cutoff <- width_cutoff
      best_accuracy <- accuracy
    }
  }
}

# Calculate overall accuracy for the test data using the best cutoffs
test_predictions <- ifelse(test$Petal.Length > best_length_cutoff & test$Petal.Width > best_width_cutoff, "virginica", "versicolor")

accuracy <- sum(test_predictions == test$Species) / length(test_predictions)
cat("Overall Accuracy for Test Data:", accuracy * 100, "%\n")

set.seed(1) 
disease <- sample(c(0,1), size=1e6, replace=TRUE, prob=c(0.98,0.02))
test <- rep(NA, 1e6)
test[disease==0] <- sample(c(0,1), size=sum(disease==0), replace=TRUE, prob=c(0.90,0.10))
test[disease==1] <- sample(c(0,1), size=sum(disease==1), replace=TRUE, prob=c(0.15, 0.85))
# Calculate the probability of a positive test
mean(test == 1)

# Calculate the probability of testing negative
negative_probability <- mean(test == 0)

# Calculate the probability of having the disease if the test is negative
disease_if_negative <- sum(disease == 1 & test == 0) / sum(test == 0)
disease_if_negative
The value of disease_if_negative should give you the correct probability that an individual has the disease if the test is negative. If you are still getting a different result, please double-check the code and ensure you are running it correctly.




User
The command gives me 0.003461 as the answer and this is not the correct answer.please modify the code to give me right answer
ChatGPT
Apologies for the oversight. Let's calculate the probability again using Bayes' theorem:
  
  R
Copy code
# Sample code provided
set.seed(1) 
disease <- sample(c(0, 1), size = 1e6, replace = TRUE, prob = c(0.98, 0.02))
test <- rep(NA, 1e6)
test[disease == 0] <- sample(c(0, 1), size = sum(disease == 0), replace = TRUE, prob = c(0.90, 0.10))
test[disease == 1] <- sample(c(0, 1), size = sum(disease == 1), replace = TRUE, prob = c(0.15, 0.85))

# Calculate the probability of testing positive
positive_probability <- mean(test == 1)

# Calculate the probability of testing negative
negative_probability <- mean(test == 0)

# Calculate the probability of having the disease if the test is negative
disease_if_negative <- sum(disease == 1 & test == 0) / sum(test == 0)
disease_if_negative

mean(disease[test==1]==1)

mean(disease[test==1]==1)/mean(disease==1)

library(dslabs)
data("heights")
heights %>% 
  mutate(height_rounded = round(height)) %>%
  group_by(height_rounded) %>%
  summarize(conditional_prob = mean(sex == "Male")) %>%
  qplot(height_rounded, conditional_prob, data = ., geom = "line")

ps <- seq(0, 1, 0.1)
heights %>% 
  mutate(g = cut(height, quantile(height, ps), include.lowest = TRUE)) %>%
  group_by(g) %>%
  summarize(p = mean(sex == "Male"), height = mean(height)) %>%
  qplot(height, p, data =.)

install.packages("tidyverse")
library(tidyverse)
library(caret)

set.seed(1)
n <- 100
Sigma <- 9*matrix(c(1.0, 0.5, 0.5, 1.0), 2, 2)
dat <- MASS::mvrnorm(n = 100, c(69, 69), Sigma) %>%
  data.frame() %>% setNames(c("x", "y"))
