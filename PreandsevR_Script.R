# Convert excel data into csv using "save as" option in excel file

# Load packages
library(tidyverse)

install.packages("conflicted")
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")

# load data
# read malnutrition data as mal_data
# used the "file.choose" option to allow the code to be useful to anyone not familiar with paths
fat_data <- read.excel(file.choose(), header = TRUE)
fat_data <- Pre_Seve
# if data has blank spaces
# mal_data <- read.csv(file.choose(), na.strings = "")

# view data
View(fat_data)

# check structure of data
str(fat_data)

# check number of rows and columns
glimpse(fat_data)

# check column names
names(fat_data)

colnames(fat_data)

## CORRELATION MATRIX


# Correlation is a statistical measure that indicates the extent to which two or more variables
# fluctuate in relation to each other. A correlation coefficient is a number between -1 and 1
# that expresses the strength and direction of the relationship between two variables.


# Correlation is a useful tool for understanding the relationship between two variables. 
# However, it is important to remember that CORRELATION DOESN'T IMPLY CAUSATION.

# In a correlation plot, 

# x = Explanatory/ independent variable
# Y = Response/ dependent variable

# Positive correlation = X increases, Y increases
# Negative correlation = X increases, Y decreases

# Most commonly correlation is the Pearson product-moment correlation

# N/B: correlation coefficient only measures the strength of linear relationships only.
# Doesn't apply for skewed data.
# What you do is log transformation of the data
# correlation
# Create a correlation matrix between all relevant numeric variables/columns in the dataset.
# create cor dataset
cor_data <- fat_data %>% 
  select(Severity, Prevalence)

# correlation matrix
cor_data_matrix <- cor(cor_data)


view(cor_data_matrix)

## OR

# calculate correlation between two variables
cor(fat_data$Prevalence, fat_data$Severity)

# Install and load ggcorrplot
install.packages("ggcorrplot")
library("ggcorrplot")

# Apply ggcorrplot function
ggcorrplot(cor_data_matrix)
ggsave(ggcorrplot(cor_data_matrix), filename = "Correlation_pre_sev.pdf",
       width = 10, height = 10, dpi = 300)
## BARPLOT

#  melt your data first over Counties.
# It will create another variable called value by default, 
# so you will need to rename it (I called it percent).
# In other words, make prevalence and severity into one column then all the percentages in one column 
library(reshape2)
barplot1 <- melt(fat_data)
names(barplot1)[3] <- "percent"

PreSevScores <- ggplot(barplot1, aes(x = Counties, y = percent, fill = variable)) +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.8)) +
  labs(y = "Disease Prevalence and Severity (%)",
       x = "Counties",
       title = "Prevalence and Severity scores for LY symptomatic plants",
       fill = "Variables (%)") +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 90, by = 10)) +  # Modify the y-axis scale
  theme(text = element_text(size = 18))  # increase font size of labels


ggsave(PreSevScores, filename = "Pre&SevScores.pdf",
       width = 10, height = 10, dpi = 300)

ggsave(PreSevScores, filename = "Pre&SevScores.png",
       width = 10, height = 10, dpi = 300)









