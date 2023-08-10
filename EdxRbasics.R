install.packages("dslabs")
library(dslabs)
data(heights)
class(heights)
class(heights$sex)
class(heights$height)
class("Male")
class(75.00000)
heights$sex[777]
heights[777,1]
max_height <- max(heights$height)
max_height
min_height <- min(heights$height)
min_height
min_height_row <- which.min(heights$height)
min_height_row
mean_height <- mean(heights$height)
mean_height
median_height <- median(heights$height)
median_height
proportion_male <- sum(heights$sex == "Male") / nrow(heights)
proportion_male
sum(heights$height > 78)
sum(heights$sex == "Female" & heights$height > 78)
