# R package dplyr for data manipulation
# there is a function calle tibble is same as data.frame but it will give the 
# extra fucntions to work with datas
# to access the file inside computer use file.choose()
scRNA_metadata <- file.choose()
library(dplyr)

my_data <- as_tibble(mtcars)

#this is how we can save a file as csv file we can change format like "rds"
#file is a name which we can name it for our file to save

write.csv(my_data, file = "cars.csv")
#filter() used to choose a row by a condition

filtered <- filter(my_data, mpg > 10, .preserve = FALSE)

#there is one subcondition called .preserver default valus is FALSE
#if it is false the empty cell we have after the filteration it will drop 
# form the list if it is TRUE the empty cell will come as it is in the final
#output

data <- filter(my_data, cyl >= 6) %>% select( hp)
data

# in this code we used filter function to filter the data with specific
# condition and from the filtered dataset we selected particular column 
#to display
# we can write the same code like following without pipe

data2 <- select(filter(my_data, cyl >= 6), hp)
data2
# we can also include it with the basic operators and booleans too

#range between some rows only specified inside mentioned row only we can see 
data3 <- my_data %>% filter(between(cyl, 3, 4))
data3

#filter across many columns using if_any() function

my_data2 <- iris
data4 <- select(filter(my_data2, if_any("Species", ~. == "versicolor")), Sepal.Length)
data4
data5 <- filter(my_data2, if_any("Species", ~. =="versicolor")) %>% select(Sepal.Length)
data5

#If we are giving the condition with the stored variable the R will take it as
#column name so while mentioning variable with condition use !!

condition <- "Sepal"
data6<- select(filter(my_data2, if_any(starts_with(!!condition))), 
               Species)
data6
#endswith helper
condition2 <- "Length"
data7<- select(filter(my_data2, if_any(ends_with(!!condition2))),
               Species)
data7

#here is_any() helper will work like 