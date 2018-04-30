#Assume "group" is the name of your grouping variable in your data set
# (so you have a few rows for each group)
# To create a vector of unique group values (when data contains many rows per group)

unique(data1$group)

# You can use the vector above to create samples, for example

#Let's assume we care about the following sample of 4 groups

subgroup = c(234, 241, 180,  39)


# The next line code creates a plot of "y" against "x" by group 
# only for the group members in your subgroup

filter(data1, group %in% subgroup) %>% 
  ggplot(aes(x=x, y=y)) + geom_point() + facet_grid(~ group)

# To create a smaller data set with unique rows per "group"
# (the top row of each "group" in the dataset is selected)
# you can also use the following code

data1[!duplicated(data1$group),]

# To obtain regression coefficients per group here are two methods
# Here I'll use "y ~ x1 + x2 +x3" as my model
# Method 1: Using a "for" loop

group.id = unique(data1$group)
n = length(group.id) #number of different categories
p = #fill this value correspond to the number of 
    # coefficient estimates per regression
coef.mat = matrix(NA, nrow = n, ncol = p )
for (i in 1:n){
  datai = data1[data1$group == group.id[i],]
  coef.mat[i,] = as.vector(coef(lm(y ~ x1 + x2 +x3, data = datai)))
}
# The object "coef.mat" is a matrix where each column correspond to 
# a given coefficient (e.g. first column are beta0 estimates)

# Method 2: Using command "lmList" from package lme4
m.list = lmList(y ~ x1 + x2 +x3|group, data1)
# the object "coef(m.list)" is equivalent to "coef.mat" from Method 1