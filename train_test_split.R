library("caTools")
library("carot")

# split data frame into train and test subset (0.8 to 0.2 ratio)
train_test_split = function(df) {
  set.seed(42)
  df = as.data.frame(df)
  sample = sample.split(df, SplitRatio = 0.8)
  train = subset(df, sample == TRUE)
  test  = subset(df, sample == FALSE)
  return(list(train, test))
}