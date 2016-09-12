\name{dietary_survey_IBS}
\alias{dietary_survey_IBS}
\docType{data}
\title{
Synthetic data using a dietary survey of patients	with irritable	bowel syndrome   (IBS)
}
\description{
The data are based on the article "A dietary survey	of patients	with irritable bowel syndrome". The mean and standard deviation of the table 1 (Foods perceived as causing or worsening irritable bowel syndrome symptoms in the IBS group and digestive symptoms in the healthy comparative group) were used to generate the synthetic data.
}
\usage{data(dietary_survey_IBS)}
\format{
  A data frame with 400 Instances and 43 attributes (including the class attribute, "class")
}
\details{
The predictors are: bread, wheat, pasta, breakfast_cereal, yeast, spicy_food, curry, chinese_takeaway, chilli, cabbage, onion, garlic, 
potatoes, pepper, vegetables_unspecified, tomato, beans_and_pulses, mushroom, fatty_foods_unspecified, sauces, chocolate,
fries, crisps, desserts, eggs, red_meat, processed_meat, pork, chicken, fish_shellfish, dairy_products_unspecified, cheese, 
cream, milk, fruit_unspecified, nuts_and_seeds, orange, apple, banana, grapes, alcohol, caffeine

The response variable ("class") consists of two groups: healthy-group (class == 0) vs. the IBS-patients (class == 1)
}
\references{
P. Hayes, C. Corish, E. O'Mahony, E. M. M. Quigley (May 2013). A dietary survey of patients with irritable bowel syndrome. 
Journal of Human Nutrition and Dietetics.
}
\examples{

data(dietary_survey_IBS)

X = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

y = dietary_survey_IBS[, ncol(dietary_survey_IBS)]
}
\keyword{datasets}
