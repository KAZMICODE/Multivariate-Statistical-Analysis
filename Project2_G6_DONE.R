#Project 2:  PCA in time series

#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
#I chose a dataset with information about percentage of coverage of vaccines  | 
#in one- year- olds, we decided to cover the window time frame of 2000 to 2020,|
#focusing in 4 vaccines: BCG(Tuberculosis), MCV1 (Measles), Pol3(Polio),       |
#DTP3 (Diphteria, tetanus , pertussis)                                         |
#we chose 4 countries: Brasil, China, India, and  Mexico                       |
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

library(CCA)
#READ THE DATA:
data<-global.vaccination.coverage_filters_ORIGINA


#CLEANING DATA
# Remove non needed columns (it contains the names of each row) 
#removing by choosing the index of column to be removed
datax <- subset(data, select = -c(2,5, 6, 7, 9, 11, 12, 13))
colnames(datax)[3] ="BCG(Tub)"
colnames(datax)[4] ="MCV1(Meas)"
colnames(datax)[5] ="Pol3(Polio)"
colnames(datax)[6] ="DTP3(Diph)"
head(datax,5)

#Cleaning the years
dataxy <- subset(datax,Year>1999 & Year<2021 )

#Now we have all the variables we are interested only from 2000 to 2020
#dropping the countries we are not interested in
dataxyz <- subset(dataxy,
                  Entity =='Brazil' | Entity =='China' |
                    Entity=='India' | Entity=='Mexico' )
rownames(dataxyz) <- 1:nrow(dataxyz)
#dataxyz is the data frame with years from 2000 to 2020,
# the 4 variables of interest and the 6 chosen countries
#head(dataxyz, 84)

#STANDARDIZING DATA

standardize_data <- function(data) {
  #--- Select the columns for standardization ---
  data_to_standardize <- data[, 3:6]
  #--- Standardize the data ---
  standardized_data <- scale(data_to_standardize)
  return(standardized_data)
}

data_brazil <- standardize_data(dataxyz[dataxyz$Entity == 'Brazil', ])
data_mexico <- standardize_data(dataxyz[dataxyz$Entity == 'Mexico', ])
data_china <- standardize_data(dataxyz[dataxyz$Entity == 'China', ])
data_india <- standardize_data(dataxyz[dataxyz$Entity == 'India', ])

rownames(data_brazil) <- 1:nrow(data_brazil)
rownames(data_china) <- 1:nrow(data_china)
rownames(data_india) <- 1:nrow(data_india)
rownames(data_mexico) <- 1:nrow(data_mexico)

#__#__#__#__#__#__#__#__#
#APPLYING THE METHOD 1  |
#__#__#__#__#__#__#__#__#

#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
#For all pairs of countries, perform Canonical Correlation Analysis (CCA)    | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

cca_brazil_china <- cancor(data_brazil, data_china)
cca_brazil_india <- cancor(data_brazil, data_india)
cca_brazil_mexico <- cancor(data_brazil, data_mexico)
cca_china_india <- cancor(data_china, data_india)
cca_mexico_china <- cancor(data_mexico, data_china)
cca_mexico_india <- cancor(data_mexico, data_india)

#Showing results of CCA for each pair
cca_brazil_china
cca_brazil_india
cca_brazil_mexico
cca_china_india
cca_mexico_china
cca_mexico_india

#______________________________________________________________________________
# For CCA, the pairs of canonical variates are calculated from the original   |
# datasets, using the canonical coeff. (also called weights or loadings)      |
# These are the coefficients that maximize the correlation between the        |
# canonical variates from each set of variables.                              |

# Given our CCA results, let's identify the pair with maximum canonical corr. |
# Once we identify the max corr., we can compute U1 and U2                    |
#______________________________________________________________________________

#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
#Comparing the canonical correlations.                                        | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

#Each max correlation at each CCA is the 1st value of cor, we put these values
#in a matrix with the name of the CCA pair it belongs to
#And get the pair with the maximal correlation.

maxcor_pairsCCA <-matrix(c(cca_brazil_china$cor[1],cca_brazil_india$cor[1],
                    cca_brazil_mexico$cor[1],cca_china_india$cor[1],
                    cca_mexico_china$cor[1],cca_mexico_india$cor[1],
                    'Brazil','Brazil','Brazil','China','Mexico','Mexico',
                    'China','India','Mexico','India','China','India'), nrow=6)

print(paste('The maximum correlation is ',max(maxcor_pairsCCA[,1]), 
'corresponding to pair', (maxcor_pairsCCA[which.max(maxcor_pairsCCA[,1]),2]),
' - ',(maxcor_pairsCCA[which.max(maxcor_pairsCCA[,1]),3])))


#______________________________________________________________________________
# Now, we have the first pair "Brazil-Mexico" with the highest canonical corr |
# with this pair we will compute the variates U1, U2.                         |
# Recall: Canonical Variate<- linear combinations of the original variables   |
# in each dataset, chosen such that their correlation is maximized            |
#______________________________________________________________________________

#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
#Canonical variates for "Brazil-Mexico"                                        | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

#Extract the Canonical Coefficients >stored in the $xcoef for Brazil and 
#  >stored in the $xcoef for Brazil and $ycoef for Mexico and are components
#   of the CCA result
#These coefficients determine how much each original variable contributes to 
#the canonical variates.

coef_brazil <- cca_brazil_mexico$xcoef[,1]
coef_mexico <- cca_brazil_mexico$ycoef[,1]

#The canonical variates U1 (Brazil) and U2 (Mexico) are computed using 
#the canonical coefficients:
U1 <- as.matrix(data_brazil) %*% coef_brazil
U2 <- as.matrix(data_mexico) %*% coef_mexico

#PRINTING U1 AND U2
#U1
#U2

#Examining the canonical variates with a scatter plot and a regression line
plot(U1, U2, main="Scatterplot of U1 vs. U2", 
     xlab="U1 (Brazil)", ylab="U2 (Mexico)")
abline(lm(U2 ~ U1), col="red")  


#______________________________________________________________________________
# Given the high can corr, there should be a strong linear trend in the       |
# scatterplot, meaning that the patterns of vaccination rates in Brazil (U1)  |
# and in Mexico(U2) have a high degree of shared structure. In the plot we    |
# can see a conglomeration of observations, this could correspond to certain  |
# years with noticeable deviations in vaccination trends compared to other    |
# years, or due an influence by the nature of the con variates themselves.    |

# Once we have U1, U2 related to Brazil and Mexico, we go back to our corr to |
# look for U3 and U4 which will be given by the second - highest canonical cor|
#______________________________________________________________________________

#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
#Comparing the canonical correlations for U3 and U4                           | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

#"delete" corr mex-brazil
maxcor_pairsCCA[which.max(maxcor_pairsCCA[,1]),1]=0
maxcor_pairsCCA
print(paste('The second maximum correlation is ',max(maxcor_pairsCCA[,1]), 
            'corresponding to pair', (maxcor_pairsCCA[which.max(maxcor_pairsCCA[,1]),2]),
            ' - ',(maxcor_pairsCCA[which.max(maxcor_pairsCCA[,1]),3])))
#nor Mexico nor Brazil involved so we can proceed

#______________________________________________________________________________
# Now, we have the second pair "China-India" with the 2nd highest cancorr     |
# with this pair we will compute the variates U3, U4.                         |
#______________________________________________________________________________

#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
#Canonical variates for "China-India"                                        | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

coef_china <- cca_china_india$xcoef[,1]
coef_india <- cca_china_india$ycoef[,1]

#The canonical variates U3 (China) and U4 (India) are computed using 
#the canonical coefficients

U3 <- as.matrix(data_china)%*% coef_china
U4 <- as.matrix(data_india)%*% coef_india

#PRINTING U3 AND U4
#U3
#U4

#Examining the canonical variates with a scatter plot and a regression line
plot(U3, U4, main="Scatterplot of U3 vs. U4", 
     xlab="U3 (China)", ylab="U4 (India)")
abline(lm(U4 ~ U3), col="red")

#______________________________________________________________________________
# Given the 2nd high can corr, it is expected that the points in the scatterp.|
# will lie close the regression line-> strong linear relation                 |
#

# Now we will get the multivariate time series Ut where each column is...     |
# and then we will compute the cross-autocovariance for each given pair (j,j')|
# at lag k, to obtain the cross-autocorrelation matriz for each pair          |
#______________________________________________________________________________

#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
#Construction of the multivariate time series Ut                              | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

Ut <- cbind(U1, U2, U3, U4)
#print
Ut

#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
#TO VISUALIZE THE PATTERNS OF THE COUNTRIES IN Ut                             | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*
# REMEMBER : Each variate represents a direction in the data space that 
#maximizes the correlation between the sets of variables from the different 
#countries.
 
library(ggplot2)
library(reshape2)

# Rename the columns to country names
colnames(Ut) <- c("Brazil", "Mexico", "China", "India")

# Convert Ut to a data frame ady
Ut <- as.data.frame(Ut)

# Create the time variable
time <- 1:nrow(Ut)

# Reshape the data for plotting with ggplot2
Ut_long <- melt(data.frame(time, Ut), id.vars = 'time')

# Now plot with ggplot2, using scale_color_manual to rename the legend labels
ggplot(Ut_long, aes(x = time, y = value, colour = variable)) +
  geom_line() +
  labs(title = 'METHOD 1 :Canonical Variate Scores Over Time',
       x = 'Time',
       y = 'Canonical Variate Score') +
  theme_minimal() +
  scale_color_viridis_d() +
  scale_color_manual(values = c("Brazil" = "yellow2", "Mexico" = "black", "China" = "red", "India" = "blue"),
                     labels = c("Brazil", "Mexico", "China", "India"))
#
#______________________________________________________________________________
#Each canonical variate captures certain characteristics of the data, and these|
#scores are a measure of how much each observation (year) is influenced by     |
#those characteristics.                                                        |
#______________________________________________________________________________

#__#__#__#__#__#__#__#__#
#APPLYING THE METHOD 2  |
#__#__#__#__#__#__#__#__#


#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
# Computing U_brazil                                                           | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

# Perform CCA between Brazil and India
cca_brazil_india <- cancor(data_brazil, data_india)
# Extract the second canonical variate for Brazil
U_brazil_india <- as.matrix(data_brazil) %*% cca_brazil_india$xcoef[, 2]
U_brazil_india
# Perform CCA between Brazil and Mexico
cca_brazil_mexico <- cancor(data_brazil, data_mexico)
# Extract the second canonical variate for Brazil
U_brazil_mexico <- as.matrix(data_brazil) %*% cca_brazil_mexico$xcoef[, 2]

# Perform CCA between Brazil and China
cca_brazil_china <- cancor(data_brazil, data_china)
# Extract the second canonical variate for Brazil
U_brazil_china <- as.matrix(data_brazil) %*% cca_brazil_china$xcoef[, 2]

# Sum up the second canonical variates to get U_brazil
U_brazil <- U_brazil_india + U_brazil_mexico + U_brazil_china
U_brazil

#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
# Computing  U_china                                                          | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

# Perform CCA between China and India
cca_china_india <- cancor(data_china, data_india)
# Extract the second canonical variate for China
U_china_india <- as.matrix(data_china) %*% cca_china_india$xcoef[, 2]
U_china_india
# Perform CCA between China and Mexico
cca_china_mexico <- cancor(data_china , data_mexico)
# Extract the second canonical variate for China
U_china_mexico <- as.matrix(data_china) %*% cca_china_mexico$xcoef[, 2]

# Perform CCA between Brazil and China
cca_china_brazil <- cancor(data_brazil, data_china)
# Extract the second canonical variate for China
U_china_brazil <- as.matrix(data_china) %*% cca_china_brazil$xcoef[, 2]

# Sum up the second canonical variates to get U_China
U_china <- U_china_india + U_china_mexico + U_china_brazil
U_china

#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
# Computing  U_india                                                          | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

# Perform CCA between Brazil and India
cca_india_brazil <- cancor(data_india, data_brazil)
# Extract the second canonical variate for India
U_india_brazil <- as.matrix(data_india) %*% cca_india_brazil$xcoef[, 2]
U_india_brazil
# Perform CCA between India and Mexico
cca_india_mexico <- cancor(data_india, data_mexico)
# Extract the second canonical variate for India
U_india_mexico <- as.matrix(data_india) %*% cca_brazil_mexico$xcoef[, 2]

# Perform CCA between India and China
cca_india_china <- cancor(data_india, data_china)
# Extract the second canonical variate for India
U_india_china <- as.matrix(data_india) %*% cca_india_china$xcoef[, 2]

# Sum up the second canonical variates to get U_India
U_india <- U_india_china + U_india_mexico + U_india_brazil
U_india


#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
# Computing  U_mexico                                                          | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
# Perform CCA between Mexico and India
cca_mexico_india <- cancor(data_mexico, data_india)
# Extract the second canonical variate for Mexico
U_mexico_india <- as.matrix(data_mexico) %*% cca_mexico_india$xcoef[, 2]
U_mexico_india
# Perform CCA between Brazil and Mexico
cca_mexico_china <- cancor(data_mexico , data_china)
# Extract the second canonical variate for Mexico
U_mexico_china <- as.matrix(data_mexico) %*% cca_mexico_china$xcoef[, 2]

# Perform CCA between Mexico and China
cca_mexico_brazil <- cancor(data_mexico, data_brazil)
# Extract the second canonical variate for Mexico
U_mexico_brazil <- as.matrix(data_mexico) %*% cca_mexico_brazil$xcoef[, 2]

# Sum up the second canonical variates to get U_Mexico
U_mexico <- U_mexico_brazil + U_mexico_china + U_mexico_india
U_mexico

#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
# Now to compute U_j (which is said as Y_j in paper ) by concatenating the    |
# values of U_brazil, U_india, U_mexico, and U_china                           |
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

# Once we have U_brazil, U_india, U_mexico, and U_china
U_j <- cbind(U_brazil, U_mexico, U_china, U_india)
U_j
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
#TO VISUALIZE THE PATTERNS OF THE COUNTRIES IN U_j                             | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

# Rename the columns to country names
colnames(U_j) <- c("Brazil", "Mexico", "China", "India")

U_j <- as.data.frame(U_j) 
time <- 1:nrow(U_j)  

# Reshape the data for plotting with ggplot2
U_j_long <- reshape2::melt(data.frame(time, U_j), id.vars = 'time')
#U_j_long
# Now plot with ggplot2
ggplot(U_j_long, aes(x = time, y = value, colour = variable)) +
  geom_line() +
  labs(title = 'METHOD 2: Canonical Variate Scores Over Time',
       x = 'Time',
       y = 'Canonical Variate Score') +
  theme_minimal() +
  scale_color_manual(values = c("Brazil" = "yellow2", "Mexico" = "black", "China" = "red", "India" = "blue"),
                     labels = c("Brazil", "Mexico", "China", "India"))
#
#______________________________________________________________________________
#The magnitude of the fluctuations indicates how volatile the canonical variate|
#scores are over time for each country. Large swings could suggest periods of |
#significant change in the vaccine data characteristics captured by that variate|                                                        |
#______________________________________________________________________________


#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
#Function for CROSS-AUTOCOVARIANCE                                            | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

compute_cross_autocovariance <- function(Ut, j, j_prime, k) {
  N <- nrow(Ut)
  mean_j <- mean(Ut[, j])
  mean_j_prime <- mean(Ut[, j_prime])
  
  gamma_hat <- sum((Ut[1:(N-k), j] - mean_j) * 
                     (Ut[(k+1):N, j_prime] - mean_j_prime)) / N
  return(gamma_hat)
}

#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
#Function for CROSS-AUTOCORRELATION                                           | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

compute_cross_autocorrelation <- function(Ut, j, j_prime, k) {
  gamma_hat_j_j_prime_k <- compute_cross_autocovariance(Ut, j, j_prime, k)
  gamma_hat_j_j_0 <- compute_cross_autocovariance(Ut, j, j, 0)
  gamma_hat_j_prime_j_prime_0 <- compute_cross_autocovariance(Ut, j_prime, j_prime, 0)
  
  rho_hat <- gamma_hat_j_j_prime_k / sqrt(gamma_hat_j_j_0 * gamma_hat_j_prime_j_prime_0)
  return(rho_hat)
}

#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
# Single function  to compute cross-autocorrelation matrix for a given lag k  |                                       | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

# Compute the cross-autocorrelation matrix for a given lag k
compute_rho_k <- function(Ut, k) {
  # Number of canonical variates (columns in Ut)
  p <- ncol(Ut)
  
  # Initialize the cross-autocorrelation matrix
  rho_hat_k <- matrix(0, p, p)
  
  # Fill in the matrix
  for(j in 1:p) {
    for(j_prime in 1:p) {
      rho_hat_k[j, j_prime] <- compute_cross_autocorrelation(Ut, j, j_prime, k)
    }
  }
  
  return(rho_hat_k)
}

#______________________________________________________________________________
# Now, we can use these functions to compute the cross-autocorrelation        |
# matrices for matrix ρ^(k) for a specific lag k by looping over all possible |
# pair (j,j') and filling in the matrix.                                      |
#______________________________________________________________________________


#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
# Computing  Cross-autocorrelation matrices by METHOD 1                                                         | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

rho_0 <- compute_rho_k(Ut, 0)
rho_1 <- compute_rho_k(Ut, 1)
rho_2 <- compute_rho_k(Ut, 2)
df_rho_0 <- as.data.frame(rho_0)
df_rho_1 <- as.data.frame(rho_1)
df_rho_2 <- as.data.frame(rho_2)

#______________________________________________________________________________
# FOR PRESENTATION IN SLIDES
library(knitr)
library(dplyr)
#install.packages("kableExtra")
library(kableExtra)

# Format the matrix to two decimal places
df_rho_0_formatted <- round(df_rho_0, digits = 2)
#df_rho_0_formatted
# Create a table
df_rho_0_table <- kable(df_rho_0_formatted, caption = "Cross-autocorrelation matrix at lag k=0 {METHOD_1}", digits = 2, format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
df_rho_0_table

# For K = 1 , "df_rho_1"
df_rho_1_formatted <- round(df_rho_1, digits = 2)
#df_rho_1_formatted
# Create a table
df_rho_1_table <- kable(df_rho_1_formatted, caption = "Cross-autocorrelation matrix at lag k=1 {METHOD_1}", digits = 2, format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
df_rho_1_table

# For K = 2 , "df_rho_2"
df_rho_2_formatted <- round(df_rho_2, digits = 2)
#df_rho_2_formatted
# Create a table
df_rho_2_table <- kable(df_rho_2_formatted, caption = "Cross-autocorrelation matrix at lag k=2 {METHOD_1}", digits = 2, format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
df_rho_2_table

#______________________________________________________________________________


#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
# Computing  Cross-autocorrelation matrices by METHOD 2                       |                                  | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

rho_hat_0 <- compute_rho_k(U_j, 0)
rho_hat_1 <- compute_rho_k(U_j, 1)
rho_hat_2 <- compute_rho_k(U_j, 2)
df_rho_hat_0 <- as.data.frame(rho_hat_0)
df_rho_hat_1 <- as.data.frame(rho_hat_1)
df_rho_hat_2 <- as.data.frame(rho_hat_2)

#______________________________________________________________________________
#These cross-autocorrelation matrices provide a measure of how much past values|
#how each of the canonical variates (which are combinations of our original    |
#vaccine data for the countries Brazil, India, Mexico, and China) is correlated| 
#with itself and with the other variates at different time lags.               |
#______________________________________________________________________________

#______________________________________________________________________________
# FOR PRESENTATION IN SLIDES

# Format the matrix to two decimal places
df_rho_hat_0_formatted <- round(df_rho_hat_0, digits = 2)
# Create a table
df_rho_hat_0_table <- kable(df_rho_hat_0_formatted, caption = "Cross-autocorrelation matrix at lag k=0 {METHOD_2}", digits = 2, format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
df_rho_hat_0_table

# For K = 1 , "df_rho_hat_1"
df_rho_hat_1_formatted <- round(df_rho_hat_1, digits = 2)
# Create a table
df_rho_hat_1_table <- kable(df_rho_hat_1_formatted, caption = "Cross-autocorrelation matrix at lag k=1 {METHOD_2}", digits = 2, format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
df_rho_hat_1_table

# For K = 2 , "df_rho_hat_2"
df_rho_hat_2_formatted <- round(df_rho_hat_2, digits = 2)
# Create a table
df_rho_hat_2_table <- kable(df_rho_hat_2_formatted, caption = "Cross-autocorrelation matrix at lag k=2 {METHOD_2}", digits = 2, format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
df_rho_hat_2_table

#______________________________________________________________________________


#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
# Checking the positive definiteness of each cross-autocorrelation function   |                                             | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*


#______________________________________________________________________________
#For lags k>0, the cross-autocorrelation matrix is not symmetric and may not be| 
#positive definite. This means it could have negative eigenvalues,  which is   |
#problematic for certain statistical analyses. A matrix needs to be positive   |
#definite for eigenvalue decomposition used in PCA, where only positive        |
#eigenvalues are meaningful.So we have two ways , either Negative eigenvalues  |
#of the correlation matrix are set to zero. This retains the structure of the  |
#matrix but ensures all eigenvalues are non-negative (APPRAOCH_1)  {OR}        |
#If we take the absolute values of the eigenvalues, (APPRAOCH_2)               |
#______________________________________________________________________________

#__#__#__#__#__#__#__#__#__#__#__#__#__#
#DEFINING A FUNCTION FOR 1st Approach |
#__#__#__#__#__#__#__#__#__#__#__#__#__#

make_positive_definite <- function(cor_matrix) {
  # Ensure the matrix is symmetric
  cor_matrix <- (cor_matrix + t(cor_matrix)) / 2
  
  # Perform eigenvalue decomposition
  eig <- eigen(cor_matrix, symmetric = TRUE)
  
  # Set negative eigenvalues to zero and keep the real parts
  eig$values <- pmax(0, Re(eig$values))
  
  # Construct the matrix using only the real parts of eigenvectors
  P <- Re(eig$vectors)
  V_prime <- diag(eig$values)
  
  # Calculate new matrix B
  B <- P %*% sqrt(V_prime)
  
  # Construct the estimated positive definite matrix
  A_hat <- B %*% t(B)
  
  # Return the adjusted matrix
  return(A_hat)
}

#__#__#__#__#__#__#__#__#__#__#__#__#__#
#DEFINING A FUNCTION FOR 2nd Approach |
#__#__#__#__#__#__#__#__#__#__#__#__#__#

make_positive_definite_abs <- function(cor_matrix) {
  # Perform eigenvalue decomposition
  eig <- eigen(cor_matrix)
  # Take the absolute values of the eigenvalues
  eig$values <- abs(eig$values)
  
  # Calculate A' using the absolute values of eigenvalues
  Lambda_prime <- diag(eig$values)
  A_prime <- eig$vectors %*% Lambda_prime %*% solve(eig$vectors)
  
  # Construct the diagonal scaling matrix D
  D_diag <- 1 / sqrt(diag(A_prime))
  D <- diag(D_diag)
  
  # Calculate the new positive definite matrix A_hat
  A_hat <- eig$vectors %*% D %*% t(eig$vectors)
  
  # Return the adjusted matrix
  return(A_hat)
}

#__#__#__#__#__#__#__#__#__#__#__#__#__#__#__#__#__#__#_
#DEFINING A FUNCTION FOR CHECKING POSITIVE DEFINITENESS|
#__#__#__#__#__#__#__#__#__#__#__#__#__#__#__#__#__#__#_

is_positive_definite <- function(mat) {
  eigenvalues <- eigen(mat)$values
  # Check if all eigenvalues are greater than a small positive threshold
  all(eigenvalues > .Machine$double.eps)
}

#____________________________________________________________
# TO CHECK FOR k=0                                          #
## And to check the originall one                           #                     
is_pd <- is_positive_definite(rho_hat_0)                    #
print(is_pd)                                                #
                                                            #
is_pd <- is_positive_definite(rho_0)                        #
print(is_pd)                                                #
                                                            #
## TO CHECK FOR k=1                                         #
                                                            #
is_pd <- is_positive_definite(rho_hat_1)                    #
print(is_pd)                                                #
is_pd <- is_positive_definite(rho_1)                        #
print(is_pd)                                                #
                                                            #
### TO CHECK FOR k=2                                        #
is_pd <- is_positive_definite(rho_hat_2)                    #
print(is_pd)                                                #
is_pd <- is_positive_definite(rho_2)                        #
print(is_pd)                                                #
                                                            #
## APPLYING FUNCTIONS                                       #
# FOR k=1                                                   #
rho_hat_1_pd <- make_positive_definite_abs(rho_hat_1)       #
rho_hat_1_pd <- make_positive_definite(rho_hat_1_pd)        #
rho_hat_1_pd_abs <- make_positive_definite_abs(rho_hat_1_pd)#
is_pd <- is_positive_definite(rho_hat_1_pd_abs)             #
print(is_pd)                                              #
## rho_1 is already positive definite                        #                                                #
                                                            #
# FOR k=2                                                   #
rho_hat_2_pd_abs <- make_positive_definite_abs(rho_hat_2)   #
is_pd <- is_positive_definite(rho_hat_2_pd_abs)             #
print(is_pd)                                                #
                                                            #
rho_2_ab <- make_positive_definite_abs(rho_2)               #
rho_2_abs <- make_positive_definite(rho_2_ab)               #
is_pd <- is_positive_definite(rho_hat_1_pd_abs)             #
print(is_pd)                                                #
#____________________________________________________________

#___________--_____________--_____________--_____________--___
#So the finalized cross-autocorrelation matrix for k=1 are    |
rho_hat_1_pd_abs #For U_j (MEHTOD_2)                          |
rho_1            #For Ut  (METHOD_1)                          |
#___________--_____________--_____________--_____________--___


#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
#APPLYING PCA ON THE CROSS-AUTOCORREALTION MATRICES OBTAINED                   | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

#______________________________________________________________________________
# We are to follow the  Eq. (3.2), to obtain the principal components PCν (t,s)|
#ν = 1,..., p, for each series s = 1,..., S, and evaluated at each time        |
#t = 1,..., N. So for this we need a 3D array of the U_j and Ut having data of |
# the 4 varaibles for each 4 country over the time period of 21 years          |
#______________________________________________________________________________

#__________-__________-__________-___
# Constructing 3D array for my data |
#__________-__________-__________-__

# Number of years and variables
num_years <- 21     # Over the time period {2000 to 2020}
num_variables <- 4  # Since we have four vaccines 
num_countries <- 4  # Since we have four countries

# Create an empty 3D array
Uf <- array(0, dim = c(num_years, num_variables, num_countries))

# Assign the data for each country to a slice of the array
Uf[,,1] <- as.matrix(data_brazil)
Uf[,,2] <- as.matrix(data_mexico)
Uf[,,3] <- as.matrix(data_china)
Uf[,,4] <- as.matrix(data_india)
Uf
# Now Ut is a 3D array with dimensions [time, variables, series]

#__________-__________-__________-_____
# Getting Results by Method_1 approach |
#__________-__________-__________-_____

# Eigen decomposition of rho_1
eigen_results <- eigen(rho_1)
eigenvectors <- eigen_results$vectors

# Number of principal components, time points, and series
num_pcs <- ncol(eigenvectors)
num_times <- dim(Uf)[1]
num_series <- dim(Uf)[3]

# Initialize the PCs array with the appropriate dimensions
PCs <- array(0, dim = c(num_times, num_pcs, num_series))

# Compute the Principal Components for each series at each time t
for (s in 1:num_series) {
  for (nu in 1:num_pcs) {
    for (t in 1:num_times) {
      # Sum the product of variables at time t and series s with the nu-th eigenvector
      PCs[t, nu, s] <- sum(Uf[t, , s] * eigenvectors[, nu])
    }
  }
}
PCs
# PCs is now a 3D array with principal component scores for each time t, 
#component nu, and series s
#______________________________________________________________________________
#COMMENT ON ABOVE METHOD USED
#There are three nested loops:
# The outermost loop iterates over the series s.
# The middle loop iterates over the principal components nu.
# The innermost loop iterates over the time points t.
# Within the innermost loop, for each series s, we compute the principal
# component score for component nu at time t by summing the products of the 
# t-th observation's variables for series s and the nu-th eigenvector.
#______________________________________________________________________________


#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
# Interpretations and Visualizations                                          | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

#___-___-___-___-___-___-
# THE 2-D Plot          |
#___-___-___-___-___-___-

# We have an array countries containing the country names
countries <- c("Brazil", "Mexico", "China", "India")

# Flatten the 3D PCs array to a 2D matrix for plotting
# We only need the first two PCs
PC1 <- PCs[, 1, ]  # First principal component
PC2 <- PCs[, 2, ]  # Second principal component

# Convert 3D slices to 2D matrix
PC1_flat <- as.vector(PC1)
PC2_flat <- as.vector(PC2)

# Create a vector that repeats the country names for each time point
country_rep <- rep(countries, each = nrow(PCs[,,1]))

# Combine into a data frame for ggplot
pc_data <- data.frame(Country = factor(country_rep), PC1 = PC1_flat, PC2 = PC2_flat)
#Defining the colours
country_colors <- c("Brazil" = "yellow2", "Mexico" = "black", "China" = "red", "India" = "blue")

# Use ggplot2 to create the scatter plot with colors
library(ggplot2)
ggplot(pc_data, aes(x = PC1, y = PC2, color = Country)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = country_colors) +  # Set custom colors
  labs(title = "Method 1 : PC1 × PC2 - all times and all series, lag k = 1",
       x = "Principal Component 1",
       y = "Principal Component 2") 
  #scale_color_brewer(palette = "Set1")  # Use a color palette that is distinct

#______________________________________________________________________________
#showing plots of the first and second principal components for all series   |
#and all times without tracking the time dependencies                        |
#______________________________________________________________________________

#___-___-___-___-___-___-
# THE 3-D Plot          |
#___-___-___-___-___-___-

# Load the library
library(scatterplot3d)

# Create a vector for the country codes (numeric)
country_codes <- rep(1:4, each = nrow(PCs))
# Flatten the 3D array into a 2D matrix for plotting
PC1 <- c(PCs[,,1][,1], PCs[,,2][,1], PCs[,,3][,1], PCs[,,4][,1])
PC2 <- c(PCs[,,1][,2], PCs[,,2][,2], PCs[,,3][,2], PCs[,,4][,2])
Time <- rep(1:nrow(PCs), 4)

# Create a data frame for plotting
PC_data <- data.frame(Time, PC1, PC2, Country = factor(country_codes))

# Define colors for each country
colors <- c("yellow2", "black","red", "blue")
names(colors) <- c("Brazil", "Mexico", "China", "India")

# Create the 3D scatter plot 
s3d <- scatterplot3d(
  z = PC_data$PC1,
  y = PC_data$PC2,
  x = PC_data$Time,
  type = "n",  # Do not plot points, use this plot as a base to add lines
  main = " METHOD 1 :PC1 × PC2 over time, all series, lag k = 1 ",  # Add your title here
  xlab = "Time",
  ylab = "PC2",
  zlab = "PC1",
  cex.main = 0.75 #Changing the size of title
)
# Assume: Brazil (solid), Mexico (dashed), China (dotted), India (dotdash)
# Define line types for each country
line_types <- c("solid", "dashed", "dotted", "dotdash")  
line_widths <- c(2, 2, 2, 2) 
# Plot lines for each country
for (i in 1:length(unique(PC_data$Country))) {
  country_data <- subset(PC_data, Country == i)
  s3d$points3d(
    x = country_data$Time,
    y = country_data$PC2,
    z = country_data$PC1,
    type = "l",  # Specify line type
    lty = i,  # Use the line type corresponding to the country (1 to 4)
    col = colors[i],
    lwd = line_widths[i]  # Set the line width
  )
}

# Add a legend manually
legend("topright", inset = c(-0.05, 0), 
       legend = names(colors),
       col = colors,
       lty = 1:4  # Specify line types corresponding to 'line_types'
)

#__________-__________-__________-_____
# Getting Results by Method_2 approach |
#__________-__________-__________-_____


# Eigen decomposition of rho_hat_1_pd_abs
eigen_results <- eigen(rho_hat_1_pd_abs)
eigenvectors <- eigen_results$vectors

# Number of principal components, time points, and series
num_pcs <- ncol(eigenvectors)
num_times <- dim(Uf)[1]
num_series <- dim(Uf)[3]

# Initialize the PCs array with the appropriate dimensions
PCs <- array(0, dim = c(num_times, num_pcs, num_series))

# Compute the Principal Components for each series at each time t
for (s in 1:num_series) {
  for (nu in 1:num_pcs) {
    for (t in 1:num_times) {
      # Sum the product of variables at time t and series s with the nu-th eigenvector
      PCs[t, nu, s] <- sum(Uf[t, , s] * eigenvectors[, nu])
    }
  }
}
PCs


#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
# Interpretations and Visualizations                                          | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

#___-___-___-___-___-___-
# THE 2-D Plot          |
#___-___-___-___-___-___-

# We have an array countries containing the country names
countries <- c("Brazil", "Mexico", "China", "India")

# Flatten the 3D PCs array to a 2D matrix for plotting
# We only need the first two PCs
PC1 <- PCs[, 1, ]  # First principal component
PC2 <- PCs[, 2, ]  # Second principal component

# Convert 3D slices to 2D matrix
PC1_flat <- as.vector(PC1)
PC2_flat <- as.vector(PC2)

# Create a vector that repeats the country names for each time point
country_rep <- rep(countries, each = nrow(PCs[,,1]))

# Combine into a data frame for ggplot
pc_data <- data.frame(Country = factor(country_rep), PC1 = PC1_flat, PC2 = PC2_flat)
#Defining the colours
country_colors <- c("Brazil" = "yellow2", "Mexico" = "black", "China" = "red", "India" = "blue")

# Use ggplot2 to create the scatter plot with colors
library(ggplot2)
ggplot(pc_data, aes(x = PC1, y = PC2, color = Country)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = country_colors) +  # Set custom colors
  labs(title = "Method 2 : PC1 × PC2 - all times and all series, lag k = 1",
       x = "Principal Component 1",
       y = "Principal Component 2") 
  #scale_color_brewer(palette = "Set1")  # Use a color palette that is distinct

#___-___-___-___-___-___-
# THE 3-D Plot          |
#___-___-___-___-___-___-

# Load the library
library(scatterplot3d)

# Create a vector for the country codes (numeric)
country_codes <- rep(1:4, each = nrow(PCs))
# Flatten the 3D array into a 2D matrix for plotting
PC1 <- c(PCs[,,1][,1], PCs[,,2][,1], PCs[,,3][,1], PCs[,,4][,1])
PC2 <- c(PCs[,,1][,2], PCs[,,2][,2], PCs[,,3][,2], PCs[,,4][,2])
Time <- rep(1:nrow(PCs), 4)

# Create a data frame for plotting
PC_data <- data.frame(Time, PC1, PC2, Country = factor(country_codes))

# Define colors for each country
colors <- c("yellow2", "black","red", "blue")
names(colors) <- c("Brazil", "Mexico", "China", "India")

# Create the 3D scatter plot 
s3d <- scatterplot3d(
  z = PC_data$PC1,
  y = PC_data$PC2,
  x = PC_data$Time,
  type = "n",  # Do not plot points, use this plot as a base to add lines
  main = " METHOD 2 :PC1 × PC2 over time, all series, lag k = 1 ",  # Add your title here
  xlab = "Time",
  ylab = "PC2",
  zlab = "PC1",
  cex.main = 0.75 #Changing the size of title
)
# Assume: Brazil (solid), Mexico (dashed), China (dotted), India (dotdash)
# Define line types for each country
line_types <- c("solid", "dashed", "dotted", "dotdash")  
line_widths <- c(2, 2, 2, 2) 
# Plot lines for each country
for (i in 1:length(unique(PC_data$Country))) {
  country_data <- subset(PC_data, Country == i)
  s3d$points3d(
    x = country_data$Time,
    y = country_data$PC2,
    z = country_data$PC1,
    type = "l",  # Specify line type
    lty = i,  # Use the line type corresponding to the country (1 to 4)
    col = colors[i],
    lwd = line_widths[i]  # Set the line width
  )
}

# Add a legend manually
legend("topright", inset = c(-0.05, 0), 
       legend = names(colors),
       col = colors,
       lty = 1:4  # Specify line types corresponding to 'line_types'
)



#___________--_____________--_____________--_____________--_____________--_____
# COMPUTATION FOR LAG K=0 , TO COMPARE THE RESULTS AND DIFFERENCES SEEN       |
#___________--_____________--_____________--_____________--_____________--_____

#__________-__________-__________-_____
# Getting Results by Method_1 approach |
#__________-__________-__________-_____

# Eigen decomposition of rho_0
eigen_results <- eigen(rho_0)
eigenvectors <- eigen_results$vectors

# Number of principal components, time points, and series
num_pcs <- ncol(eigenvectors)
num_times <- dim(Uf)[1]
num_series <- dim(Uf)[3]

# Initialize the PCs array with the appropriate dimensions
PCs <- array(0, dim = c(num_times, num_pcs, num_series))

# Compute the Principal Components for each series at each time t
for (s in 1:num_series) {
  for (nu in 1:num_pcs) {
    for (t in 1:num_times) {
      # Sum the product of variables at time t and series s with the nu-th eigenvector
      PCs[t, nu, s] <- sum(Uf[t, , s] * eigenvectors[, nu])
    }
  }
}
PCs
# PCs is now a 3D array with principal component scores for each time t, 
#component nu, and series s


#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
# Interpretations and Visualizations                                          | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

#___-___-___-___-___-___-
# THE 2-D Plot          |
#___-___-___-___-___-___-

# We have an array countries containing the country names
countries <- c("Brazil", "Mexico", "China", "India")

# Flatten the 3D PCs array to a 2D matrix for plotting
# We only need the first two PCs
PC1 <- PCs[, 1, ]  # First principal component
PC2 <- PCs[, 2, ]  # Second principal component

# Convert 3D slices to 2D matrix
PC1_flat <- as.vector(PC1)
PC2_flat <- as.vector(PC2)

# Create a vector that repeats the country names for each time point
country_rep <- rep(countries, each = nrow(PCs[,,1]))

# Combine into a data frame for ggplot
pc_data <- data.frame(Country = factor(country_rep), PC1 = PC1_flat, PC2 = PC2_flat)
country_colors <- c("Brazil" = "yellow2", "Mexico" = "black", "China" = "red", "India" = "blue")
# Use ggplot2 to create the scatter plot with colors
library(ggplot2)
ggplot(pc_data, aes(x = PC1, y = PC2, color = Country)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = country_colors) +  # Set custom colors
  labs(title = "Method 1 : PC1 × PC2 - all times and all series, lag k = 0",
       x = "Principal Component 1",
       y = "Principal Component 2") 
  #scale_color_brewer(palette = "Set1")  # Use a color palette that is distinct

#______________________________________________________________________________
#showing plots of the first and second principal components for all series   |
#and all times without tracking the time dependencies                        |
#______________________________________________________________________________

# Load the library
library(scatterplot3d)

# Create a vector for the country codes (numeric)
country_codes <- rep(1:4, each = nrow(PCs))
# Flatten the 3D array into a 2D matrix for plotting
PC1 <- c(PCs[,,1][,1], PCs[,,2][,1], PCs[,,3][,1], PCs[,,4][,1])
PC2 <- c(PCs[,,1][,2], PCs[,,2][,2], PCs[,,3][,2], PCs[,,4][,2])
Time <- rep(1:nrow(PCs), 4)

# Create a data frame for plotting
PC_data <- data.frame(Time, PC1, PC2, Country = factor(country_codes))

# Define colors for each country
colors <- c("yellow2", "black","red", "blue")
names(colors) <- c("Brazil", "Mexico", "China", "India")

# Create the 3D scatter plot 
s3d <- scatterplot3d(
  z = PC_data$PC1,
  y = PC_data$PC2,
  x = PC_data$Time,
  type = "n",  # Do not plot points, use this plot as a base to add lines
  main = " METHOD 1 :PC1 × PC2 over time, all series, lag k = 0 ",  # Add your title here
  xlab = "Time",
  ylab = "PC2",
  zlab = "PC1",
  cex.main = 0.75 #Changing the size of title
)
# Assume: Brazil (solid), Mexico (dashed), China (dotted), India (dotdash)
# Define line types for each country
line_types <- c("solid", "dashed", "dotted", "dotdash")  
line_widths <- c(2, 2, 2, 2) 
# Plot lines for each country
for (i in 1:length(unique(PC_data$Country))) {
  country_data <- subset(PC_data, Country == i)
  s3d$points3d(
    x = country_data$Time,
    y = country_data$PC2,
    z = country_data$PC1,
    type = "l",  # Specify line type
    lty = i,  # Use the line type corresponding to the country (1 to 4)
    col = colors[i],
    lwd = line_widths[i]  # Set the line width
  )
}

# Add a legend manually
legend("topright", inset = c(-0.05, 0), 
       legend = names(colors),
       col = colors,
       lty = 1:4  # Specify line types corresponding to 'line_types'
)


#__________-__________-__________-_____
# Getting Results by Method_2 approach |
#__________-__________-__________-_____


# Eigen decomposition of rho_hat_0
eigen_results <- eigen(rho_hat_0)
eigenvectors <- eigen_results$vectors
#eigenvectors
# Number of principal components, time points, and series
num_pcs <- ncol(eigenvectors)
num_times <- dim(Uf)[1]
num_series <- dim(Uf)[3]

# Initialize the PCs array with the appropriate dimensions
PCs <- array(0, dim = c(num_times, num_pcs, num_series))

# Compute the Principal Components for each series at each time t
for (s in 1:num_series) {
  for (nu in 1:num_pcs) {
    for (t in 1:num_times) {
      # Sum the product of variables at time t and series s with the nu-th eigenvector
      PCs[t, nu, s] <- sum(Uf[t, , s] * eigenvectors[, nu])
    }
  }
}
PCs


#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
# Interpretations and Visualizations                                          | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_

#___-___-___-___-___-___-
# THE 2-D Plot          |
#___-___-___-___-___-___-

# We have an array countries containing the country names
countries <- c("Brazil", "Mexico", "China", "India")

# Flatten the 3D PCs array to a 2D matrix for plotting
# We only need the first two PCs
PC1 <- PCs[, 1, ]  # First principal component
PC2 <- PCs[, 2, ]  # Second principal component

# Convert 3D slices to 2D matrix
PC1_flat <- as.vector(PC1)
PC2_flat <- as.vector(PC2)

# Create a vector that repeats the country names for each time point
country_rep <- rep(countries, each = nrow(PCs[,,1]))

# Combine into a data frame for ggplot
pc_data <- data.frame(Country = factor(country_rep), PC1 = PC1_flat, PC2 = PC2_flat)
country_colors <- c("Brazil" = "yellow2", "Mexico" = "black", "China" = "red", "India" = "blue")

# Use ggplot2 to create the scatter plot with colors
library(ggplot2)
ggplot(pc_data, aes(x = PC1, y = PC2, color = Country)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = country_colors) +  # Set custom colors
  labs(title = "Method 2 : PC1 × PC2 - all times and all series, lag k = 0",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme(legend.position = "bottom") +  # Place the legend at the bottom
  guides(color = guide_legend(nrow = 1))  # Display legend horizontally

  #scale_color_brewer(palette = "Set1")  # Use a color palette that is distinct

#___-___-___-___-___-___-
# THE 3-D Plot          |
#___-___-___-___-___-___-


# Create the 3D scatter plot 
s3d <- scatterplot3d(
  z = PC_data$PC1,
  y = PC_data$PC2,
  x = PC_data$Time,
  type = "n",  # Do not plot points, use this plot as a base to add lines
  main = " METHOD 2 :PC1 × PC2 over time, all series, lag k = 0 ",  # Add your title here
  xlab = "Time",
  ylab = "PC2",
  zlab = "PC1",
  cex.main = 0.75 #Changing the size of title
)
# Assume: Brazil (solid), Mexico (dashed), China (dotted), India (dotdash)
# Define line types for each country
line_types <- c("solid", "dashed", "dotted", "dotdash")  
colors <- c("yellow2", "black", "red", "blue")
names(colors) <- c("Brazil", "Mexico", "China", "India")
line_widths <- c(2, 2, 2, 2) 
# Plot lines for each country
for (i in 1:length(unique(PC_data$Country))) {
  country_data <- subset(PC_data, Country == i)
  s3d$points3d(
    x = country_data$Time,
    y = country_data$PC2,
    z = country_data$PC1,
    type = "l",  # Specify line type
    lty = i,  # Use the line type corresponding to the country (1 to 4)
    col = colors[i],
    lwd = line_widths[i]  # Set the line width
  )
}

# Add a legend manually
legend("topright", inset = c(-0.05, 0), 
       legend = names(colors),
       col = colors,
       lty = 1:4  # Specify line types corresponding to 'line_types'
)

#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
#                                THANK_YOU                                      | 
#__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*__*_
