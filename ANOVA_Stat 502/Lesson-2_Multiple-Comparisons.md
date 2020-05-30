-   [Comparison of Multiple Comparison
    Methods](#comparison-of-multiple-comparison-methods)
    -   [Fisher’s Least Significant Difference
        (LSD)](#fishers-least-significant-difference-lsd)
    -   [Tukey](#tukey)
    -   [Bonferroni](#bonferroni)
    -   [Scheffe](#scheffe)
    -   [Dunnett](#dunnett)

``` r
#Need for Dunnett distribution
library(nCDunnett)
```

    ## Warning: package 'nCDunnett' was built under R version 3.5.2

Stacked format of data from Lesson 1 was saved tor a csv file and read
into an R dataframe

``` r
df <- read.csv("StackedData.csv",header = F)
```

Define variables for use in subsequent steps

``` r
#Four treatments
NumOfGroups <- 4
#Number of pairwise comparisons is simply 4 choose 2
NumOfCompar <- factorial(NumOfGroups) / factorial(NumOfGroups-2) / factorial(2)
#Number of observations per treatment is 6
NumofRepl <- 6
#Set alpha to be 5%
alpha <- 0.05
```

Fit an ANOVA model and use R’s TukeyHSD function to carryout the
multiple comparison procedure

``` r
aov_Fert <- aov(V2 ~ V1, data = df)
(tuk_Res <- TukeyHSD(aov_Fert))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = V2 ~ V1, data = df)
    ## 
    ## $V1
    ##                 diff        lwr         upr     p adj
    ## F1-Control  7.600000  4.7770648 10.42293521 0.0000016
    ## F2-Control  4.866667  2.0437315  7.68960188 0.0005509
    ## F3-Control  8.200000  5.3770648 11.02293521 0.0000005
    ## F2-F1      -2.733333 -5.5562685  0.08960188 0.0598655
    ## F3-F1       0.600000 -2.2229352  3.42293521 0.9324380
    ## F3-F2       3.333333  0.5103981  6.15626854 0.0171033

Extract the mean differences from tuk\_Res and create a new data frame
called pairs\_df

``` r
pairs_df <- as.data.frame(sapply(tuk_Res, function(x) x[,1]))
```

Extract MSE and error DF from ANOVA table, and compute the standard
error for the difference between two treatment means. From Kutner et
al. (17.30b), the standard error for pairwise comparison
*D* = *μ*<sub>*i*</sub> − *μ*<sub>*i*</sub>′ is
$s\[\\hat{D}\]=\\sqrt{MSE(1/n\_i+1/n\_i')}$. The number of observations
is equal in each group, so this could be simplified to
$s\[\\hat{D}\]=\\sqrt{2MSE/n}$.

``` r
summary(aov_Fert)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## V1           3 251.44   83.81   27.46 2.71e-07 ***
    ## Residuals   20  61.03    3.05                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
myMSE <- summary(aov_Fert)[1][[1]][[3]][[2]]
errorDF <- summary(aov_Fert)[1][[1]][[1]][[2]]
se_pair <- sqrt(myMSE*((1/NumofRepl)+(1/NumofRepl)))
```

Comparison of Multiple Comparison Methods
=========================================

All of them produce simultaneous confidence intervals of the general
form: sample statistic +/- margin of error, where margin of error =
multiplier \* SE(estimate). All of the methods differ in the form of the
multiplier

Fisher’s Least Significant Difference (LSD)
-------------------------------------------

``` r
LSD_crit <- qt(1-alpha/2, df = errorDF)
MofE <- se_pair*LSD_crit
pairs_df$lowerLSD <- pairs_df$V1-MofE
pairs_df$upperLSD <- pairs_df$V1+MofE
```

Tukey
-----

``` r
qtukey_val <- qtukey(1-alpha, NumOfGroups, df = errorDF)/sqrt(2)
MofE <- se_pair*qtukey_val
pairs_df$lowerTuk <- pairs_df$V1-MofE
pairs_df$upperTuk <- pairs_df$V1+MofE
```

P-values match the TukeyHSD output (right-tail probabilities) Note, the
test statistic is $q^\*=\\sqrt{2}\\hat{D}/s\[\\hat{D}\]$ (Kunter et al.,
17.32)

``` r
ptukey(abs(pairs_df$V1*sqrt(2)/se_pair), NumOfGroups, df=errorDF,lower.tail = F)
```

    ## [1] 1.637988e-06 5.509424e-04 5.148374e-07 5.986551e-02 9.324380e-01
    ## [6] 1.710330e-02

Bonferroni
----------

``` r
Bon_crit <- qt(1-alpha/2/NumOfCompar, df = errorDF)
MofE <- se_pair*Bon_crit
pairs_df$lowerBon <- pairs_df$V1-MofE
pairs_df$upperBon <- pairs_df$V1+MofE
```

Scheffe
-------

``` r
S_crit <- sqrt((NumOfGroups-1)*qf(1-alpha, df1=NumOfGroups-1, df2=errorDF))
MofE <- se_pair*S_crit
pairs_df$lowerS <- pairs_df$V1-MofE
pairs_df$upperS <- pairs_df$V1+MofE
```

Dunnett
-------

``` r
Dun_crit <- qNCDun(p=1-alpha, nu=errorDF, rho=(rep(0.5,times=NumOfGroups-1)), delta=rep(0,times=NumOfGroups-1), two.sided=T)
MofE <- se_pair*Dun_crit
pairs_df$lowerDun <- pairs_df$V1-MofE
pairs_df$upperDun <- pairs_df$V1+MofE
pairs_df$lowerDun[4:6] <- NA
pairs_df$upperDun[4:6] <- NA
```

Comparison of the 95% simultaneous confidence intervals produced

``` r
knitr::kable(pairs_df, digits = 3)
```

|            |      V1|  lowerLSD|  upperLSD|  lowerTuk|  upperTuk|  lowerBon|  upperBon|  lowerS|  upperS|  lowerDun|  upperDun|
|------------|-------:|---------:|---------:|---------:|---------:|---------:|---------:|-------:|-------:|---------:|---------:|
| F1-Control |   7.600|     5.496|     9.704|     4.777|    10.423|     4.648|    10.552|   4.525|  10.675|     5.038|    10.162|
| F2-Control |   4.867|     2.763|     6.971|     2.044|     7.690|     1.914|     7.819|   1.792|   7.942|     2.304|     7.429|
| F3-Control |   8.200|     6.096|    10.304|     5.377|    11.023|     5.248|    11.152|   5.125|  11.275|     5.638|    10.762|
| F2-F1      |  -2.733|    -4.837|    -0.629|    -5.556|     0.090|    -5.686|     0.219|  -5.808|   0.342|        NA|        NA|
| F3-F1      |   0.600|    -1.504|     2.704|    -2.223|     3.423|    -2.352|     3.552|  -2.475|   3.675|        NA|        NA|
| F3-F2      |   3.333|     1.229|     5.437|     0.510|     6.156|     0.381|     6.286|   0.258|   6.408|        NA|        NA|
