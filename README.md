**MTest**: A Procedure for Multicollinearity Testing using Bootstrap.

## About

Functions for detecting multicollinearity. This test gives statistical support to two of the most famous methods for detecting multicollinearity in applied work: Klein’s rule and Variance Inflation Factor (VIF). See the URL for the paper associated with this package, Morales-Oñate and Morales-Oñate (2023) [doi:10.33333/rp.vol51n2.05](doi:10.33333/rp.vol51n2.05)


Consider the regression model

$$
Y_i	= \beta_{0}X_{0i} + \beta_{1}X_{1i} + \cdots+ \beta_{p}X_{pi} +u_i
$$

where $i = 1,\ldots,n$, $X_{j,i}$ are the predictors with $j = 1,\ldots,p$, $X_0 = 1$ for all $i$ and $u_i$ is the gaussian error term. 

In order to describe Klein's rule and VIF methods, we need to define *auxiliary regressions* associated to the abobe model (global). An example of an auxiliary regressions is:

$$
X_{2i} =  \gamma_{1}X_{1i} + \gamma_{3}X_{3i} + \cdots+ \gamma_{p}X_{pi} +u_i.
$$

In general, there are $p$ auxiliary regressions and the dependent variable is omitted in each auxiliary regression. Let $R_{g}^{2}$ be the coefficient of determination of the global model and $R_{j}^{2}$ the $j\text{th}$ coefficient of determination of the $j\text{th}$ auxiliary regression.


## MTest

Given a regression model, Mtest is based on computing estimates of $R_{g}^{2}$ and $R_{j}^{2}$ from $n_{b}$ bootstrap samples obtained from the dataset, $R_{g_{b}}^{2}$ and $R_{j_{b}}^{2}$ respectively. 

Therefore, in the context of MTest, the VIF rule translates into:

$$
H_0:\mu_{R_{j_{b}}^{2}}\geq 0.90,
$$
and 

$$
H_a:\mu_{R_{j_{b}}^{2}}<0.90.
$$

We seek an achieved significance level (ASL)

$$
\text{ASL} = \text{Prob}_{H_0}\{\mu_{R_{j_{b}}^{2}}\geq 0.90\}
$$

estimated by 

$$
\widehat{\text{ASL}}_{n_{b}} = \text{Card}(\mu_{R_{j_{b}}^{2}}\geq 0.90) /{n_{b}}
$$

In a similar manner, the  Klein's rule translates into:

$$
H_0:\mu_{R_{j_{b}}^{2}}\geq \mu_{R_{g_{b}}^{2}},
$$
and 

$$
H_a:\mu_{R_{j_{b}}^{2}}<\mu_{R_{g_{b}}^{2}}.
$$

We seek an achieved significance level

$$
\text{ASL} = \text{Prob}_{H_0}\{\mu_{R_{j_{b}}^{2}}\geq \mu_{R_{g_{b}}^{2}}\}
$$

estimated by 

$$
\widehat{\text{ASL}}_{n_{b}} = \text{Card}(\mu_{R_{j_{b}}^{2}}\geq\mu_{R_{g_{b}}^{2}})/{n_{b}}.
$$

It should be noted that this set up let us formulate VIF and Klein's rules in terms of statistical hypothesis testing. 




## MTest: the algorithm


$R_{g_{b}}^{2}$ and $R_{j_{b}}^{2}$ are the distributions of $R_{g}^{2}$ and $R_{j}^{2}$ induced by applying the bootstrap procedure to the dataset. Achieved significance level is computed for the VIF and Klein's rule. In the following we describe the procedure step by step:


- Create $n_{b}$ samples from original data with replacement of a given size ($n_{sam}$). 
- Compute $R_{g_{b}}^{2}$ and $R_{j_{b}}^{2}$ from each $n_{b}$ samples. This outputs a $B_{n_{b}\times (p+1)}$ matrix.
- Compute $\widehat{\text{ASL}}_{n_{b}}$ for the VIF and Klein's rule.


Note that the matrix $B_{n_{b}\times (p+1)}$ allow us to inspect results in detail and make further tests such as boxplots, pariwise Kolmogorov-Smirnov (KS) of the predictors and so on.