**MTest**: A Procedure for Multicollinearity Testing using Bootstrap.

## About

Functions for detecting multicollinearity. This test gives statistical
support to two of the most famous methods for detecting
multicollinearity in applied work: Klein’s rule and Variance Inflation
Factor (VIF). See the URL for the paper associated with this package,
Morales-Oñate and Morales-Oñate (2023) <doi:10.33333/rp.vol51n2.05>

Consider the regression model

*Y*<sub>*i*</sub> = *β*<sub>0</sub>*X*<sub>0*i*</sub> + *β*<sub>1</sub>*X*<sub>1*i*</sub> + ⋯ + *β*<sub>*p*</sub>*X*<sub>*p**i*</sub> + *u*<sub>*i*</sub>

where *i* = 1, …, *n*, *X*<sub>*j*, *i*</sub> are the predictors with
*j* = 1, …, *p*, *X*<sub>0</sub> = 1 for all *i* and *u*<sub>*i*</sub>
is the gaussian error term.

In order to describe Klein’s rule and VIF methods, we need to define
*auxiliary regressions* associated to the abobe model (global). An
example of an auxiliary regressions is:

*X*<sub>2*i*</sub> = *γ*<sub>1</sub>*X*<sub>1*i*</sub> + *γ*<sub>3</sub>*X*<sub>3*i*</sub> + ⋯ + *γ*<sub>*p*</sub>*X*<sub>*p**i*</sub> + *u*<sub>*i*</sub>.

In general, there are *p* auxiliary regressions and the dependent
variable is omitted in each auxiliary regression. Let
*R*<sub>*g*</sub><sup>2</sup> be the coefficient of determination of the
global model and *R*<sub>*j*</sub><sup>2</sup> the *j*th coefficient of
determination of the *j*th auxiliary regression.

## MTest

Given a regression model, Mtest is based on computing estimates of
*R*<sub>*g*</sub><sup>2</sup> and *R*<sub>*j*</sub><sup>2</sup> from
*n*<sub>*b*</sub> bootstrap samples obtained from the dataset,
*R*<sub>*g*<sub>*b*</sub></sub><sup>2</sup> and
*R*<sub>*j*<sub>*b*</sub></sub><sup>2</sup> respectively.

Therefore, in the context of MTest, the VIF rule translates into:

*H*<sub>0</sub> : *μ*<sub>*R*<sub>*j*<sub>*b*</sub></sub><sup>2</sup></sub> ≥ 0.90,
and

*H*<sub>*a*</sub> : *μ*<sub>*R*<sub>*j*<sub>*b*</sub></sub><sup>2</sup></sub> &lt; 0.90.

We seek an achieved significance level (ASL)

ASL = Prob<sub>*H*<sub>0</sub></sub>{*μ*<sub>*R*<sub>*j*<sub>*b*</sub></sub><sup>2</sup></sub> ≥ 0.90}

estimated by

$$
\widehat{\text{ASL}}\_{n\_{b}} = \\(\mu\_{R\_{j\_{b}}^{2}}\geq 0.90) /{n\_{b}}
$$

In a similar manner, the Klein’s rule translates into:

*H*<sub>0</sub> : *μ*<sub>*R*<sub>*j*<sub>*b*</sub></sub><sup>2</sup></sub> ≥ *μ*<sub>*R*<sub>*g*<sub>*b*</sub></sub><sup>2</sup></sub>,
and

*H*<sub>*a*</sub> : *μ*<sub>*R*<sub>*j*<sub>*b*</sub></sub><sup>2</sup></sub> &lt; *μ*<sub>*R*<sub>*g*<sub>*b*</sub></sub><sup>2</sup></sub>.

We seek an achieved significance level

ASL = Prob<sub>*H*<sub>0</sub></sub>{*μ*<sub>*R*<sub>*j*<sub>*b*</sub></sub><sup>2</sup></sub> ≥ *μ*<sub>*R*<sub>*g*<sub>*b*</sub></sub><sup>2</sup></sub>}

estimated by

$$
\widehat{\text{ASL}}\_{n\_{b}} = \\\\\mu\_{R\_{j\_{b}}^{2}}\geq\mu\_{R\_{g\_{b}}^{2}}\\/{n\_{b}}.
$$

It should be noted that this set up let us formulate VIF and Klein’s
rules in terms of statistical hypothesis testing.

## MTest: the algorithm

*R*<sub>*g*<sub>*b*</sub></sub><sup>2</sup> and
*R*<sub>*j*<sub>*b*</sub></sub><sup>2</sup> are the distributions of
*R*<sub>*g*</sub><sup>2</sup> and *R*<sub>*j*</sub><sup>2</sup> induced
by applying the bootstrap procedure to the dataset. Achieved
significance level is computed for the VIF and Klein’s rule. In the
following we describe the procedure step by step:

-   Create *n*<sub>*b*</sub> samples from original data with replacement
    of a given size (*n*<sub>*s**a**m*</sub>).
-   Compute *R*<sub>*g*<sub>*b*</sub></sub><sup>2</sup> and
    *R*<sub>*j*<sub>*b*</sub></sub><sup>2</sup> from each
    *n*<sub>*b*</sub> samples. This outputs a
    *B*<sub>*n*<sub>*b*</sub> × (*p*+1)</sub> matrix.
-   Compute $\widehat{\text{ASL}}\_{n\_{b}}$ for the VIF and Klein’s
    rule.

Note that the matrix *B*<sub>*n*<sub>*b*</sub> × (*p*+1)</sub> allow us
to inspect results in detail and make further tests such as boxplots,
pariwise Kolmogorov-Smirnov (KS) of the predictors and so on.
