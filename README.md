# Makro

## To-do list:

- Finish tasks that are not struck through
- Create a Latex file
- Upload tentative versions to git

## Part I: Writing complex functions in R

### Task 1: Functions to estimate regression models using the SSVS prior

Read the following introduction to functions in R [(Advanced R by Hadley Wickham)](http://adv-r.had.co.nz/Functions.html). In this task, you have to write a function that performs Bayesian inference in a regression model with a stochastic search variable selection prior based on the code discussed in class (see code SSVS.R). Use the economic growth dataset of Fernandez, Ley, and Steel (2001, *J. Applied Econometrics*) provided in the BMS package in R. To get this data, type *data(datafls)* after loading the BMS package.

- [x] ~~Write a function that takes the explanatory variables ***X*** as well as the endoge- nous variable ***y*** as input. In the growth dataset, the first column contains the endogenous variable whereas the remaining columns are the explanatory variables. Think about what additional inputs might be helpful! (Hint: you might want to vary *nsave* and *nburn*.) Also think carefully about the potential output of the function! (Hint: R functions can only return a single object, so use a list object.)~~

- [ ] Run the function using different values for *τ0* and *τ1*. What happens to the posterior inclusion probabilities (PIP.mean) if *τ0* is set equal to *1e−15*? Describe this finding verbally and graphically!

- [ ] The variables in ***X*** all feature a different scale. This causes problems since the simple implementation of the code sets *τ0* and *τ1* equal to fixed values that are independent of the scaling of the data. Try to standardize the data such that all columns of ***X*** (and ***y***) have mean zero and variance one.

- [ ] *(ADVANCED)* Try to implement the semi-automatic approach of George, Sun, and Ni (2008, *J. Econometrics*) in your SSVS function. This amounts to first estimating the OLS standard deviations and then scaling *τ0* and *τ1* using the corresponding OLS standard deviations.

## Part II: Model uncertainty in economic growth regressions

### Task 1: Bayesian Normal Linear Regression

Read the paper by Fernandez, Ley and Steel (2001, *J. Applied Econometrics*) as well as Chapter 11 in Gary Koop’s textbook. Consider the data for the paper by Fernandez, Ley and Steel (2001, *J. Applied Econometrics*), available from the BMS package

- [ ] Reproduce the results in Table 11.1 (in Koop) using the BMS package in R.

- [ ] Use your custom function for the SSVS prior to reproduce Table 11.1.

- [ ] How do results differ? To what extent is this related to the specific choices of *τ0* and *τ1*?
