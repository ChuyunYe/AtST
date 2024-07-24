# AtST

This is an initial version of the 'AtST' package. The goal of 'AtST' is to provide test statistics of valid mediators, based on linear structural equations, and the corresponding sequential tests with controlled False Discovery Rate (FDR). 

The statistical method of this package is based on the paper 'Adaptive-to-sub-null testing for mediation effects in structural equation models' (authors: Jiaqi Huang, Chuyun Ye, Lixing Zhu). Users can refer to this paper for more theoretical and technique details. 

## Installation

Users can install this version by running

```{r}
# install.packages('devtools')
library(devtools)
devtools::install_github("ChuyunYe/AtST")
```

## Elements in the Package

### data_list_example

Data (list). This data list is consisted of 4 elements, including 'X', 'A', 'M', 'Y':

- Element 'X': a matrix of covariate. It is arranged as a sample matrix of covariate variables. The covariate data of the i-th individual are saved in the i-th row, while the data of the j-th covariate are saved in the j-th column. 
- Element 'A': a vector of exposure variable. The exposure variable of the i-th individual is saved in the i-th entry.
- Element 'M': a matrix of mediators. It is arranged as a sample matrix of mediators. The mediator data of the i-th individual are saved in the i-th row, while the data of the j-th mediator are saved in the j-th column. 
- Element 'Y': a vector of outcome variable. The outcome variable of the i-th individual is saved in the i-th entry.

This data is provided as an example of the input data for function 'AtST_stat'. Users who are interested in using 'AtST_stat' shall arrange the data as 'data_list_example' does.

Users can can this data list by running:

```{r}
data(data_list_example)
```

### AtST_stat

Function. This function plays roles in test statistics calculation. 

With a list of data input, 'AtST_stat' provides test statistics of each mediator. The output of 'AtST_stat' is a list of 3 elements, including the AtST statistics, the T statistics of each $\beta$, and the T statistics of each $\gamma$. For more details, please refer to 'Adaptive-to-sub-null testing for mediation effects in structural equation models'. 

Moreover, the users can self-modify the tuning parameters 'an' and 'bn' in this function. If no specified 'an' or 'bn' is provided as arguments, the default value of these tuning parameters are similar with those in the paper.

### SeqHyTest_FDR

Function. This function plays roles in sequential testing by controlling FDR.

With a vector of p-value and expected FDR input, 'SeqHyTest_FDR' provides rejection cutoff and rejection result. The output of 'SeqHyTest_FDR' is a list of 2 elements, including the scalar cutoff value, and the vector of rejection indicators. Note that the entrees of this vector is one-to-one mapped to mediators. In other words, that the j-th entry of this vector is equal to 1, is equivalent to that the j-th mediator is rejected.

Moreover, this function is capable to perform two FDR controls procedures. Users can choose different procedures by setting "tuning_type='DOS'" or "tuning_type='Storey'". And for each procedure, users shall give a pre-specified tuning parameter 'tuning' as an argument. For illustration simplicity, we do not dive deep into the details of how to choose 'tuning_type' and 'tuning'. Users can refer to the paper. 

## Example

```{r}
# extracting data
data(data_list_example)
# calculating test statistics
TStat = AtST_stat(data.list.example)$TestStat
# calculating p-values
p_value.seq = 1-pchisq(q=TStat,df=1)
# sequential testing
SeqHyTest_FDR(p_value.seq=p_value.seq,tuning_type='DOS',tuning=1,FDR=0.2)
```





