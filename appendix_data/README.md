# appendix_data
Data associated with an example subset of scenarios in the appendix of Ogle _et. al. (in prep)_

## Files:


#### randomDataAppendix.csv
Randomly generated data used in the Ogle _et. al._ simulation experiment 

|  Column  |  Description  | 
|  ------  |  -----------  | 
|  Y  | randomly generated data given the “true” parameters (simulation factors) described in scenarioDataAppendix.csv; these data are used to fit the hierarchical Bayesian model. |
|  Y.cp  | Same data as in Y, but these data are used to fit the complete pooling Bayesian model.  |
|  Y.nh  | Same data as in Y, but these data are used to fit the non-hierarchical Bayesian model.  |
|  scenarioid  | unique identifier for the scenario _(d)_ associated with each pseudo dataset  |
| groupid  | unique identifer for the group-level _(k)_ associated with each observation within each pseudo dataset   |

#### scenarioDataAppendix.csv
Description of scenario settings for comparisons in the Ogle _et. al._ simulation experiment 

|  Column  |  Description  | 
|  ------  |  -----------  | 
|  scenarioid  |  unique identifier for the scenario _(d)_ associated with each pseudo dataset  |
|  K  |group size _(K)_ ; i.e., the number of group levels represented in each pseudo dataset  |
|  A  |  number of groups with the same mean  |
|  s  | standard deviation _(s)_ describing among group variation   |
|  max.nk  |  maximum sample size per group, max(nk)  |
|  num.unbal  | number of group levels with “missing data” (i.e., with nk = min(nk))  |
| dist | Distribution from with pseudo data are drawn (0 = uniform, 1 = normal) |
|  min.nk  |  minimum sample size per group, min(nk) |
 

#### groupComparisionDataAppendix.csv
True means for group comparison and unique IDs for groups and scenarios. 

|  Column  |  Description  | 
|  ------  |  -----------  | 
|  scenarioid  |   unique identifier for the scenario _(d)_ associated with each pseudo dataset  |
|  groupid1  |  unique identifier for first group in comparison of means in each scenario  |
|  groupid2  |  unique identifier for second group in comparison of means in each scenario  |
|  truth1  |  true mean for first group in comparison of means in each scenario  |
|  truth2  |  true mean for second group in comparison of means in each scenario  |
|  true.diff  |  true differences between means in each group comparison in each scenario  |