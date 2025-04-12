# TSLTO
TSLTO (Tucker decomposition-based Sparse Low-Rank high-Order Tensor Optimization model) is a model for tensor imputation and anomaly diagnosis. Specific model and evaluations can be found at paper [Traffic Flow Data Completion and Anomaly Diagnosis via Sparse and Low-Rank Tensor Optimization](https://arxiv.org/abs/2504.02245 "our paper").
## Quik Run Guide
### About Synthetic Data
If you are interested in the generation of our synthetic data, you can run `TSLTO/synthetic_dataset.m` .
### About Real-World Data
Our model can also deal with imputation and anomaly diagnosis tasks from real-world data (we only use [Guangzhou](https://zenodo.org/records/1205229 "You can get raw data Guangzhou here")).  
You first need to run
```matlab
load(guangzhou.mat);
```
then run `TSLTO/GUANGZHOU`.
