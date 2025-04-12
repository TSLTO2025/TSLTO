# TSLTO
TSLTO (Tucker decomposition-based Sparse Low-Rank high-Order Tensor Optimization model) is a model for tensor imputation and anomaly diagnosis. Specific model and evaluations can be found in paper [Traffic Flow Data Completion and Anomaly Diagnosis via Sparse and Low-Rank Tensor Optimization](https://arxiv.org/abs/2504.02245 "our paper").
## Quik Run Guide
:star: Don't forget to load `Tensor_Toolbox` and `Tensorlab` first!
### About Synthetic Data
If you are interested in generating our synthetic data, you can run `TSLTO/synthetic_dataset.m` .
### About Real-World Data
Our model can also handle imputation and anomaly diagnosis tasks in real-world data (we only use [Guangzhou](https://zenodo.org/records/1205229 "You can get raw data Guangzhou here")).  
You first need to run
```matlab
load(guangzhou.mat);
```
then run `TSLTO/GUANGZHOU`.  
### Contact Us
Junxi Man `22271014@bjtu.edu.cn`  
Yumin Lin `21261047@bjtu.edu.cn`  
*For more information, feel free to ask!* 
