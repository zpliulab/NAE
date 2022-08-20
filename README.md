# NAE: Evaluating gene regulatory network activity from dynamic expression data by regularized constraint programming #


**To evaluate the activity of regulatory events in the form of network, we propose a network activity evaluation (NAE) framework by measuring the consistency between network architecture and gene expression data across specific states based on mathematical programming.**

**Please cite our paper (Early Access) if it helps you.**
  ` C. Wang, S. Xu and Z. -P. Liu, "Evaluating gene regulatory network activity from dynamic expression data by regularized constraint programming," in IEEE Journal of Biomedical and Health Informatics, 2022, doi: 10.1109/JBHI.2022.3199243. `


![workfolw](https://github.com/zpliulab/NAE/blob/master/workflow.jpg)

Computer Environment [![MATLAB](https://img.shields.io/badge/MATLAB-R2022a-green.svg "MATLAB")](https://ww2.mathworks.cn/products/matlab.html "MATLAB"):
-
- MATLAB R2022a (requires Parallel Computing Toolbox installed)
- If your working environment cannot enable parallel pooling, please replace *parfor* with *for*


File:
-
- code: subfunction used in the NAE process
- data: all datasets.
- demo.m: run the simulation example

The core interface function of NAE is code\refer_NAE.m

    [newgra,obj,opt] = refer_NAE(data,network,Lp,parameter);

where *data* represents the gene expression profile of the time series; *network* is the network to be evaluated; *Lp* represents the norm type used; *parameter* is a structure, where parameter.i is the regularization parameter that must be set. The return value *newgra* is the inference network, *obj* is the loss value generated by each iteration, and *opt* is the parameter used.

Please write to [zpliu@sdu.edu.cn](mailto:zpliu@sdu.edu.cn) if you have any questions.
