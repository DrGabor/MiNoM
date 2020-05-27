# MiNoM
This toolbox is mainly for MIx-NOrM-based (MiNoM) scan matching algorithm in robotics. The algorithm is based on [1], and can be viewed as a boost to iterative closest points (ICP) by incorporating robust residual error modelling (REM) into the outlier rejection stage. Some major revisions are done which include point-to-plane metric and fine-tuned EM algorithm. 

The max_irls_iter in the original paper is a fixed parameter, and I found that this parameter should be gradually increased otherwise MiNoM will be pre-maturely trapped into local minimals. In other word, the on-line parameter learning (OPL) and transformation estimator (TE) should be less accurate in the early iterations in order to be robust when given the large initial transformations. 

The proposed MiNoM, especially the point-to-plane MiNoM runs fast with Matlab R2018Ra, and can be real-time if you implement it with C/C++. The speed is expected to be within 50~100ms when using Point Cloud Libarary(PCL-1.8.1) + Visual Studio 2015. 

If you want to see the fitting results on residual errors, please run testRED.m. 

PS: It is worth noting that the source code in [2] is problematic when you want to fairly compare other scan matching algorithms with sparse ICP. The functions of point2point() and point2plane() within namespace "RigidEstimator" are procedurally wrong. It is highly encouraged to re-implement or modify the "RigidEstimator".   

[1] Wang, D., Xue, J., Tao, Z., Zhong, Y., Cui, D., Du, S. and Zheng, N., 2018, October. Accurate Mix-Norm-Based Scan Matching. In 2018 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS) (pp. 1665-1671). IEEE.

[2] Bouaziz, Sofien, Andrea Tagliasacchi, and Mark Pauly. "Sparse iterative closest point." Computer graphics forum. Vol. 32. No. 5. Oxford, UK: Blackwell Publishing Ltd, 2013.
