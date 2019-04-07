# MiNoM
This toolbox is mainly for MIx-NOrM-based (MiNoM) scan matching algorithm for autonomous driving. The algorithm is based on [1], and can be viewed as a generalized sparse ICP [2] by incorporating more robust component into the objective function. Some major revisions are done which include point-to-plane metric, fine-tuned EM algorithm, and alternating directions method of multipliers (ADMM) [2] in the optimizer.  

The proposed MiNoM, especially the point-to-plane MiNoM runs fast with Matlab R2018Ra, and can be real-time if you implement it with C/C++. The speed improvement is expected to be doubled when Point Cloud Libarary(PCL-1.8.0) + Eigen-3.0 + Visual Studio 2015 are recommended. 

PS: It is worth noting that the source code in [2] is problematic when you want to fairly compare other scan matching algorithms with sparse ICP. The functions of point2point() and point2plane() within namespace "RigidEstimator" are procedurally wrong. It is highly encouraged to re-implement or modify the "RigidEstimator".   

[1] Wang, D., Xue, J., Tao, Z., Zhong, Y., Cui, D., Du, S. and Zheng, N., 2018, October. Accurate Mix-Norm-Based Scan Matching. In 2018 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS) (pp. 1665-1671). IEEE.

[2] Mavridis, P., Andreadis, A. and Papaioannou, G., 2015. Efficient sparse icp. Computer Aided Geometric Design, 35, pp.16-26.
