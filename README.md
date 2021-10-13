# CTL
A matlab implementation of paper "Credal Transfer Learning with Multi-estimation for Missing Data."

Published in IEEE ACCESS '2020’

Link: https://ieeexplore.ieee.org/abstract/document/9046814

**In case the repository or the publication was helpful in your work, please use the following to cite the original paper,**

 <pre><code> @ARTICLE{9046814,
  author={Ma, Zongfang and Liu, Zhe and Zhang, Yiru and Song, Lin and He, Jihuan},
  journal={IEEE Access}, 
  title={Credal Transfer Learning With Multi-Estimation for Missing Data}, 
  year={2020},
  volume={8},
  pages={70316-70328},
  doi={10.1109/ACCESS.2020.2983319}}
 </code></pre> 
Abstract:
Transfer learning (TL) has grown popular in recent years. It is effective to improve the classification accuracy in the target domain by using the training knowledge in the related domain (called source domain). However, the classification of missing data (or incomplete data) is a challenging task for TL because different strategies of imputation may have strong impacts on learning models. To address this problem, we propose credal transfer learning (CTL) with multi-estimation for missing data based on belief function theory by introducing uncertainty and imprecision in data imputation procedure. CTL mainly consists of three steps: Firstly, the query patterns are reasonably mapped into multiple versions in source domain to characterize the uncertainty caused by missing values. Afterwards, the multiple mapping patterns are classified in the source domain to obtain the corresponding outputs with different discounting factors. Finally, the discounted outputs, represented by the basic belief assignments (BBAs), are submitted to a new belief-based fusion system to get the final classification result for the query patterns. Three comparative experiments are given to illustrate the interests and potentials of CTL method.
