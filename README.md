# gwas_filter

<font size=8>**The *gwas_filter* is an R script that can efficiently and accurately filter genome-wide association studies (GWASs) from the GWAS Catalog Website. The selection principles of GWASs were established based on previous studies. The process of manual filtering in the GWAS Catalog was abstracted as standard algorithms. This R script was written by two programmers and tested many times. It takes six steps for gwasfilter.R to filter GWASs. There are five main self-defined functions among this R script. GWASs can be filtered based on "whether the GWAS has been replicated" "sample size" "ethnicity of the study population" and other conditions. It takes no more than 1 second for this script to filter GWASs of a single trait. This R script is user-friendly and provides an efficient and standard process to filter GWASs flexibly.

If you find gwasfilter is useful, please cite our paper: Songchun Yang, Chongyang Li, Yizhen Hu, Qiufen Sun, Jianqiao Pan, Dianjianyi Sun, Baoshan Ma, Jun Lyu, Liming Li. gwasfilter: an R script to filter genome-wide association study. Chin J Epidemiol,2021,42(10):1876-1881. DOI：10.3760/cma.j.cn112338-20200731-01003.
**

  中文介绍：
gwasfilter.R是一套能高效准确地从GWAS Catalog公开数据库中筛选全基因组关联研究（GWAS）的R脚本。参考既往研究制定GWAS的筛选原则。将人工在GWAS Catalog的筛选过程抽象为标准的算法，两名程序员共同撰写R脚本（gwasfilter.R）后，由他人多次对脚本进行测试。采用gwasfilter.R筛选GWAS包含6个步骤。该脚本内置5个主要函数，可以同时根据“是否有验证人群”、“样本量大小”和“研究人群种族”等条件对GWAS进行筛选。筛选单个性状时，程序用时不超过1秒。gwasfilter.R操作简便，筛选过程高效而标准化，可灵活应用于GWAS筛选。
  
如果您觉得该脚本对您的研究有帮助，请引用我们的论文：杨淞淳,李重阳,胡一祯,孙秋芬,潘建桥,孙点剑一,马宝山,吕筠,李立明.gwasfilter：用于筛选全基因组关联研究的R脚本[J].中华流行病学杂志,2021,42(10):1876-1881. DOI：10.3760/cma.j.cn112338-20200731-01003.

  

