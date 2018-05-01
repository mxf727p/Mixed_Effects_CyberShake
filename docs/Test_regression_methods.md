# **Mixed Effects Analysis on PGME Residuals**
### _by Xiaofeng Meng_
### _last updated on 05/01/2018_

**1. Input files**  
Ground motion data (either predicted or observed) or gournd motion residual data can be used as input file. The required format of input file is csv.  

```R
regression_flag <- 0
if (regression_flag) {
  regression_model <- 4   # 1: nlme; 2: nls; 3: nlsLM; 4: nlmer
}
L2L_flag <- 0
S2S_flag <- 0
P2P_flag <- 0
``` 
 
