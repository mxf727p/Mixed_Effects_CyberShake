# **Test different regression methods in R**
### _by Xiaofeng Meng_
### _last updated on 05/01/2018_

The GMPE is from Villani and Abrahamson, BSSA, 2015:   
_ln_(SA) = b_1 + b_2 * M + b_3 * M^2 + (b_4 + b_5 * M) * log(sqrt(b_R^2 + R_rup^2)) + b_6 * log(Vs_30 / 500)  

The test dataset if NGA_W2 in southern California:   
<img src="https://github.com/mxf727p/Mixed_Effects_CyberShake/blob/master/docs/figures/NGA_W2_SC_map.png" width="500" height="500" />
<figcaption>Fig.1 Map of southern California. Red circles are 78 earthquakes. Black triangles are 125 stations.</figcaption>


