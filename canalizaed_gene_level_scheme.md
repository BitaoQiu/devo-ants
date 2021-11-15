

1. Between-caste deviation<sub>age</sub> = <img src="http://www.sciweavers.org/tex2img.php?eq=%5Cfrac%7BMean_%7Bgyne%2Cage%7D%20-%20Mean_%7Bworker%2Cage%7D%7D%7B%5Csqrt%7BVar_%7Bgyne%2Cage%7D%20%2B%20Var_%7Bworker%2Cage%7D%7D%7D&bc=White&fc=Black&im=png&fs=12&ff=modern&edit=0" align="center" border="0" alt="\frac{Mean_{gyne,age} - Mean_{worker,age}}{\sqrt{Var_{gyne,age} + Var_{worker,age}}}" width="215" height="43" />

** So that gene with high between-caste expression difference or with low within-caste variation are with high score.**

2. C.trend = -log(Spearman.Cor(|Between-caste deviation<sub>age</sub>|, age)$p.value) # Trend of canalization.

** Calculate the correlation of age and the absolute value of between-caste deviation. Genes with increasing divergence are with high C.trend.**

** Set C.trend = 0 if Between-caste deviation<sub>pupa.young</sub> * Between-caste deviation<sub>pupa.old</sub> < 0, so that genes with opposite expression bias at late stages are with 0 C.trend. **
3. Canalization score = Between-caste deviation<sub>pupa.old</sub> * C.trend
