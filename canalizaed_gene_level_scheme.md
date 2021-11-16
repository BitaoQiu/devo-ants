1. Between-caste divergence<sub>age</sub> = <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\inline&space;\frac{Mean_{gyne,age}&space;-&space;Mean_{worker,age}}{\sqrt{Var_{gyne,age}&space;&plus;&space;Var_{worker,age}}}" title="\inline \frac{Mean_{gyne,age} - Mean_{worker,age}}{\sqrt{Var_{gyne,age} + Var_{worker,age}}}" />

>Genes with high between-caste expression difference or with low within-caste variation are with high between-caste divergence.

2. C.trend = -log(p.value of cor(|Between-caste divergence<sub>age</sub>|, age))

>Quantify the trend of between-caste divergence across developmental stages. Genes with increasing divergence are with high C.trend.

>Set C.trend = 0 if Between-caste divergence<sub>pupa.old</sub> * Between-caste divergence<sub>imago</sub> < 0, so that genes with opposite expression bias in late stages are with C.trend = 0.

3. Canalization score = Between-caste divergence<sub>pupa.old</sub> * C.trend

>Canalization score at gene level is the combination of developmental trend and terminal-stage (Old pupae) between-caste divergence.
