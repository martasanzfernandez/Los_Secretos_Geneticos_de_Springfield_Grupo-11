# ANÁLISIS DE EXPRESIÓN DIFERENCIAL

Una vez se ha realizado el análisis de expresión diferencial y se han obtenido los genes expresados diferencialmente entre pacientes con la condicion "*obesos 2*" y "*normopeso*", el último paso es relacionar estos genes con su función molecular y biológica. Por tanto, a partir de los resultados del análisis de expresió diferencial mediante *DESeq2*, se realizó un análisis de enriquecimiento.

Para ello, se utilizó el método **Gene Ontology (GO)** para **procesos biológicos (BP)**, el cual clasifica los genes en función del proceso celular en el que intervengan. Aquellos procesos con más de un gen implicado serán aquellos más representados y relacionados con el análsis. Se visualizaron los resultado mediante un *dotplot*, en el que se observa que la función *feeding behaviour* es la más representada, por un total de 3 genes.

También se realizó un segundo análisis de enriquecimiento de **Reactome ORA**, en el cual se clasifican los genes en función de las vías biológicas y rutas moleculares específicas, mostrando cómo interactúan dentro de procesos celulares concretos. De nuevo se representaron los resultados mediante un *dotplot*. 

Para ambos análisis, la gran mayoría de los procesos biológicos y rutas metabólicas están relacionadas directamente con la obesidad y la acumulación de grasa en células adiposas, lo que evidencia en mayor medida las diferencias en la expresión génica observadas entre los dos grupos del estudio.
