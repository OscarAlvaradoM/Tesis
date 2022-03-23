# Tesis

Repositorio con los archivos de código, documentación, imágenes y benchmarking, de la versión preeliminar de mi tesis de Licienciatura en Física de Óscar Anuar Alvarado Morán. 

## Planteamiento
El problema abordado en esta tesis es el del flujo de calor a través de un cubo,que puede ser modelado matemáticamente mediante la ecuación de Poisson o de Laplace cuando no hay convección, mientras que cuando sí la hay, simplemente es una ecuación estacionaria de difusión-convección (o advección). Este problema se encuentra inmerso en muchas situaciones de la física, tal como la electrostática, el problema del flujo irrotacional y en la gravitación. 

La metodología computacional que se utilizó para resolver numéricamente este modelo matemático es el que lleva por nombre "Método del Volumen Finito" (MVF), que se trata de un método integral para la solución de sistemas de ecuaciones diferenciales. Éste se basa, primero que nada, en discretizar el espacio e integrar sobre cada volumen de dicho mallado. La idea de deshacerse de una derivada con el teorema fundamental del cálculo y evaluar en la frontera. Lo que se evaluará será en forma del método de las diferencias finitas, tomando a las derivadas de cada dimensión como unos diferenciales de una cantidad. Luego se aproximan ciertos coeficientes a partir del coeficiente de difusión, lo que resulta en una ecuación para cada volumen del mallado, lo que lleva a tener un sistema de ecuaciones lineáles algebráico para todo el mallado. El sistema de ecuaiones lineales depende  en gran parte de las condiciones de frontera del problema físico. 

Una vez formado el sistema de ecuaciones lineales se resuelve mediante distintos métodos, empezando por los solvers por defecto de python y julia como los algoritmos de Jacobi, Gauss-Seidel, SOR, GMRES, GMRES reiniciado, GMRES precondicionado (con precondicionadores de Jacobi, Gauss-Seidel y SOR) y GMRES precondicionado reiniciado (con los mismo precondicionadores mencionados anteriormente. Se mide el tiempo de ejecución de cada algoritmo y se comparan entre sí.


## Organización del repositorio
Se organiza en 4 carpetas principales:

- **Src**: Contiene el código usado para la resolución del problema físico planteado por la ecuación de Laplace y de Poisson, así como los *solvers* utilizados para la comparación entre algoritmos secuenciales y paralelos.
- **Imgs**: Imágenes y diagramas usados en la tesis u obtenidos de ésta, como los resultados numéricos, analíticos, la comparación entre estos y la comparación de tiempos entre los distintos solvers.
- **Documentación**: Contiene la documentación de la tesis, en especial el mismo documento de tesis.
- **Benchmarking**: Contiene los archivos de medición de tiempos de las etapas del MVF y de los solvers implementados. Es importante recalcar que contiene las mediciones de los 10 tiempos por etapa y por solver para poder hacer estadística.

**Investigación realizada gracias al programa UNAM-PAPIIT IA-104720**
