# VoidFinder Toolkit
Repositorio para el proyecto de Void Finder de la materia de Desarrollo de Software para computo Cientifico

## Descripcion
El presente trabajo consiste en un algoritmo que tiene como proposito tomar puntos de catalogos sinteticos de galaxias para encontrar regiones de baja densidad (Cosmic Voids).

Para ello se plantea implementar una pipeline de tres pasos , modulariazada y flexible que consite en un **input**, area de **procesamiento** y **output**.

**INPUT**


En principio como input se utiliza un dataset sintetico con las siguientes caracteristicas:

- Posicion (En 3D)
- Velocidad (En 3D)
- Masa # No es indispensable
- Densidad Maxima de los Voids

De momento se plantea utilizar como input datasets correspondientes a catalogos sinteticos de galaxias dado que los mismos cuentan con las caracteristicas anteriores y no asi los catalogos observacionales. 

Concretamente cada punto en el dataset se corresponde con una galaxia.

Es posible que se agregue una seccion de preprocesamiento para el input dependiendo del catalogo sintetetico a utilizar.

**PROCESAMIENTO**

Para el procesamiento se utiliza como base un algoritmo escrito en C/C++ y que hay que refactorizar. Internamente el mismo consta de 6 partes (Ver diagrama de componentes), las cuales hay que modulariazar en primera instancia. Elementos como **Redshift-Space Distortions** y **Geometrical Distortions** son pertinentes a catalogos observacionales con los que no se trabajara aqui, sin embargo se mantendran dichos modulos para futuras aplicaciones. El modulo de **Void Velocity** es secundario.

Se buscara agregar alternativas a los otros modulos para darle flexibilidad al algoritmo segun los requerimientos del usuario.

**OUTPUT**

El output principalmente consiste de dos elementos que son cruciales para caracterizar a un *Void* como estructura, el **centro** del *void* y el **radio** del mismo.

Se utilizaran criterios distintos para definir estos parametros de salida en la seccion de **PROCESAMIENTO** como por ejemplo el numero de trazadores existentes a un determinado radio.

Los datos de salidas se guardaran en archivos, con extensiones diversas y habra posibilidad de graficar la salida en 3D.



