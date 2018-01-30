# Spatial graph embeddings and coupler curves

This program implements a method for obtaining edge lengths of a minimally rigid graph with many real spatial embeddings.
The method is based on sampling over two parameter family that preserves so called coupler curve.
See [project website](http://jan.legersky.cz/project/spatialgraphembeddings/) for the details.

Moreover, it includes Qt application for plotting coupler curves
of the 7-vertex minimally rigid graph with the maximal number of embeddings, G48.

The main functionality is provided by the package *graphEmbeddings3D*, see [Documentation](http://jan.legersky.cz/public_files/spatialGraphEmbeddings/documentation/).

## Requirements and installation
  * Python 2.7
  * For solving the system of equations corresponding to graph embeddings, polynomial homotopy continuation by the package `phcpy` is used
  ([homepages.math.uic.edu/~jan/phcpy_doc_html/](http://homepages.math.uic.edu/~jan/phcpy_doc_html/)).
  * In the sampling heuristic, clustering is done by `DBSCAN` from the package `sklearn` ([scikit-learn.org](http://scikit-learn.org/stable/install.html)).
  * For GUI application for plotting coupler curves of G48, `PyQt5` ([pypi.python.org/pypi/PyQt5](https://pypi.python.org/pypi/PyQt5)) and `matplotlib` ([matplotlib.org/](https://matplotlib.org/)) are needed.
  * For installation, just clone or download from [github.com/Legersky/SpatialGraphEmbeddings](https://github.com/Legersky/SpatialGraphEmbeddings).

## Supported graphs
  * 6 vertices: octahedron/cyclohexane (the unique 6-vertex graph with the maximal number of embeddings)
  * 7 vertices: G16a, G16b, G24, G32a, G32b, G48 (all 7-vertex graphs requiring the last Henneberg step being H2,
  the number corresponds to the number of embeddings)
  * 8 vertices: G128, G160

![graphs](http://jan.legersky.cz/public_files/spatialGraphEmbeddings/graphs_7and8vert.png "Supported graphs with 7 and 8 vertices")


If you want to compute the number of embeddings for another minimally rigid graph,
please, provide a method `constructEquations_YOUR_GRAPH` with sphere equations in *graphEmbeddings3D.graphEmbedding*, 
and modify the constructor accordingly.
For the sampling method, subgraphs suitable for sampling must be added to constructor of *graphEmbeddings3D.algRealEmbeddings*.
We appreciate if you share your changes in the [GitHub repository](https://github.com/Legersky/SpatialGraphEmbeddings).

## Tests
`python test_6vert.py` runs the sampling method for octahedron

`python test_7vert.py` verifies that there are edge lengths for G16a, G16b, G24, G32a, G32b and G48 such that all embeddings are real

`python test_8vert.py` verifies that there are edge lengths G128 and G160 have 128 real embeddings 

## Sampling
The scripts in the folder `sampling_scripts` use the proposed method for various graphs and starting edge lengths.

## Coupler curves of G48
This Qt program is launched by `python CouplerCurveG48.py`.

Functionality:
  * loading and saving edge lengths
  * plotting coupler curve of G48
  * computing number of real embeddings of G48 by PHC
  * sampling of parameters for specific subgraphs
  * iterative method for increasing the number of real embeddings
  * export to [Axel](http://axel.inria.fr/)

## Warning
The program strongly depends on PHC computation - this fails sometimes that might cause failure of the program.

## License
Copyright (C) 2018 Jan Legerský

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
