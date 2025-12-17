This implementation is part of a Bachelor thesis on "Max Cut on Planar Graphs"
It was coded and tested in Python 3.11

Additional packages required are "scipy" & "networkx" (as described in the thesis), 
  "drawsvg" for the display of graphs as vector graphics, and "tabulate" for additional output formatting.
All the packages can be installed via pip.

The module PlanarCuts contains a datastructure for planar graphs and implementations for cut algorithms on it.
This includes:
- max cut
- min cut  (both based on an algorithm by Liers and Pardella)
- approximated max cut
- approximated min cut  (both based on local search)
- random cut

To use it, use the following import:
~~~from PlanarCuts import PlGraph
~~~

Generate a planar graph for example by running
~~~PlGraph.random()
~~~

The files singleRun.py and multiRun.py serve both as an example on the methods provided by the PlanarCuts module
  and to test generating random graphs and see how the different algorithms perform.
To change the style of random graph generated, change the parameters at the start of these files.
singleRun will create a single random graph, run all avaible cut algorithms on it (the approximations multiple times
  to calculate averages), print some stats and draw vectorgraphic representations for the original graph 
  and each of the cut types
multiRun will repeatedly create random graphs with the same parameters, run max cut and several attempts of 
  approximating max cut, while calculating averages and other stats, eventually writing them to log files
