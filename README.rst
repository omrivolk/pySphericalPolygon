pySphericalPolygon
==================

Perform point in polygon operations on a sphere. 

Handles any polygon which can be viewed on a single hemisphere including self intersecting polygons and polygons around the poles. 

This module, unlike existing modules for the same operation, does not require specifying a point inside the polygon.

A spherical polygon is a closed geometric figure on the surface of a sphere which is formed by the arcs of great circles. Unlike a polygon on a 2D plain, where the outside of a polygon is of infinte area, on a sphere both the inside and outside are finite. Therefore common algorithms will require a point inside the polygon to decide which is in and which is out. This module uses a diffrent algorithm which assumes the smaller area is the inside. The disadvantage of this algorithm is that it can't handle  polygons which can't be viewd on a single hemasphire. Such polygons are extremly rear in most applications.

Psuedo code:
.. code-block::

	convert the polygon to 3D cartesian coordinates and approximate it's mean
	project the polygon onto a plain using gnomic projection with it's mean as the projection point
	project the points which needs to be checked onto the same plain
	use a 2D point in polygon algorithm 
    
This works becuse in gnomic projection great circles are projected as stright lines

Install
-----

.. code-block:: python

  pip install pySphericalPolygon
  
Usage
-----

.. code-block:: python

  import pySphericalPolygon as pysp


Create a spherical polygon with vertices [[tetah_1,phi_1],[tetah_2,phi_2]...].

For geographical purpuses tetha is latitude and phi is longtitude.

.. code-block:: python

  sp = pysp.SpericalPolygon([[0,0],
                           [80,30],
                           [10,60]])

Check if a point is inside

.. code-block:: python

	print sp.contains_points([[30,30]])
	[ True]


Check many points at once

.. code-block:: python

	print sp.contains_points([[30,30],[-30,30],[-90,40]])
	[ True False False]


Both vertices and/or points may be specified in radians

.. code-block:: python

	print sp.contains_points([[0.52359878,0.52359878],[-0.52359878,0.52359878],[-1.57079633,0.6981317]],radians=True)
	[ True False False]
	

Convention deafult is geographic:


	(-π/2 rad) -90°  ≤ tetha ≤ 90°  (π/2 rad)

	(-π   rad) -180° ≤  phi  ≤ 180° (π   rad)

But mathematic convetion is supprted too:


	(0 rad) 0° ≤ tetha ≤ 180° (π rad)

	(0 rad) 0° ≤  phi  ≤ 360° (2π rad)

.. code-block:: python
	
	  sp = pysp.SpericalPolygon([[90,0],
                           [10,30],
                           [80,60]],convention='math')
