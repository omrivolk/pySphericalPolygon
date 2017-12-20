pySphericalPolygon
==================

Perform point in polygon operations on a sphere. Handles any polygon which can be viewed on a single hemisphere including around the poles. 

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
