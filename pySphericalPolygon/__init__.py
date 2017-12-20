# -*- coding: utf-8 -*-

"""
--- pySphericalPolygon ---

Python module to handle points in spherical polygon operations.

This module, unlike existing modules for the same operation,
does not require specifying a point inside the polygon.

A spherical polygon is a closed geometric figure on the 
surface of a sphere which is formed by the arcs of great circles.
Unlike a polygon on a 2D plain, where the outside of a polygon is 
of infinte area, on a sphere both the inside and outside are finite.
Therefore common algorithms will require a point inside the polygon
to decide which is in and which is out. This module uses
a diffrent algorithm which assumes the smaller area is the inside.
The disadvantage of this algorithm is that it can't handle 
polygons which can't be viewd on a single hemasphire.
Such polygons are extremly rear in most applications.

Psuedo code:
    
    convert the polygon to 3D cartesian coordinates and approximate it's mean

    project the polygon onto a plain using gnomic projection with it's mean as the projection point

    project the points which needs to be checked onto the same plain

    use a 2D point in polygon algorithm 

This works becuse in gnomic projection great circles are projected as stright lines

By: Omry V, October 2016, omryv@post.bgu.ac.il

License: MIT 

"""

from matplotlib.path import Path
import numpy as np
import pyproj


class SpericalPolygon(object):

    '''
    A class to hold to the spherical polygon
    '''

    def __init__(self,vertices,radians=False,convention='geo'):

        '''
        Create a spherical polygon instance 

        vertices : a (N,2) list or numpy array of floats or ints


                [[ tetha_1, phi_1 ],
                 [ tetha_2, phi_2 ],
                 [ tetha_3 ,phi_3 ],
                 ...               ]

                tetha = latitude
                phi = longtitude

        radians : bool. Is the input given in radians (True) or degrees (False)

        convention : "geo" or "math". Which convention is Geogrphic or mathematical conventions

            "geo" - 

                (-π/2 rad) -90°  ≤ tetha ≤ 90°  (π/2 rad)
                (-π   rad) -180° ≤  phi  ≤ 180° (π   rad)

            "math" - 

                (0 rad) 0° ≤ tetha ≤ 180° (π rad)
                (0 rad) 0° ≤  phi  ≤ 360° (2π rad)
        '''

        #validate the input
        self.vertices_geo,self.vertices_math_rad,self.radians= self.validate_input(vertices,radians,convention)
        
        #find it's mean
        mean_tetha,mean_phi=mean_and_r_sphrical_polygon(self.vertices_math_rad)
        
        #convert the mean to geographical lat/lon
        mean_phi=180.-np.rad2deg(mean_phi)
        mean_tetha=90.-np.rad2deg(mean_tetha)

        #create a gnomic projection
        self.proj=pyproj.Proj('+proj=gnom +lat_0=%s +lon_0=%s  +R=1.'%(mean_tetha,mean_phi))
        
        tethas,phis=self.vertices_geo.T
        
        #project the vertices onto a plain
        self.projected_vertices=np.asarray(self.proj(phis,tethas,radians=radians)).T
        
        #create a 2D polygon from the projected vertices
        self.projected_polygon=Path(self.projected_vertices)
             
    def validate_input(self,points,radians,convention):

        '''
        Validate the input and return vertices in both mathematical
        convention for finding the mean and geogrphical convention
        for projecting.
        '''

        if type(radians)!=bool:
            raise ValueError('radians must be of type bool')            
                
        try:
            if type(points)!=np.ndarray:
                points=np.asarray(points,dtype=float)
            elif points.dtype!=float:
                points=points.astype(float)
            if (points.shape[0]<3) or (points.shape[1]!=2) or len(points.shape)!=2:
                raise
        except Exception as e:
            raise ValueError('Input must be a (N,2) list or numpy array of floats or ints')   
        
        
        if convention=='geo':
            
            points_geo=points+0.
            
            if radians:
                points_rad=points+0.
            else:
                points_rad=np.deg2rad(points)
                
            points_math_rad=np.zeros_like(points_rad)
            points_math_rad[:,0]=np.pi/2.-points_rad[:,0]
            points_math_rad[:,1]=np.pi-points_rad[:,1]
                
        elif convention=='math':
            if radians:
                points_math_rad=points+0.
                points_geo=np.zeros_like(points)
                points_geo[:,0]=np.pi/2.-points[:,0]
                points_geo[:,1]=np.pi-points[:,1] 
            else:
                points_math_rad=np.deg2rad(points)
                points_geo=np.zeros_like(points)
                points_geo[:,0]=90.-points[:,0]
                points_geo[:,1]=180.-points[:,1]        
        else:
            raise ValueError('convention must be either "geo" or "math"')

        return points_geo,points_math_rad,radians
    
    def contains_points(self,points,radians=False,convention='geo'):
        
        '''
        Returns a bool array which is *True* if the path contains the
        corresponding point.

        radians : bool. Is the input given in radians (True) or degrees (False)

        convention : "geo" or "math". Which convention is Geogrphic or mathematical conventions

            "geo" - 

                (-π/2 rad) -90°  ≤ tetha ≤ 90°  (π/2 rad)
                (-π   rad) -180° ≤  phi  ≤ 180° (π   rad)

            "math" - 

                (0 rad) 0° ≤ tetha ≤ 180° (π rad)
                (0 rad) 0° ≤  phi  ≤ 360° (2π rad)
        '''

        #validate the input
        if type(radians)!=bool:
            raise ValueError('radians must be of type bool')            
                
        try:
            if type(points)!=np.ndarray:
                points=np.asarray(points,dtype=float)
            elif points.dtype!=float:
                points=points.astype(float)
            if (points.shape[0]<1) or (points.shape[1]!=2) or len(points.shape)!=2:
                raise
        except Exception as e:
            raise ValueError('Input must be a (N,2) list or numpy array of floats or ints')   
        
        if convention=='geo':
            pass
        elif convention=='math':
            #convert to geogrphic convention for projection
            if radians:
                points[:,0]=np.pi/2.-points[:,0]
                points[:,1]=np.pi-points[:,1] 
            else:
                points[:,0]=90.-points[:,0]
                points[:,1]=180.-points[:,1]        
        else:
            raise ValueError('convention must be either "geo" or "math"')

        #if needed, change units (radians or degrees) to match the polygon units
        if (self.radians) and not(radians):
            points = np.deg2rad(points)
        
        elif not(self.radians) and (radians):
            points = np.rad2deg(points)
        
        tethas,phis=points.T
        
        #project the points to a plain
        projected_points=np.asarray(self.proj(phis,tethas,radians=self.radians)).T
        
        #preform a 2D point in polygon operation
        return self.projected_polygon.contains_points(projected_points)
        
        
    def plot_test(self,resolution=5):
        
        '''
        Show a plot of the earth and polygon with points tested for beeing inside.

        resolution : a positive number. The grid resolution in degrees for the test. 

        '''

        import matplotlib.pyplot as plt
        try:
            from mpl_toolkits.basemap import Basemap
        except Exception as e:
            print(e)
            raise ImportError("This function requiers Basemap, please install and try again")

        lons,lats=np.meshgrid(np.arange(-180,180.01,resolution),np.arange(-90,90.01,resolution))
        points=np.asarray([lats.flatten(),lons.flatten()]).T
        
        points_in=self.contains_points(points,radians=False)
        
        fig,ax =plt.subplots(figsize=(10,10))
        
        m=Basemap()
        
        plot_vertices=self.vertices_geo
        if self.radians:
            plot_vertices=np.rad2deg(self.vertices_geo)
        
        for i in range(len(plot_vertices)-1):
            m.drawgreatcircle(plot_vertices[i][1],plot_vertices[i][0],plot_vertices[i+1][1],plot_vertices[i+1][0],c='b',lw=3)
        m.drawgreatcircle(plot_vertices[i+1][1],plot_vertices[i+1][0],plot_vertices[0][1],plot_vertices[0][0],c='b',lw=3)
        
        m.drawcoastlines()
        m.drawmeridians(np.arange(-180,180,10))
        m.drawparallels(np.arange(-90,90,10))        
        
        for i,point in enumerate(points):
            if points_in[i]:
                c='b'
            else:
                c='r'
            m.scatter(point[1],point[0],latlon=True,s=20,c=c,lw=0)
        plt.show()

def mean_and_r_sphrical_polygon(points):
    
    ''' 
    Approximate the mean point of a spherical polygon using 3D cartesian vectors
    Coordinates are given in radinas using mathematical convention:

        0 ≤ tetha ≤ π rad
        0 ≤ phi < 2π rad

    points : a (N,2) numpy array of floats:

        [[ tetha_1, phi_1 ],
         [ tetha_2, phi_2 ],
         [ tetha_3 ,phi_3 ],
         ...               ]


    Example:

    #a square around the north pole
    >>> points = np.array([[ 0.17453293,  0.        ],
                           [ 0.17453293,  1.57079633],
                           [ 0.17453293,  3.14159265],
                           [ 0.17453293,  4.71238898]])

    >>> mean_and_r_sphrical_polygon(points)
    (0.0, 2.3561945188453621)

    '''

    tethas,phis=points.T
    
    #convert to cartesian #D
    x=np.cos(phis)*np.sin(tethas)
    y=np.sin(phis)*np.sin(tethas)
    z=np.cos(tethas)
    
    #get the mean
    mean = np.asarray([np.mean(x),np.mean(y),np.mean(z)])
    mean=mean/np.linalg.norm(mean)
    
    #find the distance to each vertice and get the max
    dx=x-mean[0]
    dy=y-mean[1]
    dz=z-mean[2]
    c=np.sqrt(dx**2+dy**2+dz**2)
    d_angle=2*np.arcsin(c.max()/2.)
    
    #check the size of the polygon
    if d_angle>np.pi/2.:
        raise ValueError('''This polygon is too big for gnomic projection (%s>%s). This library can only handle polygons which can be viewed on a single hemisphere.''' % (d_angle, np.pi/2.))
    #convert back to spherical coordinates 
    mean_phi = np.arctan2(mean[1],mean[0])
    mean_tetha = np.arccos(mean[2])
        
    return  mean_tetha,mean_phi
