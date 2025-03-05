""" rays module contains class that generates ray tracing."""
import numpy as np
import matplotlib.pyplot as plt


class Ray:
     """Represents a ray in 3D space, consisting of a 
        postional vector and a direction vector.
        Allows for adding new points and a consequently updated direction"""

     def __init__(self, pos=None, direc=None):
        """ Initializes two 3d vector instances for the position and direction of a ray. 
            Raises an exception if position is not of length 3. Likewise for direction."""
        if pos is None:
            pos = [0, 0, 0]
        if direc is None:
            direc = [0, 0, 1]   
            
        if len(pos) != 3:
            raise TypeError("Position parameter pos has incorrect size.")
        if len(direc) != 3:
            raise TypeError("Direction parameter direc has incorrect size.") 
        
        self._pos = np.array(pos, dtype=float)
        self._vertices = [self._pos.copy()] #stores list of the original points 
        self._direc = np.array(direc, dtype=float)

     def pos(self):
        """ Returns the postional components as an array pos=[x, y, z]."""
        return self._pos
    
     def direc(self):
        """ Returns the directional components as an array direc=[x,y,z]."""
        return self._direc
    
     def setpos(self, values):
        """ Sets the postional components from an array pos=[x, y, z]."""
        if len(values) != 3:
            raise Exception("Postion parameter pos has incorrect size.")
        self._pos = np.array(values, dtype=float)
        self._vertices.append(self._pos.copy()) #appends new position to list of points ray intercepts

     def setdirec(self, values):
        """Sets the directional components from an array direc=[x, y, z]."""
        if len(values) !=3:
            raise Exception("Direction parameter direc has incorrect size.")
        self._direc = np.array(values,dtype =float)   

     def append(self,pos,direc):
         """Appends the new point and direction."""
         self.setpos(pos)
         self.setdirec(direc)

     def vertices(self):
         """Returns list of points the ray has passed through."""
         return self._vertices

     def __repr__(self):
        """Returns a string representation of the object for debugging"""
        return f"(pos=array({np.array2string(self._pos, separator=',')}), dir=array({np.array2string(self._direc, separator=', ')}))"
        
     def __str__(self):
        """Returns a string representation of the objects for printing"""
        return f"({self._pos[0]:g}, {self._pos[1]:g}, {self._pos[2]:g}, {self._direc[0]:g}, {self._direc[1]:g}, {self._direc[2]:g})"
    
     def copy(self):
        """Returns a deep copy of the Ray instance"""
        return Ray(np.copy(self._pos), np.copy(self._direc)) 
         

class RayBundle:
    """Saliba."""
    def __init__(self, rmax=5., nrings=5, multi=6):
        self._rmax = rmax
        self._nrings = nrings
        self._multi = multi
        self._rays = list(self.rtrings())  # Store generated rays

    def rtrings(self):
        """Generator function for generating radius and theta coordinates."""
        coordinates_list = []
        # Generate coordinates and create rays
        for i in range(1, self._nrings + 1):
            radius = (i / self._nrings) * self._rmax
            num_points = i * self._multi
            for j in range(num_points):
                theta = (j / num_points) * 2 * np.pi
                x = radius * np.cos(theta)
                y = radius * np.sin(theta)
                ray = Ray([x, y, 0])
                coordinates_list.append(ray)

        for ray in coordinates_list:
            yield ray 

        
    def propagate_bundle(self, elements):
        for ray in self._rays:  # Use stored rays
            for optical_element in elements:
                optical_element.propagate_ray(ray)
        return self

    def track_plot(self):
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        
        for ray in self._rays:  # Use stored rays

            vertices = ray.vertices() 

            x = [arr[0] for arr in vertices]
            y = [arr[1] for arr in vertices]
            z = [arr[2] for arr in vertices]
        
            ax.grid()
            ax.plot3D(x, y, z)
        
        plt.show()

        return fig
    
    def rms(self):
        """Calculate the root mean square (RMS) of the radii of the rays in the bundle."""
        rms_gang = []
        for ray in self._rays:
            x, y = ray.pos()[:2]
            r2 = (x**2 + y**2)
            rms_gang.append(r2)
        
        rms_array = np.array(rms_gang)
        rms_array2 = np.append(rms_array,0)

        rms_value = np.sqrt(np.mean(rms_array2))
        return rms_value
    
    def spot_plot(self):
        """plots the genpolar plot at the focal point"""
        from raytracer import elements

        x_coords = []
        y_coords = []

        fig = plt.figure()
        
        for ray in self._rays:
            
            if ray.pos() is not None:
                x,y = ray.pos()[:2]
                x_coords.append(x)
                y_coords.append(y)
        
        plt.figure(figsize=(6, 6))
        plt.scatter(x_coords, y_coords, c='red', marker='o')
        plt.title('Spot Diagram at Focal Point')
        plt.xlabel('X Position')
        plt.ylabel('Y Position')
        plt.grid(True)
        plt.show()

        return fig
