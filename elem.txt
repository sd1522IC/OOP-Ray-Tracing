import numpy as np
import raytracer.rays as rays
import raytracer.physics as phys

class OpticalElement:
    
        def intercept(self, ray):
                raise NotImplementedError('intercept() needs to be implemented in derived classes')
        def propagate_ray(self, ray):
                raise NotImplementedError('propagate_ray() needs to be implemented in derived classes')
        
class SphericalRefraction(OpticalElement):
    """Defines Spherical Optical Element"""
    def __init__(self, z_0 = 10.0, aperture = 5.0, curvature = 0.2, n_1 = 1.0, n_2= 1.5):
           """Initialises the parameters defining the spherical optical element"""
           self._z_0 = z_0
           self._aperture = aperture
           self._curvature = curvature
           self._n_1 = n_1
           self._n_2 = n_2

    def intercept(self, ray = rays.Ray()):
        
        """Defines the interception point on the surface."""
        pos = ray.pos()
        k = ray.direc()
        k = k/np.linalg.norm(k)

        if self._curvature == 0:
            r = ray.pos()
            if k[2] == 0:
                return None  # If the direction vector is parallel to the plane
        
            lambda_val = (self._z_0 - r[2]) / k[2]
            termination_point = r + lambda_val * k

            return termination_point
        
        else: 
            r = pos - [0,0,self._z_0 +1/self._curvature]
            R = 1/self._curvature

            if (np.dot(r,k))**2 < np.linalg.norm(r)**2-R**2:
                return None
            
            l1 = -np.dot(r,k) - np.sqrt((np.dot(r,k)**2-(np.abs(np.linalg.norm(r)**2) - R**2)))
            l2 = -np.dot(r,k) + np.sqrt((np.dot(r,k)**2-(np.abs(np.linalg.norm(r)**2) - R**2)))
            
            if l1 > 0: #then l2 must be larger as l2>l1 always
                if self._curvature > 0:
                    intersect = pos + l1*k
                else:
                    intersect = pos + l2*k
                if intersect[0]**2 + intersect[1]**2 < self._aperture**2: #if inside aperture return point
                    return intersect
                else:
                    return None
            elif l2 < 0: #so l1 is also < 0
                return None
            else: #last case if l1<0 and l2>0
                intersect = pos + l2*k
                if intersect[0]**2 + intersect[1]**2 < self._aperture**2: #if inside aperture return point
                    return intersect
                else:
                    return None
            
    def propagate_ray(self, ray = rays.Ray()):
        """How Ray direction changes after refracting on lens,
        appending the new position and direction."""
        intercept_point = self.intercept(ray)
        
        if intercept_point is None:
            return None
 
        if self._curvature == 0:
            normal = (0, 0, 1)
        else:
            R = 1 / self._curvature
            

            center = np.array([0, 0, self._z_0 + R])
            normal = intercept_point - center
    
            normal = normal / np.sqrt(np.dot(normal, normal))

        new_direction = phys.refract(ray.direc(), normal, self._n_1, self._n_2)
        
        if new_direction is None:
            return None
        
        else:
            new_direction = new_direction/ np.sqrt(np.dot(new_direction, new_direction))
            ray.append(intercept_point, new_direction)

        return ray
    
 
    def z_0(self):
        """Returns z-axis intecept of optical element"""
        return self._z_0
    
    def aperture(self):
        """Returns aperture radius of optical element"""
        return self._aperture
    
    def curvature(self):
        """Returns curvature element of optical element"""
        return self._curvature
    
    def n_1(self):
        """Returns n_1 refreactive index element of optical element"""
        return self._n_1
    
    def n_2(self):
        """Returns n_2 refractive index element of optical element"""
        return self._n_2
    
    def focal_point(self):
        z_0 = self._z_0
        n_1 = self._n_1
        n_2 = self._n_2
        R_1 = 1/(self._curvature)
        f_recirocal= (n_1/n_2 -1)*(-1/R_1)
        f = 1/f_recirocal
        location_f = z_0 + f
        return location_f
        

class OutputPlane(OpticalElement):
    """Defines a lens where rays terminate"""
    def __init__(self, z_0):
        """Initializes the z-position of the Output plane"""
        self._z_0 = z_0

    def intercept(self, ray):
        """Finds the interception point of the ray with the Output plane"""
        r = ray.pos()
        k = ray.direc()
        
        if k[2] == 0:
            return None  # If the direction vector is parallel to the plane
        
        lambda_val = (self._z_0 - r[2]) / k[2]
        termination_point = r + lambda_val * k

        return termination_point
    
    def propagate_ray(self, ray):
        """Propagates the ray to the OutputPlane without changing its direction"""
        termination_point = self.intercept(ray)
        if termination_point is None:
            return None
        
        ray.append(termination_point, ray.direc())
        return ray
    
    