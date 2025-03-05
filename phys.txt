import raytracer.elements as elem
import numpy as np

def refract(direc,normal, n_1,n_2):
    """Snells law in vector form."""

    n = normal/np.sqrt(np.dot(normal,normal))
    i = direc/np.sqrt(np.dot(direc,direc))

    cos_theta_1 = np.dot(n,i)
    sin_theta_2 = (n_1/n_2)*np.sqrt((1-np.dot(n,i)**2))
    sin_theta_1 = (n_2/n_1)*sin_theta_2
    if sin_theta_1 >= n_2/n_1:
        return None
          
    t = n_1/n_2*i - ((n_1/n_2)*cos_theta_1 + np.sqrt(1-(sin_theta_2)**2))*n
    t_norm = t/ np.sqrt(np.dot(t,t))
    return np.array(t_norm)
