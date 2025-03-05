
Lenses:

class BiConvex:
    """Defines Bi-Convex Optical Element"""

    def __init__(self, z_0, curvature1, curvature2, n_inside, n_outside, aperture, use_plus_case=True):
            """Initializes the parameters defining the Bi-Convex optical element"""
            self._z_0 = z_0
            self._curvature1 = curvature1
            self._curvature2 = curvature2
            self._n_inside = n_inside
            self._n_outside = n_outside
            self._aperture = aperture
            self._use_plus_case = use_plus_case
            
            self._thickness = self.calculate_thickness()
            self._wavelength = 0.000588

    def calculate_thickness(self):
        """Calculates the thickness based on z_0 and aperture"""
        root = np.sqrt((self._aperture / 2) ** 2 + self._z_0 ** 2)
        if self._use_plus_case:
            return 2 * (-self._z_0 + root)
        else:
            return 2 * (-self._z_0 - root)

    def propagate_ray(self, ray):
        """Propagates ray through system"""
        convex_lens1 = SphericalRefraction(z_0=self._z_0, aperture=self._aperture, curvature=self._curvature1, n_1=self._n_outside, n_2=self._n_inside)
        convex_lens1.propagate_ray(ray)
        convex_lens2 = SphericalRefraction(z_0=self._z_0 + self._thickness, aperture=self._aperture, curvature=self._curvature2, n_1=self._n_inside, n_2=self._n_outside)
        convex_lens2.propagate_ray(ray)

    def focal_point(self):
            """Finds the focal point for the biconvex lens using the thick lens equation"""
            n_1 = self._n_outside
            n_2 = self._n_inside
            R_1 = 1 / self._curvature1
            R_2 = 1 / self._curvature2
            t = self._thickness

            # Using the lensmaker's equation for a thick lens
            f_reciprocal = ((n_2 - n_1) * ((1 / R_1) - (1 / R_2) + ((n_2 - n_1) * t) / (n_2 * R_1 * R_2))) / n_1
            f = 1 / f_reciprocal
            location_f = self._z_0 + self._thickness +f
            
            print(self._thickness)
            
            print((153-(self._thickness +self._z_0))/(self._thickness +self._z_0))
            return location_f

    def plot_delta_function_with_points(self):
        """Plots the delta function along with given data points"""

        mod_y_vals = [1, 2, 3, 4, 5]
        origin_z_vals = [153.4628,153.3702,153.2153, 152.9968, 152.7134]
        
        f_0_initial = self.focal_point()
        n = self._n_inside
        r_1 = 1 / self._curvature1
        r_2 = 1 / self._curvature2
        q = (r_2 + r_1) / (r_2 - r_1)
        s_1 = self._z_0 + self._thickness / 2
        s_2 = s_1  # Approximation
        p_initial = (s_1 - s_2) / (s_1 + s_2)

        # Define the K function
        def K(f_0, n, q, p):
            term1 = (n + 2) / (n - 1) * q**2
            term2 = 4 * (n + 1) * q * p
            term3 = (3 * n + 2) * (n - 1) * p**2
            term4 = n**3 / (n - 1)
            return (1 / (4 * f_0)) * (1 / (n * (n - 1))) * (term1 + term2 + term3 + term4)

        def delta_function(z, f_0, p):
            k = K(f_0, n, q, p)
            return np.sqrt(2 * (f_0 - z) / k) 

        popt, _ = curve_fit(lambda z, f_0, p: delta_function(z, f_0, p), origin_z_vals, mod_y_vals, p0=[f_0_initial, p_initial])
        f_0_opt, p_opt = popt
        fitted_y_vals = delta_function(np.array(origin_z_vals), f_0_opt, p_opt)

        z_values = np.linspace(150, 160, 100)
        k_theoretical = K(f_0_initial, n, q, p_initial)
        y_theoretical = np.sqrt(2 * (f_0_initial - z_values) / k_theoretical)

        z_values_extended = np.linspace(min(origin_z_vals) - 5, max(origin_z_vals) + 5, 200)
        fitted_y_vals_extended = delta_function(z_values_extended, f_0_opt, p_opt)
        origin_z_vals_adjusted = np.abs(np.array(origin_z_vals) - f_0_opt)
        z_values_adjusted = np.abs(z_values - f_0_initial)
        z_values_extended_adjusted = np.abs(z_values_extended - f_0_opt)

        fig, ax = plt.subplots()
        ax.plot(mod_y_vals, origin_z_vals_adjusted, 'bo', label='Data Points') 
        ax.plot(fitted_y_vals_extended, z_values_extended_adjusted, 'r-', label=f'Fitted Curve (f_0={f_0_opt:.4f}, s_prime={0.4413})')
        ax.plot(y_theoretical, z_values_adjusted, 'g--', label=f'Theoretical Curve (f_0={f_0_initial:.4f}, s_prime={0.4535}) )') 
        ax.grid()
        ax.set_title('Magnitude of Spherical Aberation Against Distance for Select Bundle', fontsize=20)
        ax.set_xlabel('Initial Ray Distance from Optical Axis (mm)', fontsize=15) 
        ax.set_ylabel('Distance Between Spherical Aberration Intersection Theoretical Focal Point', fontsize=15) 
        ax.legend(loc='upper right')


Analysis:


import matplotlib.patches as patches
plt.rcParams['font.family'] = 'Times New Roman'

def task17_ext():

    # Case 3: Curvature1 > 0 and Curvature2 < 0 (biconvex)
    rb_bc = RayBundle(rmax=5, nrings=5, multi=6)
    bc = BiConvex(z_0=100.0, curvature1=0.02, curvature2=-0.02, n_inside=1.5168, n_outside=1, aperture=50)
    focal_point_bc = bc.focal_point()
    op_bc = OutputPlane(focal_point_bc)

    rb_bc.propagate_bundle([bc, op_bc])
    spot_plot_bc = rb_bc.spot_plot()
    track_plot_bc = rb_bc.track_plot()

    fig, pt = plt.subplots()
    pt.plot(focal_point_bc, 0, 'kx', markersize=7, label='Intercept Point', zorder=5)

    # Plot the tracks of the rays in 2D
    for ray in rb_bc._rays:
        vertices = ray.vertices()
        y = [v[1] for v in vertices]
        z = [v[2] for v in vertices]
        pt.plot(z, y)
    
    pt.legend(['Focal Point'] ,loc='upper right', bbox_to_anchor=(0.88, 1), fontsize='large')
    pt.grid()
    pt.set_title('Ray Paths through Biconvex Lens to Theoretical Focal Point', fontsize = 20)
    pt.set_xlabel('Ray Z Position (mm)', fontsize =15)
    pt.set_ylabel('Ray Y Position (mm)', fontsize =15)

def task18_ext():
    """Saka"""

    bc = BiConvex(z_0=100.0, curvature1=0.02, curvature2=-0.02, n_inside=1.5168, n_outside=1, aperture=50)
    focal_point_bc = bc.focal_point()
    op_bc = OutputPlane(focal_point_bc)
    
    ray1 = rays.Ray([0.1, 0.1, 0])
    ray2 = rays.Ray([0, 0, 0])
    ray3 = rays.Ray([-0.1, -0.1, 0])
    ray_squad = [ray1, ray2, ray3]

    colors = ['r', 'g', 'b']
    ray_labels = [
        'Ray1 = [0.1,0.1,0]', 'Ray2 = [0, 0, 0]', 'Ray3 = [-0.1, -0.1, 0]'
    ]

    fig, pt = plt.subplots()

    for idx, r in enumerate(ray_squad):
        bc.propagate_ray(r)
        op_bc.propagate_ray(r)
        vertices = r.vertices()
        
        f = bc.focal_point()

        y = [arr[1] for arr in vertices]
        z = [arr[2] for arr in vertices]
        color = colors[idx]

        pt.plot(z, y, color=color, label=ray_labels[idx])

    pt.plot(focal_point_bc, 0, 'kx', markersize=7, label='Focal Point', zorder=5)
    pt.grid()
    pt.set_title('Paraxial Ray Paths through Bifocal Lens to Theoretical Focal Point', fontsize =20)
    pt.set_xlabel('Ray Z Position (mm)', fontsize =15)
    pt.set_ylabel('Ray Y Position (mm)', fontsize =15)
    pt.legend(loc='upper right', bbox_to_anchor=(1, 1), fontsize=15)

    return fig, f

def task19_ext():
    """Odegaard"""
    bc = BiConvex(z_0=100.0, curvature1=0.02, curvature2=-0.02, n_inside=1.5168, n_outside=1, aperture=50)
    focal_point_bc = bc.focal_point()
    op_bc = OutputPlane(250)
    
    ray1 = rays.Ray([0, 5, 0])
    ray2 = rays.Ray([0, 4, 0])
    ray3 = rays.Ray([0, 3, 0])
    ray4 = rays.Ray([0, 2, 0])
    ray5 = rays.Ray([0, 1, 0])
    ray6 = rays.Ray([0, 0, 0])
    ray7 = rays.Ray([0, -1, 0])
    ray8 = rays.Ray([0, -2, 0])
    ray9 = rays.Ray([0, -3, 0])
    ray10 = rays.Ray([0, -4, 0])
    ray11 = rays.Ray([0, -5, 0]) 

    ray_squad = [ray1, ray2, ray3, ray4, ray5, ray6, ray7, ray8, ray9, ray10, ray11]

    # Define ray labels
    ray_labels = [
        'Ray1 = [0, 5, 0]', 'Ray2 = [0, 4, 0]', 'Ray3 = [0, 3, 0]', 
        'Ray4 = [0, 2, 0]', 'Ray5 = [0, 1, 0]', 'Ray6 = [0, 0, 0]', 
        'Ray7 = [0, 1, 0]', 'Ray8 = [0, 2, 0]', 'Ray9 = [0, -3, 0]', 
        'Ray10 = [0, -4, 0]', 'Ray11 = [0, -5, 0]']

    colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    def get_color(value):
        return colors[int(value) % len(colors)]

    fig, pt = plt.subplots()
    pt.plot(focal_point_bc, 0, 'kx', markersize=10, label='Theoretical Focal Point', zorder=5)

    # Plot rays with specific colors based on modulus
    for idx, r in enumerate(ray_squad):
        bc.propagate_ray(r)
        op_bc.propagate_ray(r)
        vertices = r.vertices()
        y = [arr[1] for arr in vertices]
        z = [arr[2] for arr in vertices]
        value = abs(y[0])  # Using the absolute value of the initial y position for modulus
        color = get_color(value)
        pt.plot(z, y, color=color, label=ray_labels[idx]) 

    # Add grid, title, and labels
    pt.grid()
    pt.set_title(r'Select Ray Bundle Paths through BiConvex Lens to Output Plane (Zoom)', fontsize =20)
    pt.set_xlabel('Ray Z Position (mm)', fontsize =15)
    pt.set_ylabel('Ray Y Position (mm)', fontsize =15)
    pt.legend(loc='upper right', bbox_to_anchor=(1, 1), fontsize=13, ncol=3)

    return fig

def task20_ext():
    """Rice""" 
    bc = BiConvex(z_0=100.0, curvature1=0.02, curvature2=-0.02, n_inside=1.5168, n_outside=1, aperture=50)
    bc.plot_delta_function_with_points()

def task21_ext(): 
    """Gabriel"""
    bc = BiConvex(z_0=100.0, curvature1=0.02, curvature2=-0.02, n_inside=1.5168, n_outside=1, aperture=50)
    focal_point_bc = bc.focal_point()
    op_bc = OutputPlane(250)
    
    ray1 = rays.Ray([0, 3.8516, 0])
    ray2 = rays.Ray([0, -3.8516, 0])
    ray_squad = [ray1, ray2]
    ray_labels = ['Ray1 = [0, + (initial RMS ), 0]', 'Ray2 = [0, - (initial RMS), 0]']

    colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    def get_color(value):
        return colors[int(value) % len(colors)]

    fig, pt = plt.subplots()
    pt.plot(focal_point_bc, 0, 'kx', markersize=10, label='Theoretical Focal Point', zorder=5)
    pt.plot(153.036, 0, 'gx', markersize=10, label='Simulated Focal Point', zorder=5)   
    for idx, r in enumerate(ray_squad):
        bc.propagate_ray(r)
        op_bc.propagate_ray(r)
        vertices = r.vertices()
        y = [arr[1] for arr in vertices]
        z = [arr[2] for arr in vertices]
        value = abs(y[0])  
        color = 'b'
        pt.plot(z, y, color=color, label=ray_labels[idx])

    pt.grid()
    pt.set_title(r'Initial Rms of Bundle Ray Paths through BiConvex Lens to Output Plane', fontsize=20)
    pt.set_xlabel('Ray Z Position (mm)', fontsize=15)
    pt.set_ylabel('Ray Y Position (mm)', fontsize=15)
    pt.set_yticks([3.8516, -3.8516])
    pt.set_yticklabels(['rms radius = 3.8516', 'rms radius = -3.8516'], fontsize=11)

    S = 100 + (6.155281280883031 / 2)
    
    pt.plot([0, S], [0, 0], 'r--', label='S length')
    pt.legend(loc='upper right', bbox_to_anchor=(1, 1), fontsize=13, ncol=1)
    return fig


def task22_ext():
    rb_bc = RayBundle(rmax=9, nrings=5, multi=6)
    rms = rb_bc.rms()
    print(rms)
    
    fig, ax = plt.subplots()  
    rb_bc.spot_plot()
    ax = plt.gca()
    
    circle = patches.Circle((0, 0), radius =3.8516, edgecolor='blue', facecolor='none', linestyle='--', linewidth=2, label='Radius 3.8516')
    ax.add_patch(circle)
    ax.legend(loc='upper right', bbox_to_anchor=(1, 1), fontsize=13, ncol=1)

#task17_ext()

    #task18_ext()

    #task19_ext()

    #task20_ext()

    #task21_ext()

    #task22_ext()
