import numpy as np

" CONSTANTS "
# carbon-carbon bond length in aromatic systems 
a = np.sqrt(3) * 1.44e-10  # m
# mass of carbon atom
mass_C = 1.99e-26  # kg
h_bar = 1.054571817e-34 # Js
eV2J = 1.60218e-19 # eV/J


class SWCNT:
    """
    A class to represent a spin-active nuclear bath.
    Takes in positions and randomly assigns isotopes.

    Attributes
    ----------
    atoms: 
        ase's Atoms object to get geometry
    spin_table: 
        contains atomic number, atomic mass, percent abundance, nuclear_spin
        of species present in self.atoms
    bath_geometry: 
        contains atomic symbol, isotope, xyz positions, nuclear spin(I),
        and the distance from e- Spin

    Methods
    -------
    apply_rcutoff(e_spin_position, cutoff_radius):
        remove sites beyond a certain radius from e- spin
        use when e- spin has been defined
    
    Usage Guide
    -------
        SWCNT=bath.from_file('CONTCAR_NO2',(1,1,3))
        SWCNT.bath_geometry
        SWCNT.apply_r_cutoff(origin=[0,0,0], r=40)
        SWCNT.bath_geometry
    """
    def __init__(self, n: int, m: int):
        """
        Structural properties of (n,m) SWCNT.

        Parameters
        ----------
            n:
                circumferntial folds. Pure (n,0) yields zigzag SWCNT.
            m:
                longtitudinal folds. m introduces chirality to zigzag SWCNTs.
        """
        self.n=n
        self.m=m
        #for unit cell length of SWCNT
        self.L = get_unit_cell_length(self.n, self.m) #m
        #get linear mass density
        self.rho= get_n_atoms_per_cell(self.n,self.m) / self.L #Kg/m
        #TODO: hook up to a database and make E11 and exciton localization class properties 
        # confinement length (m)
        self.sigma = 3e-9 #m #exciton
        # energy of the excited state (Hz)
        self.omega = 1.41 * eV2J / h_bar # Hz
        # sound velocity in SWCNT (m/s)
        # approx. as velocity in graphene (Xu, B. App. Phys. Lett., 96(18), 183108.)
        self.v_s = 19900 # m/s
        # deformation potential (J)
        # approx. as D in graphene (Yu, J. JCP, 103(15), 6697–6705.)
        self.D_s = 2.243047191e-18 # J


    def get_unit_cell_length(n: int, m: int):
        """
        for present (n,m) use a formula to give unit cell length
        
        Parameters
        ----------
                n:
                    circumferntial folds. Pure (n,0) yields zigzag SWCNT.
                m:
                    longtitudinal folds. m introduces chirality to zigzag SWCNTs.
        """
        return 3 * a / (get_dr(n,m)) * np.sqrt(n**2 + n*m + m**2) #m


    def get_n_atoms_per_cell(n: int, m: int):
        """
        for present (n,m) use a formula to give number of atoms per cell
        
        Parameters
        ----------
                n:
                    circumferntial folds. Pure (n,0) yields zigzag SWCNT.
                m:
                    longtitudinal folds. m introduces chirality to zigzag SWCNTs.
        """
        return 4* np.sqrt(n**2 + n*m + m**2) / get_dr(n,m) #Kg/cell

        
    def get_dr(n: int, m: int):
        """
        for L, rho calculations. d_r defines redundancy along SWCNT vector
        
        Parameters
        ----------
                n:
                    circumferntial folds. Pure (n,0) yields zigzag SWCNT.
                m:
                    longtitudinal folds. m introduces chirality to zigzag SWCNTs.
        """
        return np.gcd(2*m + n, 2*n + m)  #unitless
    
        #provided user input set exciton confinement
    def set_sigma(self,sigma):
        self.sigma=sigma

    #provided user input set first excited state
    def set_omega(self,omega):
        self.omega=omega
