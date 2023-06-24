# SWCNT exciton T2 calculation

This repo was first developed by Diego Barrutia and now maintained by Jing Trerayapiwat. Diego participated in the 2022 Boston University Physics Research Opportunity (BU-PRO) program under the supervision of Jing and professor Sahar Sharifzadeh. Here, we implemented a model to calculate quantum decoherence of an exciton in single-walled carbon nanotube as a result of scattering with lateral acoustic phonons. For more details about the model see: Galland, C., Högele, A., Türeci, H. E., & Imamoglu, A. (2008). Non-Markovian decoherence of localized nanotube excitons by acoustic phonons. Physical Review Letters, 101(6), 1–4. https://doi.org/10.1103/PhysRevLett.101.067402 

SWCNT class calculates and stores intrinsic physical properties of SWCNT. run_X is the main program which plots exciton decoherence vs time.

TODO: link the SWCNT class to a database and obtain SWCNT exciton confinement length instead of asking user for an input.

# Theoretical Background

Polarization is linearly dependent to the linear susceptibility, $\chi(t)$, which describes the electron’s response to a perturbation. For an acoustic lateral phonon, the frequency, $\omega (q) = v \cdot q$, is linearly dependent on the phonon momentum, $q$, with a coefficient $v$, the sound velocity. Within the independent boson model, the linear susceptibility of the QD in response to a $\delta$-function-shaped laser pulse at $t = 0$ is decomposed into a temperature-dependent and temperature-independent parts,

$$\chi(t)  = -ie^{-i\bar{\Omega}t}\chi_{T}(t)\chi_{0}(t)$$,     

$$\chi_{T}(t)  \propto i\mathrm{exp}(\sum_{q}|\gamma(q)|^{2}[-n(q)|e^{-i\omega(q)t}-1|^{2}])$$,    
$$\chi_{0}(t)  \propto i\mathrm{exp}(\sum_{q}|\gamma(q)|^{2}[e^{-i\omega(q)t}-1])$$,
$\gamma$(q) = g/$\omega$, and $n(q)=(e^{\frac{h\omega(q)}{k_BT}}-1)^{-1}$ is the phonon occupation number. The polaron-shifted transition energy, $\bar{\Omega}$ is approximated as the bare energy of the excited state without phonon corrections. We describes the exciton-phonon matrix elements as a result of coupling, $g(q)$, along the repeating $z$ direction of SWCNT as a product of $G(q)$, the deformation coupling, and $F(q)$, the form factor, ie.  $g=G(q)F(q)$, with

$$F(q) = \int dz|\Psi^{exc}(z)|^2e^{iqz}$$,
$$G(q) = \frac{D_sq}{\sqrt{2\rho L \hbar \omega(q)}}$$,
$$\Psi(z) = \pi^{-\frac{1}{4}}\sigma^{-\frac{1}{2}} e^{-\frac{z^2}{2\sigma^2}}$$,
with $\Psi$, the exciton wave function fitted with a Gaussian. $L$ is the unit cell length, $\rho$ the linear mass density, and $D_s$ the deformation potential, a measure of exciton coupling strength to the stretching mode. Upon simplification we obtain,

$$\gamma(q) =  \frac{D_sq}{\sqrt{2\rho L \hbar v^3q}} e^{-\frac{q^2\sigma^2}{4}}$$.

To describes the decay of polarization we take the norm of the linear susceptibility and obtain the following form:

$$|\chi_0| = \prod_{q}e^{-|\gamma|^2(1-cos(wt))} = \frac{L}{\pi} \int_{q=0}^{BZ}e^{-|\gamma|^2(1-cos(wt))}$$,
$$|\chi_{T}| = \prod_{q}e^{-2n|\gamma|^2(1-cos(wt))} = \frac{L}{\pi} \int_{q=0}^{BZ}e^{-2n|\gamma|^2(1-cos(wt))}$$,
$$|\chi| = |\chi_{0}\chi_{T}| = \prod_{q}e^{-|\gamma|^2(1-cos(wt)(1+2n))} = \frac{L}{\pi} \int_{q=0}^{BZ}e^{-|\gamma|^2(1-cos(wt)(1+2n))}$$,
where $\frac{\pi}{L}$, defines the edge of the Brillouin Zone.  The above equations can be easily implemented using Scipy and Numpy.

The linear mass density increases with larger SWCNT diameter. The linear mass density of a (n,m) SWCNT chirality can be calculated from (Dresselhaus, G., Dresselhaus, M. S., & Saito, R. (1998). Physical Properties of Carbon Nanotubes (1st ed.). Imperial College Press.):
$$\rho = m_C\frac{ N_{atoms/cell}}{L}$$,
$$N_{atoms/cell} = \frac{4(n^2+nm+m^2)}{d_R}$$,
$$L = \frac{3a}{d_R} \sqrt{n^2+nm+m^2},
$$d_R = gcd(2m+n,2n+m)$$,
with a = C-C aromatic bond length, 1.55 \AA{}, $m_C$ the mass of a carbon atom. The dependence of $\rho$ on (n,m) means it is possible to tune $T_2$ as a function of chirality.




