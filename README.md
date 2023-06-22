# SWCNT exciton T2 calculation

This work is done by Diego Barrutia and Jing Trerayapiwat. Diego participated in the 2022 Boston University Physics Research Opportunity (BU-PRO) program under the supervision of Professor Sahar Sharifzadeh. Here, we developed a model for quantum coherence time between carbon nanotube excitons and acoustic phonons from "Non-Markovian decoherence of localized nanotube excitons by acoustic phonons" (link below). The highlight of this experience came at the summer’s end when we derived a novel relationship between the chirality of the carbon nanotube and coherence time. 

Research paper link: https://arxiv.org/abs/0802.2046

equation.py: Using the Brillouin Zone and Riemann Sum to convert sigma to an integral, we are able to graph |X(t)| vs Time (ps). This is the first graph. Additionally , the long term of the X(t) is derived form the Indepent Boson Model (IBM) and plotted. We use different chirality values to understand its effects on the decoherence behavior.

Chirality: Using the magnitude of unit vector and unit cell length, we are able to calculate for chirality. This chirality is later applied in equation.py.

Norm Comparison (Coding and Analytical).py: To double check out work, we derive the norm analytically and compared it with the code. As we can see in the model produced, both graph fit confirming the reliability of the code. 

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




