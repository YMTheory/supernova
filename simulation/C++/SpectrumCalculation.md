# How to calculate PDF

### Two dimensional PDF considering the massive neutrino effect

- If no masses:

$$
\frac{\mathrm{d}N}{\mathrm{d}T_p} = C \int_0^\infin\mathrm{d}T'G(T'_p; T_p,\delta T_p)\times \int_0^\infin\mathrm{d}E_\nu F(E_\nu) \frac{\mathrm{d}\sigma}{\mathrm{d}T_p'(E_\nu, T_p')}.
$$



- Time of flight:

$$
\Delta t = 5.14 \times \left(\frac{m_\nu}{\mathrm{eV}}\right)^2 \times \left(\frac{10\,\mathrm{MeV}}{E_\nu} \right)^2 \times \left(\frac{D}{10\,\mathrm{kpc}}.\right)
$$

with fixed neutrino mass $m_\nu$, the time of flight only depends on neutrino energy.



- Detected event counts in unit time and kinetic energy interval $\left[T_1, T_2\right]$

$$
N(t,T_1,T_2) = \int_{T_1}^{T_2} \mathrm{d}T \int_{t_\mathrm{start}}^{t} \int_{E_\nu^\mathrm{min}}^{E_\nu^\mathrm{max}} F(E_\nu, t') \times \frac{\mathrm{d}\sigma}{\mathrm{d}T}(E_\nu, T)\times \delta\left[t - (t'+\Delta t)\right] \mathrm{d}E_\nu \mathrm{d}t'.
$$

- Detected event counts in unit time and visible energy interval  $\left[E_1, E_2\right]$:
  $$
  \begin{align}
  N(t, T_1(E_1), T_2(E_2)) = N(t, g(E_1), g(E_2)) = &\int_{g(E_1)}^{g(E_2)} g'(E)\mathrm{d}E  \int_{t_\mathrm{start}}^{t} \int_{E_\nu^\mathrm{min}}^{E_\nu^\mathrm{max}} F(E_\nu, t') \nonumber \\
  &\times \frac{\mathrm{d}\sigma}{\mathrm{d}T}(E_\nu, T)\times \delta\left[t - (t'+\Delta t)\right] \mathrm{d}E_\nu \mathrm{d}t'.
  \end{align}
  $$
  

  - Absorb the delta-function into the integral (eES, pES...):
    $$
    \frac{N(t, E)}{\mathrm{d}t\mathrm{d}E} = \int_{E_\nu^\mathrm{min}}^{E_\nu^\mathrm{max}} \left[ F(E_\nu, t-\Delta t(E_\nu))
    \times \frac{\mathrm{d}\sigma}{\mathrm{d}T}(E_\nu, T(E)) \times g'(E)\right] \mathrm{d}E_\nu .
    $$

  - For simplified IBD case, no neutron recoiling considered ($E_\nu \approx E+1.8$ MeV):
    $$
    \frac{N(t, E)}{\mathrm{d}t\mathrm{d}E} =  F(E+1.8\,\mathrm{MeV}, t-\Delta t(E+1.8\,\mathrm{MeV}))
    \times \sigma(E+1.8\,\mathrm{MeV}) \times g'(E).
    $$
    

    

