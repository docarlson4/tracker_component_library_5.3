# Land vs Sea Clutter
---

## The Core Issue: Distribution Shape Depends on Resolution

In land clutter modeling, you typically treat σ° as a fixed (possibly random) backscatter coefficient and scale by cell area — resolution affects the *mean clutter power* but not the *distributional shape*. Sea clutter breaks this assumption entirely.

The K-distribution shape parameter ν is physically tied to the **effective number of independent scatterers per resolution cell**. Roughly:

$$\nu \propto \lambda_s \cdot A_c = \lambda_s \cdot R \cdot \frac{c\tau}{2} \cdot \theta_{az}$$

where λₛ is scatterer surface density. Consequences:

- **High resolution (small cell):** few scatterers → small ν → heavy-tailed, spiky clutter
- **Low resolution (large cell):** many scatterers → large ν → CLT applies → approaches Rayleigh
- ν **must be anchored to your radar parameters** — you cannot treat it as a free distributional choice independent of geometry

This is the critical inadequacy: land clutter lets you decouple the distribution from resolution; sea clutter does not.

---

## Two-Timescale Temporal Structure

Sea clutter is a compound process with explicitly different dynamics on two scales:

| Component | Physical origin | Timescale |
|-----------|----------------|-----------|
| **Speckle** | Microstructure interference | ~ms (∝ 1/Δf_clutter) |
| **Texture** | Gravity wave modulation | ~seconds (wave period) |

Your land clutter model likely has a single slow decorrelation timescale (scan-to-scan). Sea clutter requires both: fast AR speckle *modulated by* a slowly evolving gamma texture — exactly the structure in your Ward/Tough/Watts and Rosenberg implementations.

The speckle correlation time is approximately:

$$\tau_{speckle} \approx \frac{\lambda}{2\pi \sigma_v}$$

where σᵥ is the velocity spread of scatterers — which itself depends on sea state, wind, and grazing angle.

---

## Spatial Correlation Across Range Gates

Land clutter texture is tied to terrain features — spatially structured but essentially frozen. Sea clutter texture is a *propagating* spatial process (wave field), so:

- Adjacent range gates have correlated texture driven by wave coherence length
- The texture correlation function in range maps to a physical distance on the sea surface
- A CPI long enough to span multiple wave periods sees non-stationary texture

This is why your sea clutter implementations use FIR filtering of gamma variates per range gate with a specified correlation — you're synthesizing a spatially correlated texture field, not independent draws per cell.

---

## Grazing Angle Changes the Distribution Shape, Not Just the Mean

In land clutter, σ°(ψ) parameterizes the mean backscatter but the distributional family is roughly stable. In sea clutter:

- At **low grazing angles**, the spike mechanism (breaking waves, Bragg resonance) dominates — ν collapses toward 0.1–0.5, extremely non-Rayleigh
- At **moderate grazing**, ν ~ 1–5, moderately heavy-tailed
- At **high grazing (near vertical)**, many scatterers per cell, approaches Rayleigh

So grazing angle affects both the mean *and* the shape parameter simultaneously.

---

## Doppler Structure

Sea clutter has a **non-zero mean Doppler** (wave orbital velocity projected onto LOS) plus a spread tied to the wave velocity spectrum. This matters for:

- CFAR design in the Doppler domain
- STAP — sea clutter is not at DC the way land clutter approximates
- MTI cancellation — a notch centered at zero Doppler will *not* null sea clutter effectively

---

## Summary of What Your Land Models Lack

| Property | Land clutter (typical) | Sea clutter |
|----------|----------------------|-------------|
| Resolution dependence | Only mean power | *Shape parameter* ν and mean |
| Temporal structure | Single slow timescale | Two-timescale (speckle + texture) |
| Spatial correlation | Static terrain-driven | Propagating wave field |
| Grazing angle effect | Mean σ° only | Both mean and ν |
| Doppler | Near-zero | Non-zero mean + spread |
| Stationarity | Quasi-static | Non-stationary within CPI |

The compound model you've already implemented (gamma texture × Rayleigh speckle) is structurally correct for sea clutter — the key is that ν must be treated as a function of radar parameters (range, pulsewidth, beamwidth, grazing) rather than a fixed distributional assumption made independently of the simulation geometry.

---

# 

## Velocity Spread in Sea Clutter — A Deep Dive

---

### Physical Origins

The clutter Doppler spectrum is not a single line but a spread distribution arising from several distinct mechanisms that superpose:

**1. Orbital wave motion**
Surface gravity waves produce elliptical particle orbits. The surface velocity has both horizontal and vertical components. The radar sees the projection onto the LOS:

$$v_{orb}(x,t) = a\omega_w \cos(k_w x - \omega_w t) \cos\psi$$

where $a$ is wave amplitude, $\omega_w$ the wave angular frequency, $k_w$ the wavenumber, and $\psi$ the grazing angle. This gives a *mean* Doppler shift (not zero) plus a modulation — but over a resolution cell containing many waves, this contributes spread.

**2. The wave dispersion relation**
Ocean waves are not monochromatic. The Pierson-Moskowitz or JONSWAP spectrum gives a distribution of wavenumbers $k$ with dispersion:

$$\omega(k) = \sqrt{gk \tanh(kd)}$$

(deep water: $\omega = \sqrt{gk}$, so phase velocity $c_p = \sqrt{g/k}$). Each spectral component contributes a slightly different Doppler shift to the clutter spectrum. The width of the ocean wave spectrum maps *directly* onto the clutter Doppler spread.

**3. Bragg scattering — the dominant mechanism at low grazing**
At low grazing angles, the dominant backscatter mechanism is resonant Bragg scattering from short capillary/gravity waves whose wavelength satisfies:

$$\lambda_{Bragg} = \frac{\lambda_{radar}}{2\cos\psi}$$

For X-band (3 cm) at grazing 5°: $\lambda_{Bragg} \approx 1.5$ cm — very short waves that ride on top of long swell. These Bragg waves have:

$$v_{Bragg} = \pm c_p(\lambda_{Bragg}) = \pm\sqrt{\frac{g\lambda_{Bragg}}{2\pi}}$$

giving two spectral lines (forward/receding). But these short waves are *advected* by the underlying long wave orbital velocity field, which is itself a random process. The result: each Bragg line is *broadened* by the velocity variance of the long wave surface current. This is the Valenzuela/Plant mechanism and it's the primary reason $\sigma_v$ is large relative to what you'd estimate from the Bragg phase velocity alone.

**4. Breaking waves and whitecaps**
Non-Bragg scattering from breaking crests produces broadband, non-stationary Doppler contributions with high velocity spread and large amplitude spikes — the sea spikes that drive the heavy tail of the K-distribution. These are episodic, highly non-Gaussian in both amplitude and velocity, and essentially impossible to model deterministically.

**5. Wind-driven drift**
The surface Stokes drift and wind-driven current add a quasi-DC offset to the mean Doppler, but also contribute spread through their spatial variability across the resolution cell.

---

### Spectral Model

The clutter Doppler PSD is commonly modeled as Gaussian in velocity:

$$S(v) = \frac{\sigma_c^2}{\sqrt{2\pi}\sigma_v} \exp\left(-\frac{(v - \bar{v})^2}{2\sigma_v^2}\right)$$

with mean velocity $\bar{v}$ (nonzero — typically 0.3–1.5 m/s toward/away from radar depending on wind direction and sea state) and spread $\sigma_v$. In frequency:

$$S(f) = \frac{\sigma_c^2}{\sqrt{2\pi}\sigma_f} \exp\left(-\frac{(f - \bar{f})^2}{2\sigma_f^2}\right), \quad \sigma_f = \frac{2\sigma_v}{\lambda}$$

Representative values for X-band:

| Sea state | $\sigma_v$ (m/s) | $\bar{v}$ (m/s) |
|-----------|-----------------|-----------------|
| SS 1–2 | 0.1–0.3 | 0.2–0.5 |
| SS 3–4 | 0.3–0.7 | 0.5–1.0 |
| SS 5–6 | 0.7–1.5 | 1.0–2.0 |

At L-band, $\sigma_v$ is similar in m/s but maps to a *much narrower* fractional Doppler bandwidth — this is why lower frequencies are more tolerant of sea clutter for GMTI.

---

### Speckle Decorrelation Time Revisited

The temporal autocorrelation of the speckle component is the Fourier transform of $S(f)$ — for a Gaussian spectrum:

$$\rho_{speckle}(\tau) = \exp\left(-2\pi^2 \sigma_f^2 \tau^2\right) = \exp\left(-\frac{8\pi^2 \sigma_v^2 \tau^2}{\lambda^2}\right)$$

The 1/e decorrelation time:

$$\tau_{speckle} = \frac{\lambda}{2\sqrt{2}\pi \sigma_v}$$

Plugging in X-band ($\lambda = 3$ cm), $\sigma_v = 0.5$ m/s:

$$\tau_{speckle} \approx \frac{0.03}{2\sqrt{2}\pi \times 0.5} \approx 6.7 \text{ ms}$$

At PRF = 1 kHz, that's roughly **7 pulses** — meaning sea clutter speckle decorrelates *within a coherent processing interval*. This has direct consequences:

- You cannot assume clutter is fully correlated across the CPI (as land clutter often is modeled)
- Nor is it fully decorrelated pulse-to-pulse (as thermal noise is)
- The partial correlation structure shapes the clutter covariance matrix $\mathbf{R}_c$ that STAP and adaptive beamformers must estimate and invert

---

### Implications for the Clutter Covariance Matrix

The space-time clutter covariance (for a ULA with $N$ elements, $M$ pulses) is not rank-1 as it would be for a single point scatterer. The temporal covariance block for a single range gate is:

$$[\mathbf{R}_{temp}]_{mn} = \sigma_c^2 \cdot \rho_{tex}(|m-n|T_{PRI}) \cdot \rho_{spk}(|m-n|T_{PRI})$$

where $\rho_{tex}$ is the slow gamma texture correlation and $\rho_{spk}$ is the fast Gaussian-spectrum speckle correlation above. The product of two distinct timescale processes gives a covariance that is *not* described by any single AR model — this is the fundamental difficulty for clutter whitening.

The effective clutter rank (Brennan's rule generalization) is inflated relative to land clutter because the Doppler spread means clutter energy is not confined to the clutter ridge but smeared around it in angle-Doppler space.

---

### Connection to Your Implementations

In your Ward/Rosenberg MATLAB implementations, the FIR filter applied to the complex Gaussian speckle channel has a frequency response shaped to match $S(f)$ above. The filter bandwidth (in normalized frequency) is:

$$B_{norm} = \frac{2\sigma_v}{\lambda \cdot PRF}$$

This is where $\sigma_v$, $\lambda$, and PRF must be mutually consistent — if you synthesize a data cube with a particular $\sigma_v$ but then run STAP or MTI designed around a different assumed spread, the SINR predictions will be off. It's worth verifying your filter design explicitly encodes the physical $\sigma_v$ rather than an ad-hoc normalized bandwidth.

The mean Doppler shift $\bar{v}$ is equally important: it should offset the speckle filter's center frequency away from DC by $\bar{f} = 2\bar{v}/\lambda$, otherwise your simulated sea clutter sits at zero Doppler and a naive MTI canceler will suppress it — producing optimistic false-alarm performance that won't hold against real sea clutter data.