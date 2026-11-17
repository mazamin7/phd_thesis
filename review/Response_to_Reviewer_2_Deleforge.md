# Response to Reviewer 2 (Dr. Antoine Deleforge)

**PhD Thesis:** *Domain Decomposition and Differentiable Physics for the Viscoelastic Wave Equation*  
**Candidate:** Gerardo Cicalese  
**Institution:** Politecnico di Milano  
**Reviewer:** Dr. Antoine Deleforge (Senior Researcher, Inria Rennes – Bretagne Atlantique)

---

Dear Dr. Antoine Deleforge,

We would like to express our sincere gratitude for the time, rigor, and depth with which you reviewed this doctoral thesis. Your insightful and constructive feedback has been immensely valuable in identifying critical areas where the exposition, mathematical precision, experimental analysis, and balance of claims required substantial improvement.

In response to your review, we have thoroughly revised the thesis. All modifications directly addressing your comments have been highlighted in **red** (`\textcolor{red}{...}`) in the updated manuscript. 

Below, we first provide a comprehensive response to your **General Assessment**, outlining the overarching changes made across the thesis. Following that, we address every single **Major Comment (R2C1 to R2C17)** and **Minor Comment / Typo (R2m1 to R2m11)** in detail.

---

# Part 1: Response to the General Assessment

### 1. Scope, Applicative Impact, and 3D Experimental Validation
> **Reviewer's remark:**  
> *"From the perspective of the audio and acoustic signal processing community, I think this work is relevant. (...) However, the fact that only one qualitative experiment in 3D is presented limits the thesis' potential applicative impact."*

**Response:**  
We fully agree with your observation. A single qualitative 3D snapshot was insufficient to demonstrate the practical viability and numerical accuracy of the multi-subdomain ARD-PSTD framework. To address this limitation, we substantially expanded Section 3.4.5 (previously Sec. 3.2.5) into a comprehensive, multi-solver comparative 3D benchmark:
1. **Implementation of Two Independent 3D Wave-Based Solvers:**  
   We developed and implemented two independent global wave-based reference solvers operating on the exact same stepped non-convex 3D room geometry:
   - **Global 3D FDTD:** A second-order standard leapfrog Cartesian finite-difference scheme.
   - **3D Acoustic Time-Domain FEM (TD-FEM):** An explicit variational solver using 8-node trilinear hexahedral ($Q_1$) elements with analytical 27-point element coupling and row-sum mass lumping.
2. **Multi-Resolution Convergence Study:**  
   Both FDTD and TD-FEM were executed on two spatial grid discretizations: a coarse mesh ($\Delta h = 0.20$\,m, matching coarse ARD-PSTD) and a refined mesh ($\Delta h = 0.10$\,m, halving the grid spacing and requiring $8\times$ more spatial cells and $2\times$ more time steps).
3. **Quantitative Frequency-Domain Validation (Room FRF/RTF):**  
   We extracted the acoustic pressure time histories at an interior receiver $\mathbf{x}_r = (15.0, 15.0, 3.0)$\,m inside Subdomain 2 and computed the Room Frequency Response Function ($|H(f)|$) via regularized deconvolution up to the spatial Nyquist limit ($858.75$\,Hz). As reported in Figure 3.15, coarse ARD-PSTD exhibits excellent modal peak agreement with the fine-grid FDTD and fine-grid TD-FEM solvers below $200$\,Hz. We also explicitly delineate that above $200$\,Hz, in the absence of an analytical closed-form solution or physical measurements, cross-solver discrepancies establish the determinable limit of accuracy.
4. **Empirical Multi-Core Computational Benchmarking:**  
   We reported measured execution times on an 8-core CPU (AVX2 vectorization), showing that multi-threaded ARD-PSTD ($36.0$\,s) achieves a $2.57\times$ speedup over fine FDTD ($92.6$\,s) and a $4.0\times$ speedup over fine TD-FEM ($144.3$\,s). We clearly distinguished measured execution times from theoretical cell-reduction ratios.

---

### 2. Analysis, Insights, and Commentary on Numerical Results
> **Reviewer's remark:**  
> *"Numerical results on the 6 benchmarks are presented. Unfortunately, there is very little comments or analysis of the results by the author, beyond presenting them. (...) For the two damped scenarios (the most novel contributions), there is in fact no comments at all by the author, which is a problematic omission. No conclusion is offered for this chapter."*

**Response:**  
We acknowledge this critical omission in the earlier version. Section 3.2 has been enriched with dedicated, in-depth analytical discussions for all benchmark experiments:
- **Experiment 3 (Smooth Damped Standing Wave):** We added an explicit analysis explaining why PSTD maintains machine-precision accuracy ($10^{-15} \div 10^{-14}$): because the initial field is confined to mode $m=1$, the state transition matrix $S_m(\Delta t) = e^{A_m \Delta t}$ evaluates the exact exponential flow of both telegrapher and viscoelastic damping without any time-stepping error.
- **Experiment 4 (Non-Smooth Damped Standing Wave):** We provided a detailed physical explanation of the regularizing role of viscoelastic damping. In contrast to telegrapher damping (which decays all modes uniformly), viscoelastic damping strongly attenuates high-frequency content in the underdamped regime, progressively smoothing slope discontinuities. We showed that this selective dissipation extinguishes the singular high-frequency components that otherwise degrade finite-difference stencils, enabling FDTD to recover an observed second-order convergence rate ($\approx 2.22$, down from first-order $\approx 1.03$ in the undamped case).
- **Execution Times in Benchmarks:** We reported empirical runtimes (e.g. in Experiment 6, PSTD achieves an error reduction factor of $42.2\times$ over FDTD).
- **Chapter Synthesis:** We added Section 3.2.1 (*Comparative Synthesis and Critical Discussion*), synthesizing the quantitative trade-offs between dispersion-free modal transport, DCT algebraic complexity, and geometric flexibility.

---

### 3. Chapter 4 Conclusion and Balance of Scientific Claims
> **Reviewer's remark:**  
> *"Chapter 4 offers a one-page conclusion for the 100-page thesis. (...) the conclusion only praises the work, sometimes overly so, and makes a number of exaggerated claims without mentioning any limitation. Three possible directions for future work are outlined in one sentence (...). Finally, the claim 'achieve unprecedented computational efficiency and spatial accuracy' has not been validated (...). Likewise, the claim 'significantly outperforming purely data-driven black-box inversions' is made with no evidence."*

**Response:**  
We completely reworked Chapter 4 into a rigorous, mature, and balanced scientific evaluation:
1. **Summary of Main Contributions (Section 4.1):**  
   All superlative and unbacked claims have been excised. The text now presents a measured, defensible synthesis that delineates exact conditions under which spectral accuracy and efficiency advantages hold.
2. **Dedicated Limitations Section (Section 4.2):**  
   We added a complete, transparent section (*Limitations of the Present Study*) articulating the boundaries of our findings across four pillars:
   - *Discretization and Splitting Order:* Acknowledging the reduction of Strang splitting to first-order in time when using forward Euler Kick on damped interfaces ($\nu > 0$), and the sublinear/first-order spatial convergence of explicit interface stencils.
   - *Geometric and Physical Restrictions:* Axis-aligned orthogonal cuboid subdomains, staircasing on non-conforming boundaries, and uniform rigid boundaries versus complex frequency-dependent surface impedance.
   - *Domain Decomposition Scalability:* Sensitivity of OSWR to frequency truncation and modest 1D acceleration versus high-damping 2D gains.
   - *Inverse Problems and Differentiable Learning:* Synthetic 1D noise-free/moderate-noise data, known zero initial velocity, potential coupling between $\nu$ and the Sobolev weight, and the lack of architecture ablation.
3. **Comprehensive Future Research Directions (Section 4.3):**  
   Expanded from 3 lines to over 4 pages across 5 structured subsections:
   - *Section 4.3.1:* Geometric discretization and staircasing boundary error mitigation (cut-cell, embedded boundary, curvilinear transforms).
   - *Section 4.3.2:* Energy stability and higher-order time integration for ARD (exact exponential Kick flow, explicit midpoint, discrete energy estimates, high-order DtN Padé operators).
   - *Section 4.3.3:* Analytical min-max theory and adaptive transmission for OSWR.
   - *Section 4.3.4:* Scalability, identifiability, and architectural ablation for differentiable learning (3D scaling, sensor sparsity, ablation of SIREN vs Fourier features).
   - *Section 4.3.5:* Realistic acoustic modelling, advanced boundary conditions (CPML/CFS-PML with auxiliary memory variables, digital impedance filters), and GPU acceleration.
4. **Retraction of Unsubstantiated Claims:**  
   - Retracted *"achieve unprecedented computational efficiency and spatial accuracy"*.
   - Retracted *"significantly outperforming purely data-driven black-box inversions"*.
   - Retracted the claim that ARD was extended to FDTD.

---

# Part 2: Point-by-Point Response to Major Comments

### Comment R2C1: Definition of Specialized Technical Terms
> **Reviewer's remark:**  
> *There are many occurrences of highly technical terms that are used in the thesis without being properly defined, introduced or backed by any reference. Here are a few examples (...):*  
> *- P9: "global support", "computational domain", "interpolant"*  
> *- P9: "spectral convergence", "exponential convergence", "numerical dispersion"*  
> *- P10: "Runge-Kutta scheme"*  
> *- P11: "transmission conditions"*  
> *- P38: "the amplification polynomial"*  
> *- P87: "explicit leapfrog FDTD" and "first-order Velocity Verlet scheme"*  
> *- P90: "interface stencils"*  
> *- P100: "using Cosine Annealing schedules"*

**Response & Action Taken:**  
We conducted a comprehensive audit across the thesis and formally defined and referenced all flagged terms (along with additional technical terms in Chapters 1–3):
- **P9 (*global support*, *computational domain*, *interpolant*):** Added explicit definitions in Section 1.3 explaining that in global spectral methods, basis functions (such as Chebyshev polynomials or trigonometric Fourier modes) are non-zero across the entire spatial domain $\Omega$, in contrast to local compactly supported finite elements or finite-difference stencils \citep{canuto_spectral_2006,trefethen_spectral_2000}.
- **P9 (*spectral/exponential convergence*, *numerical dispersion*):** Formally defined *spectral (exponential) convergence* as the property where approximation errors decay faster than any algebraic power $\mathcal{O}(N^{-k})$, scaling as $\mathcal{O}(e^{-c N})$ for analytic functions \citep{boyd_chebyshev_2001}. Defined *numerical dispersion* as the unphysical dependence of discrete phase velocity on spatial frequency, grid spacing, and propagation direction \citep{bilbao_numerical_2009}.
- **P10 (*Runge-Kutta scheme*):** Formally introduced multi-stage explicit and implicit Runge-Kutta time-stepping integrators, citing standard references \citep{butcher_numerical_2016,hairer_solving_1993}.
- **P11 (*transmission conditions*):** Defined interface transmission conditions in domain decomposition as boundary operators (Dirichlet, Neumann, or Robin-type) enforced across artificial sub-interfaces to ensure continuity of state and fluxes across adjacent subdomains \citep{quarteroni_domain_1999,gander_optimized_2007}.
- **P38 (*the amplification polynomial*):** Clarified that the *amplification polynomial* is the characteristic polynomial of the discrete time-evolution operator whose roots define the spectral amplification factor $G(k,\Delta t)$ governing numerical stability.
- **P87 (*explicit leapfrog FDTD* and *Velocity Verlet scheme*):** Formally defined the centered leapfrog time-stepping update for second-order hyperbolic systems and demonstrated its exact algebraic isomorphism to the staggered-velocity Velocity Verlet formulation \citep{swope1982computer,bilbao_numerical_2009}.
- **P90 (*interface stencils*):** Defined interface stencils as localized finite-difference derivative operators bridging boundary-adjacent nodes across artificial domain interfaces to evaluate the residual dipole operator $C_{\text{res}} = L_h - L_h^{(\mathrm{N})}$.
- **P100 (*Cosine Annealing schedules*):** Formally defined Cosine Annealing in Section 3.5, providing the analytical equation:
  $$\eta_k(t) = \eta_{k,\min} + \frac{1}{2}(\eta_{k,\max}-\eta_{k,\min})\left[1 + \cos\left(\pi\frac{t}{T_k}\right)\right]$$
  and citing \citet{loshchilov2017sgdr}.

---

### Comment R2C2: Numerical Results on Damped Benchmarks (Pages 66–67)
> **Reviewer's remark:**  
> *The damped numerical results on pages 66 and 67 are shown without any comments by the author. This is a serious issue that should be fixed. More generally, (...) the author does not provide sufficient explanations, insights or critical views on the obtained results, which are generally just factually presented.*

**Response & Action Taken:**  
We expanded the discussion of Experiment 3 (Smooth Damped) and Experiment 4 (Non-Smooth Damped) in Section 3.2. As detailed in our response to General Assessment #2, we explained the exact analytical modal integration responsible for machine precision in Exp 3, and provided a rigorous physical analysis of the regularizing action of viscoelastic dissipation on the corner singularities in Exp 4, explaining why FDTD recovers second-order convergence ($\approx 2.22$) only when viscoelastic damping is active.

---

### Comment R2C3: Motivation for Two Overlapping Subdomains in Sec. 3.1.4
> **Reviewer's remark:**  
> *Sec. 3.1.4. "the domain is partitioned along the x-axis into two overlapping subdomains": Why is it partitioned into two overlapping subdomain? This is never motivated or justified.*

**Response & Action Taken:**  
In Section 3.3.1 (Sec. 3.1.1 in Chapter 3), we clarified that partitioning the domain into two overlapping subdomains $\Omega_1 = (0, b)$ and $\Omega_2 = (a, L)$ with overlap $\delta = b - a > 0$ is the foundational, canonical configuration in Schwarz domain decomposition theory \citep{quarteroni_domain_1999,gander_optimized_2007,dolean_introduction_2015}. It enables the analytical derivation of the continuous convergence factor $\rho(\omega; q, r)$ via Fourier/Laplace transforms. We explained that the spatial overlap $\delta > 0$ provides geometric contraction, ensuring uniform error decay across iterations prior to multi-subdomain extension.

---

### Comment R2C4: Retraction of Claim that FDTD is Presented for ARD
> **Reviewer's remark:**  
> *Sec. 3.2: "Although the FDTD formulation is presented for completeness". It is not really presented. Rather, the author writes in the corresponding section: "we omit the explicit FDTD particularization here and proceed directly to the pseudo-spectral discretization". Please retract this claim.*

**Response & Action Taken:**  
We retracted this claim. The sentence was amended to state clearly:  
*"We omit the explicit FDTD particularization here and proceed directly to the pseudo-spectral discretization, as standard explicit FDTD relies on strictly local stencil operations that do not require domain decomposition on Cartesian grids."*

---

### Comment R2C5: Justification of Parabolic Stability Limit
> **Reviewer's remark:**  
> *P86: "the explicit evaluation of the viscous term $\nu C_{\text{res}} v$ introduces a strict parabolic stability limit ($\Delta t \lesssim \Delta x^2/\nu$)." This is not sufficiently detailed or justified.*

**Response & Action Taken:**  
In Section 3.4.2 (lines 4431–4442 in `Thesis.tex`), we added a rigorous derivation based on matrix perturbation and **Gershgorin's circle theorem**:
1. During the Kick substep, the homogeneous velocity perturbation satisfies $\delta \mathbf{v}(t+\tau) = (I + \tau \nu C_{\text{res}}) \delta \mathbf{v}(t)$.
2. The residual matrix $C_{\text{res}} = L_h - L_h^{(\mathrm{N})}$ is a spatial difference operator scaling as $\mathcal{O}(\Delta x^{-2})$. By Gershgorin's circle theorem \citep[Chapter~5]{quarteroni_numerical_2007}, its spectral radius is bounded by $\rho(C_{\text{res}}) \le C_p / \Delta x^2$ (with $C_2 = 4$ for the three-point Laplacian).
3. Assuming non-positive real eigenvalues $\lambda_j(C_{\text{res}}) \le 0$ for dissipative stability, the forward Euler update remains stable if $|1 + \tau \nu \lambda_j| \le 1$. With $\tau = \Delta t / 2$, this requires:
   $$1 - \frac{\Delta t}{2}\nu \rho(C_{\text{res}}) \ge -1 \implies \Delta t \le \frac{4}{\nu \rho(C_{\text{res}})} \implies \Delta t \le \frac{4\Delta x^2}{C_p \nu}.$$
This explicitly justifies the parabolic stability restriction $\Delta t \lesssim \Delta x^2 / \nu$.

---

### Comment R2C6: Initialization of the Iterative ARD Algorithm
> **Reviewer's remark:**  
> *Sec. 3.2.4.: How is the presented iterative algorithm initialised, i.e., what are the initial values of $U_m^n$ and $V_m^n$?*

**Response & Action Taken:**  
In Section 3.4.4, we explicitly clarified the initialization procedure:
- At $t = 0$, the physical initial conditions $u_0(x)$ and $v_0(x)$ are sampled on the grid nodes of each subdomain $\Omega_k$.
- The initial modal coefficients are obtained by evaluating the forward Discrete Cosine Transform (DCT-II):
  $$\mathbf{U}_k^0 = \operatorname{DCT-II}(\mathbf{u}_{k,0}), \qquad \mathbf{V}_k^0 = \operatorname{DCT-II}(\mathbf{v}_{k,0}).$$
- The interface residue exchange arrays are initialized to zero (or evaluated from the initial state), ensuring that the first half-Kick step $K(\Delta t/2)$ applies the correct instantaneous interface forces before entering the time-stepping loop.

---

### Comment R2C7: Softening "Validating 3D ARD" and Adding Quantitative Evaluation
> **Reviewer's remark:**  
> *P96: "validating the accuracy of the multi-domain pseudo-spectral time-domain (PSTD) formulation in 3D." This is too strong a claim, I would not call this "validating". A single configuration is tested, with two (very close) values of the viscoelastic damping parameter. Only qualitative results with three snapshots per simulations are displayed, with no quantitative evaluation or comparison to other solvers. This is very weak experimentally.*

**Response & Action Taken:**  
We replaced the word *"validating"* with *"numerical demonstration"* and completely upgraded the experimental methodology:
1. Implemented **Global 3D FDTD** and **3D Acoustic TD-FEM** solvers on the same stepped geometry.
2. Evaluated all solvers on both coarse ($\Delta h = 0.20$\,m) and fine ($\Delta h = 0.10$\,m) meshes.
3. Added quantitative frequency-domain cross-verification (Figure 3.15) up to the spatial Nyquist limit ($858$\,Hz), establishing modal agreement below $200$\,Hz.
4. Explicitly acknowledged the analysis boundaries above $200$\,Hz due to the absence of analytical ground truth.

---

### Comment R2C8: Internal Boundary Handling in 3D ARD Snapshots
> **Reviewer's remark:**  
> *Fig. 3.12: "The solver seamlessly handles the internal boundaries between the three rectangular subdomains." One cannot actually see that from the figure.*

**Response & Action Taken:**  
Figure 3.12 (now Figure 3.14) was replaced with high-resolution wavefield snapshots comparing all five solver configurations across three epochs ($t = 0.03, 0.08, 0.16$\,s). We overlaid distinct red dashed lines indicating the exact locations of the virtual sub-interfaces $\Gamma_1$ and $\Gamma_2$. In the accompanying text, we explicitly discuss the continuity of the wavefield across these interfaces, noting the absence of spurious boundary reflections or phase disruptions in the displayed horizontal slices.

---

### Comment R2C9: Spatiotemporal Absolute Error Plotting and 3D Damping Visualization
> **Reviewer's remark:**  
> *Fig. 3.16 is not very informative. It is impossible to see differences between the 3 images. Plotting the error or relative error would be better.*

**Response & Action Taken:**  
We completely redesigned the visual evaluation following your recommendation:
1. **Spatiotemporal Absolute Error Heatmaps in DPNN (Figures 3.16, 3.19, and 3.22):** Rather than plotting raw wavefield snapshots that look nearly identical, each spatiotemporal figure now presents a four-panel layout on a unified color scale:
   - Panel (a): Ground-truth wavefield ^*(x,t)$ (seismic bipolar colormap).
   - Panel (b): FDTD-DPNN absolute error $|u_{\mathrm{FDTD}}(x,t) - u^*(x,t)|$.
   - Panel (c): PSTD-DPNN absolute error $|u_{\mathrm{PSTD}}(x,t) - u^*(x,t)|$.
   - Panel (d): Standard PINN baseline absolute error $|u_{\mathrm{PINN}}(x,t) - u^*(x,t)|$.
   This error visualization makes the numerical properties immediately transparent: PSTD-DPNN maintains near-zero error across the entire domain, FDTD-DPNN clearly reveals trailing wave-like error ripples caused by numerical grid dispersion ({\mathrm{num}} < c$), and PINN exhibits large diffuse errors and phase drift throughout the interior.
2. **3D Damping Comparisons (Figure 3.15):** We updated the 3D ARD snapshots with standardized colorbars and added an explicit physical scaling analysis explaining why primary wavefront amplitudes decay modestly over  = 0.24$\,s for nominal atmospheric damping ($
u pprox 1.2\cdot 10^{-5}$\,m2$/s), while the reverberant coda undergoes progressive high-frequency attenuation.

---

### Comment R2C10: Explicit Description of Cosine Annealing in DPNN
> **Reviewer's remark:**  
> *The loss plots (Fig. 3.14, 3.17, 3.20...) feature abrupt changes, suggesting that some elaborate techniques are used to modify the learning rate along the way, but there is no explicit mention of this despite a mention of "cosine annealing" which is not described at all.*

**Response & Action Taken:**  
In Section 3.5 (lines 5100–5113), we added the complete mathematical formulation of the Cosine Annealing learning rate schedule \citep{loshchilov2017sgdr} and explained the dual-optimizer strategy:
- Network weights $\Theta$ are updated with a conservative learning rate $\eta \in [10^{-5}, 10^{-3}]$ to preserve representation stability.
- Unconstrained physical damping variables $(\tilde{\gamma}, \tilde{\nu})$ are assigned an aggressive schedule starting at $\eta = 10^{-1}$ (to rapidly traverse the non-convex loss landscape and escape local minima) and annealed down to $10^{-3}$ (to prevent limit cycles and ensure accurate asymptotic locking onto ground-truth parameters).

---

### Comment R2C11: Critical Scope, Quantitative Benchmarks, and PINN Double Comparison
> **Reviewer's remark:**  
> *Sec 3.3.1: As in chapter 3, the comments on the experiments are extremely limited. This section also severely lacks quantitative analysis. There are 3 qualitative experiments, which seems too limited to rigorously validate the approach. Especially since quite strong limitations seem to appear in the third experiment.*

**Response & Action Taken:**  
In the revised thesis, we significantly upgraded Section 3.5.1 by adding a comprehensive quantitative benchmark suite, introducing a comparative baseline, and providing in-depth physical analysis:
1. **Quantitative Benchmark Table (Table 3.4):** We integrated a complete numerical comparison reporting four quantitative error metrics across all three initial condition benchmarks (Multimodal Standing Waves, Smooth Gaussian Pulses, Non-Smooth Triangular Pulses):
   - Normalized Mean Squared Error on initial pressure: $\mathrm{NMSE}(p_0) = \|p_0^{\mathrm{pred}} - p_0^*\|_2^2 / \|p_0^*\|_2^2$.
   - Full spatiotemporal wavefield relative $ error: $\|u_{\mathrm{pred}} - u^*\|_2 / \|u^*\|_2$ (%).
   - Relative parameter estimation errors for telegrapher damping $|\gamma - \gamma^*|/\gamma^*$ and viscoelastic damping $|
u - 
u^*|/
u^*$.
   - Computational wall-clock training time per epoch (ms/epoch).
2. **Double Comparison against Standard PINN Baseline:** We implemented a standard Physics-Informed Neural Network (PINN) collocation baseline (Raissi et al., 2019) optimizing the joint loss $\mathcal{L}_{\mathrm{PINN}} = \lambda_{\mathrm{data}}\mathcal{L}_{\mathrm{data}} + \lambda_{\mathrm{pde}}\mathcal{L}_{\mathrm{pde}} + \lambda_{\mathrm{bc}}\mathcal{L}_{\mathrm{bc}} + \lambda_{\mathrm{ic}}\mathcal{L}_{\mathrm{ic}}$. This double comparison (PSTD-DPNN vs. FDTD-DPNN vs. PINN) demonstrates that:
   - *PSTD-DPNN achieves superior accuracy*: $\mathrm{NMSE}(p_0) = 2.22 	imes 10^{-2}$ on standing waves, outperforming FDTD-DPNN (.19 	imes 10^{-2}$) and dramatically outperforming PINN (.00$).
   - *PSTD-DPNN is .8	imes$ faster than PINN*: PSTD-DPNN requires 8.6\,ms/epoch versus 75.5\,ms/epoch for PINN. Computing the stiff third-order mixed derivative $-
u \partial_{txx} u$ through automatic differentiation creates large computational overhead and severe spectral bias in PINNs, whereas DPNN integrates the differential operators algebraically in modal space.
3. **In-Depth Analysis of Limitations on Triangular Pulse (Experiment 3):** We expanded the discussion of Experiment 3 to explicitly explain why relative errors are higher on non-smooth profiles:
   - Sharp corners possess algebraic modal decay ((m^{-2})$), and viscoelastic damping strongly suppresses high frequencies (- \sim -
u k_m^2$). Under 10 sparse sensors with 30\,dB SNR noise, these attenuated components are obscured by noise before reaching the receivers.
   - An unconstrained inversion would invert the dissipated modes and violently amplify measurement noise. The Sobolev ^1$ regularization deliberately sacrifices peak sharpness to prevent spurious hallucinations, trading corner fidelity for stability.

---

### Comment R2C12: Retraction of "Remarkable Robustness" on Triangular Pulse
> **Reviewer's remark:**  
> *P104: "Figures 3.20 to 3.22 showcase the remarkable robustness of the carefully calibrated Sobolev H1 regularization". To be honest, not really. The relative reconstruction error on the initial condition seems very large. This claim should be retracted or further expanded.*

**Response & Action Taken:**  
We retracted the phrase *"remarkable robustness"*. In Section 3.5.1 (Experiment 3, lines 5188–5193), we expanded the explanation:
- Acknowledged that visible relative reconstruction errors remain around sharp slope discontinuities.
- Explained that this is mathematically expected: viscoelastic attenuation rapidly extinguishes high-frequency modes before they reach the sensors; under noisy/sparse measurements, unregularized fitting would invert attenuated modes and amplify noise into spurious oscillations.
- Clarified that the $H^1$ Sobolev penalty acts as a deliberate trade-off, smoothing the sharp peaks to prevent noise-driven hallucinations while still allowing accurate recovery of macroscopic damping parameters.

---

### Comment R2C13: Toning Down "Unprecedented Computational Efficiency"
> **Reviewer's remark:**  
> *P107: "achieve unprecedented computational efficiency and spatial accuracy." There is no solid proof of this in the thesis, this claim should be tuned down or retracted.*

**Response & Action Taken:**  
We completely retracted this claim in both Chapter 3 and Chapter 4. It was replaced with measured multi-threaded benchmarks: coarse ARD-PSTD is $1.63\times$ slower than coarse FDTD per step due to 3D DCT-II transforms, but executes $2.57\times$ faster than fine FDTD and $4.0\times$ faster than fine TD-FEM while completely eliminating numerical dispersion for resolved frequencies.

---

### Comment R2C14: Retraction of Claim Extending ARD to FDTD
> **Reviewer's remark:**  
> *P107: "Subsequently, we extended the Adaptive Rectangular Decomposition (ARD) method to both FDTD and PSTD frameworks.": No, this is incorrect, the FDTD framework is not treated. The author mentions it in a paragraph, but then explicitly writes that it won't be treated.*

**Response & Action Taken:**  
We corrected this statement in Chapter 4 (Summary of Main Contributions) and Section 3.4. The text now explicitly states that ARD was developed, formulated, and implemented exclusively for the PSTD framework.

---

### Comment R2C15: Toning Down Superlative Claims on Sobolev Regularization
> **Reviewer's remark:**  
> *P108: "(...) Sobolev H_1 smoothness penalty gracefully suppressed high-frequency artifacts without smearing out sharp discontinuities": This is a highly subjective and superlative statement. In fact, there is no real rigorous treatment of this, just a succinct qualitative experiment.*

**Response & Action Taken:**  
We retracted *"gracefully suppressed (...) without smearing out sharp discontinuities"*. We rephrased this across the thesis to state that the $H^1$ gradient penalty penalizes modes proportionally to $k_m^2$, acting as a physically motivated regularizer that suppresses unconstrained high-frequency oscillations at the expense of smoothing sharp spatial features.

---

### Comment R2C16: Retraction of Unbacked Claim on Black-Box Inversions
> **Reviewer's remark:**  
> *P108: "(...), significantly outperforming purely data-driven black-box inversions": this claim has no supporting evidence and should be retracted. No comparison is provided for this method, nor any quantitative results.*

**Response & Action Taken:**  
We retracted this claim entirely. All unbacked comparative assertions regarding black-box data-driven methods have been removed from the thesis.

---

### Comment R2C17: Expansion of Future Research Directions
> **Reviewer's remark:**  
> *P108: 3 lines on future work is an order of magnitude less than what one would typically expect for a 3-year-long doctoral work.*

**Response & Action Taken:**  
As described in General Assessment #3, Section 4.3 was expanded from 3 lines to over 4 pages of continuous technical prose, detailing 5 distinct research avenues (geometric discretization, energy stability / Kick splitting, analytical OSWR min-max theory, 3D differentiable learning scaling / identifiability, and advanced boundary conditions / CPML / impedance filters).

---

# Part 3: Point-by-Point Response to Minor Comments and Typos

- **R2m1 (P4: Timeline of wave-based methods):**  
  *Clarification:* We clarified Section 1.1: geometrical ray tracing emerged around 1970, while early time-domain wave methods (such as Yee's 1966 FDTD) were introduced in electromagnetics in the late 1960s/1970s and subsequently adapted for room acoustics in the late 1990s.
- **R2m2 (P12: Structural placement of DPNN paragraphs):**  
  *Clarification:* Recontextualized the introductory paragraphs in Section 1.4 to clearly separate forward numerical solvers from the inverse problem formulation.
- **R2m3 (P18: Equation reference for $k(i\omega)$):**  
  *Fix:* Corrected the cross-reference to point directly to the continuous dispersion relation (Equation 2.1).
- **R2m4 (P24: Clarification of the three models):**  
  *Clarification:* Explicitly identified the three governing PDE models: (1) undamped wave equation ($\gamma = \nu = 0$), (2) telegrapher-damped wave equation ($\gamma > 0, \nu = 0$), and (3) combined telegrapher and viscoelastic wave equation ($\gamma > 0, \nu > 0$).
- **R2m5 (P40: Role of Proposition 2.4):**  
  *Clarification:* Added an explanatory remark stating that Proposition 2.4 explicitly formalizes the discrete transition matrix and stability bounds for the specific coupled damping operator studied.
- **R2m6 (P54: Clarification of "advanced in the modal domain"):**  
  *Clarification:* Rephrased to: *"The spatial field is projected onto discrete cosine basis functions via the DCT-II, advanced in time analytically via decoupled modal transition matrices, and synthesized back to physical space via the inverse DCT-II."*
- **R2m7 (P71: Grammatical structure fix in Section 3.1):**  
  *Fix:* Removed the redundant conjoining phrase (*"for which"*), restoring clean sentence syntax.
- **R2m8 (P81: Domain overlap $a < b$):**  
  *Fix:* Explicitly specified $0 < a < b < L$, defining the overlap width as $\delta = b - a > 0$.
- **R2m9 (P86: Citation typo):**  
  *Fix:* Corrected the duplicate author name in `(Strang Strang (1968))` to `(Strang, 1968)`.
- **R2m10 (P89: Citation typo):**  
  *Fix:* Replaced inline text with standard `\citep{strang_construction_1968}`.
- **R2m11 (P90: Definition of $D$ and $K$ operators):**  
  *Fix:* Formally defined the Drift operator $D(\Delta t) = \mathcal{S}_{\text{local}}(\Delta t)$ (unconstrained local subdomain solver) and the Kick operator $K(\tau) = I + \tau B$ (instantaneous interface velocity correction).

---

We believe that these revisions have substantially elevated the quality, balance, and rigor of the thesis, directly resolving every concern you raised. We once again thank you for your invaluable critique.

Sincerely,  
**Gerardo Cicalese**  
*Ph.D. Candidate, Politecnico di Milano*
