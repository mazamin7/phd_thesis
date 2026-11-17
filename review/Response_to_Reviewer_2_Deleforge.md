# Response to Reviewer 2 (Dr. Antoine Deleforge)

> **Archived working draft — do not submit.** This Markdown copy is stale and contains obsolete experimental values and claims. The current response is [Response_to_Reviewer_2_Deleforge.tex](Response_to_Reviewer_2_Deleforge.tex).

**PhD Thesis:** *Domain Decomposition and Differentiable Physics for the Viscoelastic Wave Equation*  
**Candidate:** Gerardo Cicalese  
**Advisor:** Prof. Ilario Mazzieri  
**Co-Advisors:** Prof. Gabriele Ciaramella, Prof. Fabio Antonacci  
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
We agree that presenting a single qualitative 3D snapshot comparison limited the practical evaluation of the proposed solvers. In the revised dissertation, we significantly expanded Section 3.2.6 (3D ARD) and Section 4.2 (Limitations):
1. **Comprehensive Multi-Solver 3D Benchmark:** We expanded the 3D room simulation ($L_x \times L_y \times L_z = 6.0\,\mathrm{m} \times 4.0\,\mathrm{m} \times 3.0\,\mathrm{m}$, $T = 0.24$\,s) into a systematic benchmark comparing five configurations: coarse FDTD ($h = 0.10$\,m), fine FDTD ($h = 0.05$\,m), coarse TD-FEM, fine TD-FEM, and our 3-subdomain ARD-PSTD solver.
2. **Quantitative Execution Runtimes:** We added measured 8-core CPU runtimes: ARD-PSTD requires 102.4\,s per run ($1.63\times$ slower than coarse FDTD per step due to 3D DCT-II transforms, but $2.57\times$ faster than fine FDTD at 263.8\,s, and $4.0\times$ faster than fine TD-FEM at 409.6\,s).
3. **Damping Parameter Scaling:** We added wavefield comparisons across four damping regimes: lossless ($\nu = 0$), low damping ($\nu \approx 2.4\cdot 10^{-6}$\,m$^2$/s), moderate damping ($\nu \approx 1.2\cdot 10^{-5}$\,m$^2$/s, corresponding to nominal atmospheric absorption at 20$^\circ$C and 50% relative humidity), and high damping ($\nu \approx 5.9\cdot 10^{-5}$\,m$^2$/s). Over $T = 0.24$\,s with $k_{\max} \approx 15.7$\,m$^{-1}$, primary wavefronts undergo modest decay while the reverberant coda experiences progressive high-frequency attenuation.
4. **Explicit Architectural Limitations:** In Section 4.2, we explicitly bounded the operational domain: the present 3D implementation is restricted to axis-aligned rectangular subdomains, rigid boundaries, and nominal audio frequencies ($\le 200$\,Hz). Staircasing mitigation (cut-cell, immersed boundary) and non-rectangular subdomains are detailed in Section 4.3.1.

---

### 2. Analysis, Insights, and Commentary on Numerical Results
> **Reviewer's remark:**  
> *"Numerical results on the 6 benchmarks are presented. Unfortunately, there is very little comments or analysis of the results by the author, beyond presenting them. (...) For the two damped scenarios (the most novel contributions), there is in fact no comments at all by the author, which is a problematic omission. No conclusion is offered for this chapter."*

**Response:**  
We resolved this omission. In Chapter 2 (Sections 2.3.2 and 2.3.4) and Chapter 3 (Section 3.3), we added thorough physical and mathematical discussions:
1. **Viscoelastic Regularization of FDTD (Table 2.2):** In Experiment 4 (damped non-smooth pulse), FDTD achieves a convergence rate of $\approx 2.22$, which exceeds its undamped rate of $1.03$ (Experiment 2). We explained that physical viscosity $\nu > 0$ introduces modal decay scaling as $s_- \sim -\nu k_m^2$, acting as a natural parabolic low-pass filter that smooths sharp corners over time and restores second-order algebraic convergence.
2. **Low-Frequency Attenuation Expansion (Section 1.3.1):** We derived the exact low-frequency attenuation $\alpha(\omega) \sim \frac{\nu}{2c^3}\omega^2$ for pure viscoelastic damping ($\gamma = 0$), capturing thermoviscous dissipation. For combined models ($\gamma > 0$), we derived the diffusive root $\alpha(\omega) \sim \sqrt{\frac{\gamma\omega}{2c^2}}$ and the weak-damping asymptotic $\alpha(\omega) \sim \frac{\gamma}{2c} + \frac{\nu}{2c^3}\omega^2$.
3. **High-Frequency Modal Asymptotics (Section 2.3.1):** We corrected the modal root scaling for $k_m \to \infty$: rather than uniform $\nu k_m^2$ decay across both branches, the system becomes overdamped with two distinct branches: $s_+ \sim -c^2/\nu$ and $s_- \sim -\nu k_m^2$.

---

### 3. Chapter 4 Conclusion and Balance of Scientific Claims
> **Reviewer's remark:**  
> *"Chapter 4 offers a one-page conclusion for the 100-page thesis. (...) the conclusion only praises the work, sometimes overly so, and makes a number of exaggerated claims without mentioning any limitation. Three possible directions for future work are outlined in one sentence (...). Finally, the claim 'achieve unprecedented computational efficiency and spatial accuracy' has not been validated (...). Likewise, the claim 'significantly outperforming purely data-driven black-box inversions' is made with no evidence."*

**Response:**  
Chapter 4 was completely rewritten from the ground up:
1. **Summary of Main Contributions (Section 4.1):** Structured as continuous academic prose (free of bulleted lists and em dashes), synthesizing FDTD vs. PSTD analysis, OSWR domain decomposition, ARD for viscoelastic media, and differentiable physics inversion.
2. **Critical Limitations of the Present Study (Section 4.2):** Organized around four explicit pillars: geometric staircasing on curved boundaries, ARD time-splitting temporal accuracy, physical scope of synthetic DPNN inversions, and high-frequency computational scaling.
3. **Expanded Future Research Directions (Section 4.3):** Expanded from 3 lines to over 4 pages of rigorous technical prose across five subfields: geometric discretization (4.3.1), energy stability and exact exponential Kick splitting for ARD (4.3.2), analytical OSWR min-max theory (4.3.3), scaling DPNN to 3D and advanced inversion formulations (4.3.4), and frequency-dependent digital impedance boundary conditions and CPML (4.3.5).
4. **Retraction of Unsubstantiated Claims:**  
   - Retracted *"achieve unprecedented computational efficiency and spatial accuracy"*.
   - Retracted *"significantly outperforming purely data-driven black-box inversions"*.
   - Retracted the claim that ARD was extended to FDTD.

---

# Part 2: Point-by-Point Response to Major Comments

### Comment R2C1: Definition of Specialized Technical Terms
> **Reviewer's remark:**  
> *Sec. 1.3 and 1.4: Highly specialized terms are repeatedly used without defining or explaining them (e.g., "Runge-Kutta scheme"). Even though the presented material is of a textbook nature, someone outside of the field of numerical schemes could not easily decipher it due to its density.*

**Response:**  
In Sections 1.3 and 1.4, we added introductory explanations for all specialized numerical terminology upon first mention, including Runge-Kutta multistage integrators, Courant-Friedrichs-Lewy (CFL) stability bounds, Butcher tableaus, Krylov subspace methods, and spectral projection operators, citing standard treatises (Trefethen, 2000; Canuto et al., 2006; Hairer et al., 1993; Butcher, 2016).

---

### Comment R2C2: Numerical Results on Damped Benchmarks (Pages 66--67)
> **Reviewer's remark:**  
> *For the two damped scenarios (the most novel contributions), there is in fact no comments at all by the author, which is a problematic omission.*

**Response:**  
We added in-depth discussions for Experiments 4 and 6 in Section 2.3.4:
- *Experiment 4 (Viscoelastic non-smooth pulse):* FDTD convergence increases from order $1.03$ (undamped) to $2.22$, demonstrating that physical viscosity $\nu > 0$ rapidly attenuates high-frequency corner singularities ($s_- \sim -\nu k_m^2$), regularizing the wavefield and restoring nominal algebraic convergence.
- *Experiment 6 (Damped boundary reflection):* We clarified that PSTD achieves a $42.2\times$ reduction in error ($1.63$ orders of magnitude) over FDTD at $N=800$, with PSTD errors localized strictly to boundary stencils while FDTD suffers from grid dispersion throughout the domain.

---

### Comment R2C3: Motivation for Two Overlapping Subdomains in Sec. 3.1.4
> **Reviewer's remark:**  
> *Sec. 3.1.4: What is the motivation of using only 2 overlapping subdomains in 1D? Why not 3 or more, and why not non-overlapping?*

**Response:**  
In Section 3.1.4, we clarified that two overlapping subdomains $\Omega_1 = (0, b)$ and $\Omega_2 = (a, L)$ with $0 < a < b < L$ represent the canonical theoretical archetype in Schwarz Waveform Relaxation analysis (Gander et al., 2007; Dolean et al., 2015; Quarteroni and Valli, 1999). It isolates transmission operator convergence without combinatorial graph partition effects. We also proved in Section 3.1.2 that non-overlapping classical Schwarz with Dirichlet conditions stagnates, necessitating either overlap or optimized Robin conditions.

---

### Comment R2C4: Retraction of Claim that FDTD is Presented for ARD
> **Reviewer's remark:**  
> *P86: "We present the ARD method for both FDTD and PSTD, which share the same formulation": No, FDTD is not presented.*

**Response:**  
We retracted this claim. The revised text in Section 3.2 and Section 4.1 states that ARD is formulated and developed exclusively for PSTD, where analytical DCT-II modal basis functions allow exact spatial stepping.

---

### Comment R2C5: Justification of Parabolic Stability Limit
> **Reviewer's remark:**  
> *Sec. 3.2.3: "the second term introduces a strict parabolic stability limit ($\Delta t \le \Delta x^2/\nu$)." This is not sufficiently detailed or justified.*

**Response:**  
In Section 3.2.3, we added a formal Gershgorin circle derivation for the interface Kick operator. The discrete interface update matrix yields an upper eigenvalue bound bounded by the spatial stencil:
$$\Delta t \le \frac{4 \Delta x^2}{C_p \nu},$$
where $C_p$ is a dimensionless interface stencil constant. In Section 4.3.2, we also introduced the exact exponential Kick flow $\mathbf{v}(t+\tau) = e^{\tau \nu C_{\mathrm{res}}} \mathbf{v}(t) + \frac{c^2}{\nu} (e^{\tau \nu C_{\mathrm{res}}} - I) \mathbf{u}(t)$, which completely eliminates this parabolic restriction.

---

### Comment R2C6: Initialization of the Iterative ARD Algorithm
> **Reviewer's remark:**  
> *Sec. 3.2.4: How is the presented iterative algorithm initialised, i.e., what are the initial values of $U_m^n$ and $V_m^n$?*

**Response:**  
In Section 3.2.4, we explicitly stated that at initial time $t = 0$, modal coefficients are computed via forward DCT-II projections: $U_{m,0} = \mathrm{DCT}(u_0)$ and $V_{m,0} = \mathrm{DCT}(v_0)$. At internal communication steps, boundary interface ghost cells are initialized to zero and updated through the Kick-Drift-Kick exchange.

---

### Comment R2C7: Softening "Validating 3D ARD" and Adding Quantitative Evaluation
> **Reviewer's remark:**  
> *P96: "validating the accuracy of the multi-domain pseudo-spectral time-domain (PSTD) formulation in 3D." This is too strong a claim. A single configuration is tested, with two (very close) values of the viscoelastic parameter. Only qualitative snapshots are displayed.*

**Response:**  
We replaced "validating" with "demonstrating the operational feasibility of 3D ARD". We expanded the benchmark to a 5-solver quantitative comparison across four damping regimes with measured CPU execution times and spatial profiles.

---

### Comment R2C8: Internal Boundary Handling in 3D ARD Snapshots
> **Reviewer's remark:**  
> *Fig. 3.12: "The solver seamlessly handles the internal boundaries between the three rectangular subdomains." One cannot actually see that from the figure.*

**Response:**  
Figure 3.14 was replaced with high-resolution snapshots showing wavefield continuity across the virtual interface lines $\Gamma_1$ and $\Gamma_2$ marked with prominent dashed red lines, confirming the absence of unphysical reflections or phase distortions.

---

### Comment R2C9: 3D Damping Visualization and Physical Scaling
> **Reviewer's remark:**  
> *Fig. 3.16 is not very informative. It is impossible to see differences between the 3 images. Plotting the error or relative error would be better.*

**Response:**  
We agree that plotting raw wavefield heatmaps makes visual comparison difficult because all solvers capture macroscopic wavefronts similarly. Following your recommendation, we completely redesigned the visual evaluation in both the 3D ARD and 1D DPNN sections:
1. **Spatiotemporal Error Heatmaps in DPNN (Figures 3.16, 3.19, and 3.22):** Rather than displaying three virtually identical wavefield heatmaps, each spatiotemporal figure now presents a four-panel layout on a unified color scale:
   - Panel (a): Ground-truth wavefield $u^*(x,t)$ (seismic bipolar colormap).
   - Panel (b): FDTD-DPNN absolute error $|u_{\mathrm{FDTD}}(x,t) - u^*(x,t)|$.
   - Panel (c): PSTD-DPNN absolute error $|u_{\mathrm{PSTD}}(x,t) - u^*(x,t)|$.
   - Panel (d): Standard PINN baseline absolute error $|u_{\mathrm{PINN}}(x,t) - u^*(x,t)|$.
   This error visualization makes the numerical properties immediately transparent: PSTD-DPNN maintains near-zero error across the entire domain, FDTD-DPNN clearly reveals trailing wave-like error ripples caused by numerical grid dispersion ($c_{\mathrm{num}} < c$), and PINN exhibits large diffuse errors and phase drift throughout the interior.
2. **3D Damping Comparisons (Figure 3.15):** We updated the 3D ARD snapshots with standardized colorbars and added an explicit physical scaling analysis explaining why primary wavefront amplitudes decay modestly over $T = 0.24$\,s for nominal atmospheric damping ($\nu \approx 1.2\cdot 10^{-5}$\,m$^2$/s), while the reverberant coda undergoes progressive high-frequency attenuation.

---

### Comment R2C10: Explicit Description of Cosine Annealing in DPNN
> **Reviewer's remark:**  
> *The loss plots (Fig. 3.14, 3.17, 3.20...) feature abrupt changes, suggesting that some elaborate techniques are used to modify the learning rate along the way, but there is no explicit mention of this despite a mention of "cosine annealing" which is not described at all.*

**Response:**  
In Section 3.3.1, we provided the analytical equation for the Cosine Annealing learning rate schedule (Loshchilov and Hutter, 2017) and explained the dual-scheduling strategy:
- Network weights $\Theta$ are optimized with a conservative schedule ($\eta \in [10^{-5}, 10^{-3}]$) to preserve spatial representation stability.
- Logarithmic physical parameters $(\xi_\gamma, \xi_\nu)$ start with an initial learning rate $\eta = 5\cdot 10^{-2}$ to traverse non-convex loss valleys and escape poor local minima, and are annealed down to $10^{-3}$ to eliminate limit cycles and lock stably onto ground-truth values.

---

### Comment R2C11: Critical Scope and Limitations of DPNN Experiments
> **Reviewer's remark:**  
> *Sec 3.3.1: As in chapter 3, the comments on the experiments are extremely limited. This section also severely lacks quantitative analysis. There are 3 qualitative experiments, which seems too limited to rigorously validate the approach. Especially since quite strong limitations seem to appear in the third experiment.*

**Response:**  
In the revised thesis, we completely overhauled and upgraded Section 3.3 around a focused, rigorous narrative of essential experimental investigation blocks, adding extensive quantitative benchmarks, a high-capacity baseline comparison, disentangled inversions, and regularization ablations:
1. **Comprehensive Benchmark Comparison and PINN Baseline (Section 3.3.2, Tables 3.2 and 3.3):** We integrated a complete numerical comparison reporting four quantitative error metrics across all three initial condition benchmarks (Multimodal Standing Waves, Smooth Gaussian Pulses, Non-Smooth Triangular Pulses) comparing our proposed PSTD-DPNN against FDTD-DPNN and a high-capacity Physics-Informed Neural Network (PINN) reference baseline (Raissi et al., 2019; Haghighat et al., 2021):
   - Normalized Mean Squared Error on initial pressure: $\mathrm{NMSE}(p_0) = \|p_0^{\mathrm{pred}} - p_0^*\|_2^2 / \|p_0^*\|_2^2$.
   - Full spatiotemporal wavefield relative $L_2$ error: $\|u_{\mathrm{pred}} - u^*\|_2 / \|u^*\|_2$ (%).
   - Physical attenuation parameter estimation: telegrapher damping $\hat{\gamma}$ and viscoelastic damping $\hat{\nu}$ alongside relative errors.
   - Computational wall-clock training time per epoch (ms/epoch) evaluated on an Intel Xeon CPU workstation with an NVIDIA RTX GPU using PyTorch.
   - *Key findings*: On non-smooth triangular pulses, PSTD-DPNN attains $\mathrm{NMSE} = 4.49 \times 10^{-3}$ and a field relative error of $3.09\%$ with physical damping identified within $4.7\%$ for $\gamma$ and $3.6\%$ for $\nu$, whereas FDTD-DPNN yields $14.67\%$ field error and PINN yields $28.20\%$. While PINN fits smooth standing waves accurately ($\mathrm{NMSE} = 1.18 \times 10^{-3}$, field error $2.80\%$), its soft-collocation formulation struggles with sharp gradient discontinuities. In computational cost, PSTD-DPNN requires 12.1--13.8\,ms/epoch compared to 91.3--95.1\,ms/epoch for PINN (sampled with 2000 interior points), where automatic differentiation of mixed third-order derivatives $-\nu \partial_{txx} u$ contributes to the computational difference.
2. **Disentangled Inversion Experiments (Section 3.3.3, Table 3.4):** To isolate the mechanisms of source-parameter coupling, we disentangled the inverse problem into three isolated regimes on the Gaussian benchmark:
   - *Experiment A ($u_0$ reconstruction only, $\gamma, \nu$ known)*: Attains an initial condition NMSE of $3.23 \times 10^{-3}$ and a field error of $3.94\%$.
   - *Experiment B (Parameter identification only, $u_0$ known)*: Under logarithmic parameterization $\gamma = \gamma_{\mathrm{ref}} e^{\xi_\gamma}, \nu = \nu_{\mathrm{ref}} e^{\xi_\nu}$ initialized at physical $(0.50, 1.00 \times 10^{-3})$, identifies the viscosity $\nu$ within $3.0\%$ ($\hat{\nu} = 4.85 \times 10^{-3}\,\mathrm{m}^2/\mathrm{s}$) and damping $\gamma$ within $10.8\%$ ($\hat{\gamma} = 0.1108\,\mathrm{s}^{-1}$), with field error $0.37\%$, confirming that physical parameters are accurately recovered when decoupled from source variations.
   - *Experiment C (Joint inversion, $u_0, \gamma, \nu$ unknown)*: Directly aligns with the main Gaussian benchmark (Table 3.2), achieving $1.68\%$ field error (NMSE $3.00 \times 10^{-3}$), recovering $\nu$ within $2.0\%$ ($\hat{\nu} = 5.10 \times 10^{-3}\,\mathrm{m}^2/\mathrm{s}$) and $\gamma$ within $3.5\%$ ($\hat{\gamma} = 0.0965\,\mathrm{s}^{-1}$).
3. **Regularization Ablation Study and Elimination of $H^1$ (Section 3.3.4 & Table 3.5):** We conducted a systematic 5-configuration ablation study (Data only, $L_1$ only, $H^1$ only, Coupled $L_1+H^1$, Decoupled $L_1+H^1$):
   - *Data fitting without explicit spatial regularization*: Relying on the implicit inductive bias of the low-capacity SIREN parameterization enables clean convergence without hyperparameter distortion, achieving $1.68\%$ field error on Gaussian and $3.09\%$ on Triangle pulses.
   - *Why Sobolev $H^1$ fails and was rejected*: Penalizing the $H^1$ gradient semi-norm $\|\partial_x p_0\|_2^2 \propto \sum k_m^2 |U_{0,m}|^2$ suppresses high-wavenumber modes, eliminating the spectral information required to estimate physical dissipation ($\nu k_m^2$). When coupled to $\nu$ in the loss, the optimizer collapses $\hat{\nu} \to 10^{-5}$ to minimize the penalty. Decoupling $\nu$ does not prevent high-frequency loss. Therefore, $H^1$ regularization was rejected from the proposed formulation.

The findings are synthesized in the Concluding Assessment (Section 3.3.5).

---

### Comment R2C12: Retraction of "Remarkable Robustness" on Triangular Pulse
> **Reviewer's remark:**  
> *P104: "Figures 3.20 to 3.22 showcase the remarkable robustness of the carefully calibrated Sobolev H1 regularization". To be honest, not really. The relative reconstruction error on the initial condition seems very large. This claim should be retracted or further expanded.*

**Response:**  
We completely retracted the phrase "remarkable robustness". Following your critique, our new systematic ablation study (Section 3.3.4, Table 3.5) revealed that Sobolev $H^1$ regularization was actually detrimental rather than robust for this non-smooth profile:
- The proposed unregularized formulation achieves an initial condition NMSE of only $4.49 \times 10^{-3}$ and a full-field relative error of $3.09\%$, recovering $\gamma$ within $4.7\%$ and $\nu$ within $3.6\%$.
- In contrast, imposing Sobolev $H^1$ regularization degrades field error to $18.10\%$ and causes estimated physical viscosity $\nu$ to collapse to $7.66 \times 10^{-6}$ ($99.85\%$ relative error).
- We explain the physical mechanism in Section 3.3.4: penalizing the $H^1$ gradient semi-norm $\|\partial_x p_0\|_2^2 \propto \sum k_m^2 |U_{0,m}|^2$ introduces excessive smoothing on slope discontinuities and penalizes high spatial wavenumbers, systematically eliminating the spectral information required to estimate physical dissipation. Furthermore, when coupled to $\nu$ in the loss, the optimizer trades off physical viscosity to minimize the regularization penalty.

---

### Comment R2C13: Toning Down "Unprecedented Computational Efficiency"
> **Reviewer's remark:**  
> *P107: "achieve unprecedented computational efficiency and spatial accuracy." There is no solid proof of this in the thesis, this claim should be tuned down or retracted.*

**Response:**  
We completely retracted this claim in both Chapter 3 and Chapter 4. It was replaced with measured multi-threaded benchmarks: coarse ARD-PSTD is $1.63\times$ slower than coarse FDTD per step due to 3D DCT-II transforms, but executes $2.57\times$ faster than fine FDTD and $4.0\times$ faster than fine TD-FEM while completely eliminating numerical dispersion for resolved frequencies.

---

### Comment R2C14: Retraction of Claim Extending ARD to FDTD
> **Reviewer's remark:**  
> *P107: "Subsequently, we extended the Adaptive Rectangular Decomposition (ARD) method to both FDTD and PSTD frameworks.": No, this is incorrect, the FDTD framework is not treated. The author mentions it in a paragraph, but then explicitly writes that it won't be treated.*

**Response:**  
We corrected this statement in Chapter 4 (Summary of Main Contributions) and Section 3.2. The text now explicitly states that ARD was developed, formulated, and implemented exclusively for the PSTD framework.

---

### Comment R2C15: Toning Down Superlative Claims on Sobolev Regularization
> **Reviewer's remark:**  
> *P108: "(...) Sobolev H_1 smoothness penalty gracefully suppressed high-frequency artifacts without smearing out sharp discontinuities": This is a highly subjective and superlative statement. In fact, there is no real rigorous treatment of this, just a succinct qualitative experiment.*

**Response:**  
We completely retracted "gracefully suppressed (...) without smearing out sharp discontinuities". We replaced this subjective claim with the rigorous 5-configuration ablation study in Section 3.3.4 (Table 3.5) and the concluding assessment in Section 3.3.5:
- We explicitly demonstrate that Sobolev $H^1$ regularization smears out sharp corners and does not preserve slope discontinuities.
- We show that $H^1$ penalization introduces excessive smoothing ($7.42\%$ error on Gaussian, $18.10\%$ on triangle) and creates an algorithmic trap where the optimizer drives physical viscosity $\nu \to 0$ to minimize the gradient penalty.
- We conclude that explicit regularization is not universally beneficial: for smooth profiles, an $L^1$ penalty provides modest benefits ($2.42\%$ error), whereas for triangular profiles, the unregularized formulation achieves the lowest error ($3.09\%$) and the most accurate recovery of physical attenuation parameters. Consequently, $H^1$ regularization was rejected from the proposed DPNN formulation.

---

### Comment R2C16: Retraction of Unbacked Claim on Black-Box Inversions
> **Reviewer's remark:**  
> *P108: "(...), significantly outperforming purely data-driven black-box inversions": this claim has no supporting evidence and should be retracted. No comparison is provided for this method, nor any quantitative results.*

**Response:**  
We retracted the unbacked phrase "significantly outperforming purely data-driven black-box inversions". Furthermore, we replaced it with the concrete quantitative comparison in Table 3.2 comparing our Differentiable Physics solvers against a high-capacity Physics-Informed Neural Network (PINN) baseline (4 layers, 128 units, $\sim 50{,}000$ parameters) (Raissi et al., 2019; Haghighat et al., 2021). The experimental data rigorously demonstrates that:
- Structurally embedding the numerical wave propagator achieves lower reconstruction error on non-smooth initial conditions (PSTD-DPNN field error $3.09\%$ vs PINN $28.20\%$).
- PSTD-DPNN provides a $6.6\times$ to $7.7\times$ speedup per training epoch compared to soft-collocation PINNs, due to evaluating differential operators algebraically in modal space rather than through automatic differentiation of stiff third-order mixed derivatives.
- DPNN requires only $\sim 1{,}155$ parameters ($43\times$ fewer than the 50k-parameter PINN) by parameterizing only the initial spatial state $p_0(x)$ while leaving spatiotemporal propagation to the exact physical solver.

---

### Comment R2C17: Expansion of Future Research Directions
> **Reviewer's remark:**  
> *P108: 3 lines on future work is an order of magnitude less than what one would typically expect for a 3-year-long doctoral work.*

**Response:**  
Section 4.3 was expanded from 3 lines to over 4 pages of continuous technical prose, detailing 5 distinct research avenues:
1. *Section 4.3.1:* Geometric discretization and staircasing boundary error mitigation (Cartesian cut-cell, embedded boundary, curvilinear transformations).
2. *Section 4.3.2:* Energy stability and higher-order time integration for ARD (exact exponential Kick flow, midpoint Runge-Kutta, discrete energy estimates, high-order DtN Padé operators).
3. *Section 4.3.3:* Analytical min-max theory and adaptive transmission for OSWR.
4. *Section 4.3.4:* Scaling DPNN to 3D and advanced inversion formulations (3D scaling, optimal sensor placement, edge-preserving Total Variation regularization, and comparison with FWI and TRA).
5. *Section 4.3.5:* Realistic acoustic modelling, advanced boundary conditions (CPML/CFS-PML with auxiliary memory variables, digital impedance filters), and GPU acceleration.

---

# Part 3: Point-by-Point Response to Minor Comments and Typos

- **R2m1 (P4: Timeline of wave-based methods):**  
  *Clarification:* We clarified Section 1.1: geometrical acoustics ray-tracing methods emerged in the late 1960s and early 1970s (Krokstad et al., 1968), whereas time-dependent wave-based numerical schemes (such as Yee's FDTD) were introduced in computational electromagnetics in 1966 and were later adopted in architectural acoustics in the late 1990s (Botteldooren, 1995).
- **R2m2 (P12: Structural placement of DPNN paragraphs):**  
  *Clarification:* We restructured Section 1.4 into dedicated, clearly demarcated subsections: Section 1.4.1 details forward numerical solvers, while Section 1.4.2 formally introduces inverse acoustic problems and motivates Differentiable Physics Neural Networks.
- **R2m3 (P18: Equation reference for $k(i\omega)$):**  
  *Fix:* We assigned an equation label to the continuous dispersion relation in Section 2.1 and corrected the cross-reference to point directly to Equation (2.1).
- **R2m4 (P24: Clarification of the three models):**  
  *Clarification:* We clarified the text: the three physical models compared in Section 2.1 are (1) the lossless undamped wave equation ($\gamma = \nu = 0$), (2) the telegrapher model ($\gamma > 0, \nu = 0$), and (3) the combined telegrapher and viscoelastic model ($\gamma > 0, \nu > 0$).
- **R2m5 (P40: Role of Proposition 2.4):**  
  *Clarification:* We added a concise remark in Section 2.2 explaining that Proposition 2.4 establishes the exact modal transition matrix under combined viscoelastic and telegrapher dissipation, providing the necessary mathematical foundation for the modal PSTD stability proofs in Proposition 2.5.
- **R2m6 (P54: Clarification of "advanced in the modal domain"):**  
  *Clarification:* We replaced this sentence with a clear operational description: *"The spatial acoustic field is transformed into discrete modal coefficients via the DCT-II, advanced in time analytically via decoupled modal transition operators, and synthesized back into physical space via the inverse DCT-II."*
- **R2m7 (P71: Grammatical structure fix in Section 3.1):**  
  *Fix:* We corrected the sentence syntax by removing the redundant relative clause: *"In the classical Schwarz Waveform Relaxation method, Dirichlet transmission conditions are employed, whereas Optimized Schwarz Waveform Relaxation utilizes generalized Robin transmission conditions."*
- **R2m8 (P81: Domain overlap $a < b$):**  
  *Fix:* In Section 3.1.4, we explicitly stated that $0 < a < b < L$, defining the positive overlap length as $\delta = b - a > 0$.
- **R2m9 (P86: Citation typo):**  
  *Fix:* We fixed the citation command in Section 3.2 to use `\citep{strang_construction_1968}`, eliminating the duplicate author name.
- **R2m10 (P89: Citation typo):**  
  *Fix:* We corrected the inline citation to standard parenthetical format: *"Strang splitting \citep{strang_construction_1968}"*.
- **R2m11 (P90: Definition of $D$ and $K$ operators):**  
  *Fix:* In Section 3.2.2, we provided formal mathematical definitions: $D(\Delta t) = \mathcal{S}_{\text{local}}(\Delta t)$ denotes the exact modal time-evolution operator for the uncoupled subdomain interior, while $K(\tau) = I + \tau B$ denotes the discrete Kick operator integrating the singular interface transmission operator.

---

We believe that these extensive revisions, enhanced quantitative benchmarks, double comparison against the PINN baseline, and rigorous error analyses directly and comprehensively address every concern you raised. We thank you once again for your constructive guidance.

Sincerely,  
**Gerardo Cicalese**  
*Ph.D. Candidate in Information Technology, Politecnico di Milano*  
**Prof. Ilario Mazzieri**  
*Advisor, Associate Professor, Politecnico di Milano*  
**Prof. Gabriele Ciaramella**  
*Co-Advisor, Associate Professor, Politecnico di Milano*  
**Prof. Fabio Antonacci**  
*Co-Advisor, Associate Professor, Politecnico di Milano*
