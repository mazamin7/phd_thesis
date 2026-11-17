import os

tex_content = r"""\documentclass[11pt,a4paper]{article}

\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage[margin=2.5cm]{geometry}
\usepackage{amsmath,amssymb,amsfonts}
\usepackage{booktabs}
\usepackage{xcolor}
\usepackage{enumitem}
\usepackage{titlesec}
\usepackage{fancyhdr}
\usepackage[most]{tcolorbox}
\usepackage[round,authoryear]{natbib}
\usepackage[hidelinks,colorlinks=true,linkcolor=blue!70!black,citecolor=blue!70!black,urlcolor=blue!70!black]{hyperref}

% Page styling
\setlength{\headheight}{14pt}
\pagestyle{fancy}
\fancyhf{}
\lhead{\textsc{Response to Reviewer 2 --- Dr.~Antoine Deleforge}}
\rhead{\textsc{Gerardo Cicalese}}
\cfoot{\thepage}

\definecolor{commentbg}{RGB}{245, 247, 250}
\definecolor{commentframe}{RGB}{70, 95, 130}
\definecolor{responseframe}{RGB}{40, 110, 70}

\newtcolorbox{revcomment}[2][]{
    enhanced,
    breakable,
    colback=commentbg,
    colframe=commentframe,
    fonttitle=\bfseries,
    title={#2},
    boxrule=0.8pt,
    arc=2mm,
    top=2mm, bottom=2mm, left=3mm, right=3mm,
    before skip=3mm, after skip=2mm,
    #1
}

\newenvironment{response}{
    \par\vspace{1mm}\noindent\textbf{\textcolor{responseframe}{Response:}}\par\vspace{1mm}
}{
    \par\vspace{3mm}
}

\begin{document}

\title{\textbf{\Large Formal Response to Reviewer 2 (Dr.~Antoine Deleforge)}\\[0.3em]
\large Re: Ph.D. Dissertation Review\\\emph{``Domain Decomposition and Differentiable Physics for the Viscoelastic Wave Equation''}}
\author{\textbf{Candidate:} Gerardo Cicalese\\\textbf{Advisor:} Prof.~Fabio Antonacci\\\textbf{Institution:} Dipartimento di Elettronica, Informazione e Bioingegneria, Politecnico di Milano}
\date{November 2026}
\maketitle

\noindent Dear Dr.~Antoine Deleforge,\\[0.5em]
We express our deepest gratitude for your comprehensive, constructive, and rigorous review of the doctoral dissertation \emph{``Domain Decomposition and Differentiable Physics for the Viscoelastic Wave Equation''}. Your critical feedback has pinpointed key areas requiring physical nuance, formal mathematical correction, toned-down scientific claims, and substantial experimental expansion.

In response to your guidance, we conducted an exhaustive revision across the manuscript. In accordance with the revision guidelines, \textbf{all textual modifications, newly derived mathematical analyses, and expanded experimental discussions in the revised dissertation are highlighted in \textcolor{red}{red}} (\texttt{\textbackslash textcolor\{red\}\{...\}}).

Below, we provide a point-by-point response to your General Assessment, all 17 Major Comments, and all 11 Minor Comments and Typos.

\tableofcontents
\newpage

% ==============================================================================
\section{Part 1: Response to the General Assessment}
% ==============================================================================

\subsection{Scope, Applicative Impact, and 3D Experimental Validation}
\begin{revcomment}{General Assessment: Applicative Impact and 3D Validation}
``However, the fact that only one qualitative experiment in 3D is presented limits the thesis' potential applicative impact (...) For the audio and acoustic signal processing community, a lot remains to be explored in terms of improving and accelerating commonly used simulators.''
\end{revcomment}

\begin{response}
We agree that presenting a single qualitative 3D snapshot comparison limited the practical evaluation of the proposed solvers. In the revised dissertation, we significantly expanded Section~3.4 (3D ARD) and Section~4.2 (Limitations):
\begin{enumerate}
    \item \textbf{Comprehensive Multi-Solver 3D Benchmark:} We expanded the 3D room simulation ($L_x \times L_y \times L_z = 6.0\,\mathrm{m} \times 4.0\,\mathrm{m} \times 3.0\,\mathrm{m}$, $T = 0.24$\,s) into a systematic benchmark comparing five configurations: coarse FDTD ($h = 0.10$\,m), fine FDTD ($h = 0.05$\,m), coarse TD-FEM, fine TD-FEM, and our 3-subdomain ARD-PSTD solver.
    \item \textbf{Quantitative Execution Runtimes:} We added measured 8-core CPU runtimes: ARD-PSTD requires 102.4\,s per run ($1.63\times$ slower than coarse FDTD per step due to 3D DCT-II transforms, but $2.57\times$ faster than fine FDTD at 263.8\,s, and $4.0\times$ faster than fine TD-FEM at 409.6\,s).
    \item \textbf{Damping Parameter Scaling:} We added wavefield comparisons across four damping regimes: lossless ($\nu = 0$), low damping ($\nu \approx 2.4\cdot 10^{-6}$\,m$^2$/s), moderate damping ($\nu \approx 1.2\cdot 10^{-5}$\,m$^2$/s, corresponding to nominal atmospheric absorption at 20$^\circ$C and 50\% relative humidity), and high damping ($\nu \approx 5.9\cdot 10^{-5}$\,m$^2$/s). Over $T = 0.24$\,s with $k_{\max} \approx 15.7$\,m$^{-1}$, primary wavefronts undergo modest decay while the reverberant coda experiences progressive high-frequency attenuation.
    \item \textbf{Explicit Architectural Limitations:} In Section~4.2, we explicitly bounded the operational domain: the present 3D implementation is restricted to axis-aligned rectangular subdomains, rigid boundaries, and nominal audio frequencies ($\le 200$\,Hz). Staircasing mitigation (cut-cell, immersed boundary) and non-rectangular subdomains are detailed in Section~4.3.1.
\end{enumerate}
\end{response}

\subsection{Analysis, Insights, and Commentary on Numerical Results}
\begin{revcomment}{General Assessment: Commentary on Damped Results}
``Unfortunately, there is very little comments or analysis of the results by the author, beyond presenting them (...) For the two damped scenarios (the most novel contributions), there is in fact no comments at all by the author, which is a problematic omission.''
\end{revcomment}

\begin{response}
We resolved this omission. In Chapter~2 (Sections~2.3.2 and~2.3.4) and Chapter~3 (Section~3.5), we added thorough physical and mathematical discussions:
\begin{enumerate}
    \item \textbf{Viscoelastic Regularization of FDTD (Table~2.2):} In Experiment~4 (damped non-smooth pulse), FDTD achieves a convergence rate of $\approx 2.22$, which exceeds its undamped rate of $1.03$ (Experiment~2). We explained that physical viscosity $\nu > 0$ introduces modal decay scaling as $s_- \sim -\nu k_m^2$, acting as a natural parabolic low-pass filter that smooths sharp corners over time and restores second-order algebraic convergence.
    \item \textbf{Low-Frequency Attenuation Expansion (Section~1.3.1):} We derived the exact low-frequency attenuation $\alpha(\omega) \sim \frac{\nu}{2c^3}\omega^2$ for pure viscoelastic damping ($\gamma = 0$), capturing thermoviscous dissipation. For combined models ($\gamma > 0$), we derived the diffusive root $\alpha(\omega) \sim \sqrt{\frac{\gamma\omega}{2c^2}}$ and the weak-damping asymptotic $\alpha(\omega) \sim \frac{\gamma}{2c} + \frac{\nu}{2c^3}\omega^2$.
    \item \textbf{High-Frequency Modal Asymptotics (Section~2.3.1):} We corrected the modal root scaling for $k_m \to \infty$: rather than uniform $\nu k_m^2$ decay across both branches, the system becomes overdamped with two distinct branches: $s_+ \sim -c^2/\nu$ and $s_- \sim -\nu k_m^2$.
\end{enumerate}
\end{response}

\subsection{Chapter 4 Conclusion and Balance of Scientific Claims}
\begin{revcomment}{General Assessment: Concluding Chapter Balance}
``Finally, the conclusion is rather succinct and only 3 lines on future work is an order of magnitude less than what one would typically expect for a 3-year-long doctoral work. Some claims are too strong and should be qualified or retracted.''
\end{revcomment}

\begin{response}
Chapter~4 was completely rewritten from the ground up:
\begin{enumerate}
    \item \textbf{Summary of Main Contributions (Section~4.1):} Structured as continuous academic prose (free of bulleted lists and em dashes), synthesizing FDTD vs.\ PSTD analysis, OSWR domain decomposition, ARD for viscoelastic media, and differentiable physics inversion.
    \item \textbf{Critical Limitations of the Present Study (Section~4.2):} Organized around four explicit pillars: geometric staircasing on curved boundaries, ARD time-splitting temporal accuracy, physical scope of synthetic DPNN inversions, and high-frequency computational scaling.
    \item \textbf{Expanded Future Research Directions (Section~4.3):} Expanded from 3 lines to over 4 pages of rigorous technical prose across five subfields: geometric discretization (4.3.1), energy stability and exact exponential Kick splitting for ARD (4.3.2), analytical OSWR min-max theory (4.3.3), 3D scalable DPNN inversion (4.3.4), and frequency-dependent digital impedance boundary conditions and CPML (4.3.5).
\end{enumerate}
\end{response}

\newpage
% ==============================================================================
\section{Part 2: Point-by-Point Response to Major Comments}
% ==============================================================================

\begin{revcomment}{Comment R2C1: Definition of Specialized Technical Terms}
Sec. 1.3 and 1.4: Highly specialized terms are repeatedly used without defining or explaining them (e.g., ``Runge-Kutta scheme''). Even though the presented material is of a textbook nature, someone outside of the field of numerical schemes could not easily decipher it due to its density.
\end{revcomment}

\begin{response}
In Sections~1.3 and~1.4, we added introductory explanations for all specialized numerical terminology upon first mention, including Runge-Kutta multistage integrators, Courant-Friedrichs-Lewy (CFL) stability bounds, Butcher tableaus, Krylov subspace methods, and spectral projection operators, citing standard treatises \citep{trefethen_spectral_2000,canuto_spectral_2006,hairer_solving_1993,butcher_numerical_2016}.
\end{response}

\begin{revcomment}{Comment R2C2: Numerical Results on Damped Benchmarks (Pages 66--67)}
For the two damped scenarios (the most novel contributions), there is in fact no comments at all by the author, which is a problematic omission.
\end{revcomment}

\begin{response}
We added in-depth discussions for Experiments~4 and~6 in Section~2.3.4:
\begin{itemize}
    \item \emph{Experiment~4 (Viscoelastic non-smooth pulse):} FDTD convergence increases from order $1.03$ (undamped) to $2.22$, demonstrating that physical viscosity $\nu > 0$ rapidly attenuates high-frequency corner singularities ($s_- \sim -\nu k_m^2$), regularizing the wavefield and restoring nominal algebraic convergence.
    \item \emph{Experiment~6 (Damped boundary reflection):} We clarified that PSTD achieves a $42.2\times$ reduction in error ($1.63$ orders of magnitude) over FDTD at $N=800$, with PSTD errors localized strictly to boundary stencils while FDTD suffers from grid dispersion throughout the domain.
\end{itemize}
\end{response}

\begin{revcomment}{Comment R2C3: Motivation for Two Overlapping Subdomains in Sec. 3.1.4}
Sec. 3.1.4: What is the motivation of using only 2 overlapping subdomains in 1D? Why not 3 or more, and why not non-overlapping?
\end{revcomment}

\begin{response}
In Section~3.1.4, we clarified that two overlapping subdomains $\Omega_1 = (0, b)$ and $\Omega_2 = (a, L)$ with $0 < a < b < L$ represent the canonical theoretical archetype in Schwarz Waveform Relaxation analysis \citep{gander_optimized_2007,dolean_introduction_2015,quarteroni_domain_1999}. It isolates transmission operator convergence without combinatorial graph partition effects. We also proved in Section~3.1.2 that non-overlapping classical Schwarz with Dirichlet conditions stagnates, necessitating either overlap or optimized Robin conditions.
\end{response}

\begin{revcomment}{Comment R2C4: Retraction of Claim that FDTD is Presented for ARD}
P86: ``We present the ARD method for both FDTD and PSTD, which share the same formulation'': No, FDTD is not presented.
\end{revcomment}

\begin{response}
We retracted this claim. The revised text in Section~3.2 and Section~4.1 states that ARD is formulated and developed exclusively for PSTD, where analytical DCT-II modal basis functions allow exact spatial stepping.
\end{response}

\begin{revcomment}{Comment R2C5: Justification of Parabolic Stability Limit}
Sec. 3.2.3: ``the second term introduces a strict parabolic stability limit ($\Delta t \le \Delta x^2/\nu$).'' This is not sufficiently detailed or justified.
\end{revcomment}

\begin{response}
In Section~3.2.3, we added a formal Gershgorin circle derivation for the interface Kick operator. The discrete interface update matrix yields an upper eigenvalue bound bounded by the spatial stencil:
\begin{equation}
\Delta t \le \frac{4 \Delta x^2}{C_p \nu},
\end{equation}
where $C_p$ is a dimensionless interface stencil constant. In Section~4.3.2, we also introduced the exact exponential Kick flow $\mathbf{v}(t+\tau) = e^{\tau \nu C_{\mathrm{res}}} \mathbf{v}(t) + \frac{c^2}{\nu} (e^{\tau \nu C_{\mathrm{res}}} - I) \mathbf{u}(t)$, which completely eliminates this parabolic restriction.
\end{response}

\begin{revcomment}{Comment R2C6: Initialization of the Iterative ARD Algorithm}
Sec. 3.2.4: How is the presented iterative algorithm initialised, i.e., what are the initial values of $U_m^n$ and $V_m^n$?
\end{revcomment}

\begin{response}
In Section~3.2.4, we explicitly stated that at initial time $t = 0$, modal coefficients are computed via forward DCT-II projections: $U_{m,0} = \mathrm{DCT}(u_0)$ and $V_{m,0} = \mathrm{DCT}(v_0)$. At internal communication steps, boundary interface ghost cells are initialized to zero and updated through the Kick-Drift-Kick exchange.
\end{response}

\begin{revcomment}{Comment R2C7: Softening ``Validating 3D ARD'' and Adding Quantitative Evaluation}
P96: ``validating the accuracy of the multi-domain pseudo-spectral time-domain (PSTD) formulation in 3D.'' This is too strong a claim. A single configuration is tested, with two (very close) values of the viscoelastic parameter. Only qualitative snapshots are displayed.
\end{revcomment}

\begin{response}
We replaced ``validating'' with ``demonstrating the operational feasibility of 3D ARD''. We expanded the benchmark to a 5-solver quantitative comparison across four damping regimes with measured CPU execution times and spatial profiles.
\end{response}

\begin{revcomment}{Comment R2C8: Internal Boundary Handling in 3D ARD Snapshots}
Fig. 3.12: ``The solver seamlessly handles the internal boundaries between the three rectangular subdomains.'' One cannot actually see that from the figure.
\end{revcomment}

\begin{response}
Figure~3.14 was replaced with high-resolution snapshots showing wavefield continuity across the virtual interface lines $\Gamma_1$ and $\Gamma_2$ marked with prominent dashed red lines, confirming the absence of unphysical reflections or phase distortions.
\end{response}

\begin{revcomment}{Comment R2C9: 3D Damping Visualization and Physical Scaling}
Fig. 3.16 is not very informative. It is impossible to see differences between the 3 images. Plotting the error or relative error would be better.
\end{revcomment}

\begin{response}
We agree with the Reviewer that plotting raw wavefield heatmaps makes visual comparison difficult because all solvers capture the macroscopic wavefronts similarly. Following your recommendation, we have completely redesigned the visual evaluation in both the 3D ARD and 1D DPNN sections:
\begin{enumerate}
    \item \textbf{Spatiotemporal Error Heatmaps in DPNN (Figures~3.16, 3.19, and~3.22):} Rather than displaying three virtually identical wavefield heatmaps, each spatiotemporal figure now presents a four-panel layout on a unified color scale:
    \begin{itemize}
        \item Panel~(a): Ground-truth wavefield $u^*(x,t)$ (seismic bipolar colormap).
        \item Panel~(b): FDTD-DPNN absolute error $|u_{\mathrm{FDTD}}(x,t) - u^*(x,t)|$.
        \item Panel~(c): PSTD-DPNN absolute error $|u_{\mathrm{PSTD}}(x,t) - u^*(x,t)|$.
        \item Panel~(d): Standard PINN baseline absolute error $|u_{\mathrm{PINN}}(x,t) - u^*(x,t)|$.
    \end{itemize}
    This error visualization makes the numerical properties immediately transparent: PSTD-DPNN maintains near-zero error across the entire domain, FDTD-DPNN clearly reveals trailing wave-like error ripples caused by numerical grid dispersion ($c_{\mathrm{num}} < c$), and PINN exhibits large diffuse errors and phase drift throughout the interior.
    \item \textbf{3D Damping Comparisons (Figure~3.15):} We updated the 3D ARD snapshots with standardized colorbars and added an explicit physical scaling analysis explaining why primary wavefront amplitudes decay modestly over $T = 0.24$\,s for nominal atmospheric damping ($\nu \approx 1.2\cdot 10^{-5}$\,m$^2$/s), while the reverberant coda undergoes progressive high-frequency attenuation.
\end{enumerate}
\end{response}

\begin{revcomment}{Comment R2C10: Explicit Description of Cosine Annealing in DPNN}
The loss plots (Fig. 3.14, 3.17, 3.20...) feature abrupt changes, suggesting that some elaborate techniques are used to modify the learning rate along the way, but there is no explicit mention of this despite a mention of ``cosine annealing'' which is not described at all.
\end{revcomment}

\begin{response}
In Section~3.5, we provided the analytical equation for the Cosine Annealing learning rate schedule \citep{loshchilov2017sgdr} and explained the dual-scheduling strategy:
\begin{itemize}
    \item Network weights $\Theta$ are optimized with a conservative schedule ($\eta \in [10^{-5}, 10^{-3}]$) to preserve spatial representation stability.
    \item Logarithmic physical parameters $(\xi_\gamma, \xi_\nu)$ start with an initial learning rate $\eta = 5\cdot 10^{-2}$ to traverse non-convex loss valleys and escape poor local minima, and are annealed down to $10^{-3}$ to eliminate limit cycles and lock stably onto ground-truth values.
\end{itemize}
\end{response}

\begin{revcomment}{Comment R2C11: Critical Scope and Limitations of DPNN Experiments}
Sec 3.3.1: As in chapter 3, the comments on the experiments are extremely limited. This section also severely lacks quantitative analysis. There are 3 qualitative experiments, which seems too limited to rigorously validate the approach. Especially since quite strong limitations seem to appear in the third experiment.
\end{revcomment}

\begin{response}
We thank the Reviewer for this crucial observation. In the revised thesis, we completely overhauled and upgraded Section~3.5 around a focused, rigorous narrative of three essential experimental investigation blocks, adding extensive quantitative benchmarks, a high-capacity baseline comparison, disentangled inversions, and regularization ablations:
\begin{enumerate}
    \item \textbf{Comprehensive Benchmark Comparison and PINN Baseline (Section~3.5.1, Tables~3.2 and~3.3):} We integrated a complete numerical comparison reporting four quantitative error metrics across all three initial condition benchmarks (Multimodal Standing Waves, Smooth Gaussian Pulses, Non-Smooth Triangular Pulses) comparing our proposed PSTD-DPNN against FDTD-DPNN and a high-capacity Physics-Informed Neural Network (PINN) reference baseline \citep{raissi_pinn_2019,haghighat_sciann_2021}:
    \begin{itemize}
        \item Normalized Mean Squared Error on initial pressure: $\mathrm{NMSE}(p_0) = \|p_0^{\mathrm{pred}} - p_0^*\|_2^2 / \|p_0^*\|_2^2$.
        \item Full spatiotemporal wavefield relative $L_2$ error: $\|u_{\mathrm{pred}} - u^*\|_2 / \|u^*\|_2$ (\%).
        \item Physical attenuation parameter estimation: telegrapher damping $\hat{\gamma}$ and viscoelastic damping $\hat{\nu}$ alongside relative errors.
        \item Computational wall-clock training time per epoch (ms/epoch) evaluated on an Intel Xeon CPU workstation with an NVIDIA RTX GPU using PyTorch.
        \item \emph{Key findings}: On non-smooth triangular pulses, PSTD-DPNN attains $\mathrm{NMSE} = 4.49 \times 10^{-3}$ and a field relative error of $3.09\%$ with physical damping identified within $4.7\%$ for $\gamma$ and $3.6\%$ for $\nu$, whereas FDTD-DPNN yields $14.67\%$ field error and PINN yields $28.20\%$. While PINN fits smooth standing waves accurately ($\mathrm{NMSE} = 1.18 \times 10^{-3}$, field error $2.80\%$), its soft-collocation formulation struggles with sharp gradient discontinuities. In computational cost, PSTD-DPNN requires 12.1--13.8\,ms/epoch compared to 91.3--95.1\,ms/epoch for PINN (sampled with 2000 interior points), where automatic differentiation of mixed third-order derivatives $-\nu \partial_{txx} u$ contributes to the computational difference.
    \end{itemize}
    \item \textbf{Disentangled Inversion Experiments (Section~3.5.2, Table~3.4):} To isolate the mechanisms of source-parameter coupling, we disentangled the inverse problem into three isolated regimes on the Gaussian benchmark:
    \begin{itemize}
        \item \emph{Experiment A ($u_0$ reconstruction only, $\gamma, \nu$ known)}: Attains an initial condition NMSE of $3.23 \times 10^{-3}$ and a field error of $3.94\%$.
        \item \emph{Experiment B (Parameter identification only, $u_0$ known)}: Under logarithmic parameterization $\gamma = \gamma_{\mathrm{ref}} e^{\xi_\gamma}, \nu = \nu_{\mathrm{ref}} e^{\xi_\nu}$ initialized at physical $(0.50, 1.00 \times 10^{-3})$, identifies the viscosity $\nu$ within $3.0\%$ ($\hat{\nu} = 4.85 \times 10^{-3}\,\mathrm{m}^2/\mathrm{s}$) and damping $\gamma$ within $10.8\%$ ($\hat{\gamma} = 0.1108\,\mathrm{s}^{-1}$), with field error $0.37\%$, confirming that physical parameters are accurately recovered when decoupled from source variations.
        \item \emph{Experiment C (Joint inversion, $u_0, \gamma, \nu$ unknown)}: Directly aligns with the main Gaussian benchmark (Table~3.2), achieving $1.68\%$ field error (NMSE $3.00 \times 10^{-3}$), recovering $\nu$ within $2.0\%$ ($\hat{\nu} = 5.10 \times 10^{-3}\,\mathrm{m}^2/\mathrm{s}$) and $\gamma$ within $3.5\%$ ($\hat{\gamma} = 0.0965\,\mathrm{s}^{-1}$).
    \end{itemize}
    \item \textbf{Regularization Ablation Study and Elimination of $H^1$ (Section~3.5.3 \& Table~3.5):} We conducted a systematic 5-configuration ablation study (Data only, $L_1$ only, $H^1$ only, Coupled $L_1+H^1$, Decoupled $L_1+H^1$):
    \begin{itemize}
        \item \emph{Data fitting without explicit spatial regularization}: Relying on the implicit inductive bias of the low-capacity SIREN parameterization enables clean convergence without hyperparameter distortion, achieving $1.68\%$ field error on Gaussian and $3.09\%$ on Triangle pulses.
        \item \emph{Why Sobolev $H^1$ fails and was rejected}: Penalizing the $H^1$ gradient semi-norm $\|\partial_x p_0\|_2^2 \propto \sum k_m^2 |U_{0,m}|^2$ suppresses high-wavenumber modes, eliminating the spectral information required to estimate physical dissipation ($\nu k_m^2$). When coupled to $\nu$ in the loss, the optimizer collapses $\hat{\nu} \to 10^{-5}$ to minimize the penalty. Decoupling $\nu$ does not prevent high-frequency loss. Therefore, $H^1$ regularization was rejected from the proposed formulation.
    \end{itemize}
\end{enumerate}
The findings are synthesized in the Concluding Assessment (Section~3.5.4).
\end{response}

\begin{revcomment}{Comment R2C12: Retraction of ``Remarkable Robustness'' on Triangular Pulse}
P104: ``Figures 3.20 to 3.22 showcase the remarkable robustness of the carefully calibrated Sobolev H1 regularization''. To be honest, not really. The relative reconstruction error on the initial condition seems very large. This claim should be retracted or further expanded.
\end{revcomment}

\begin{response}
We completely retracted the phrase ``remarkable robustness''. Following your critique, our new systematic ablation study (Section~3.5.3, Table~3.5) revealed that Sobolev $H^1$ regularization was actually detrimental rather than robust for this non-smooth profile:
\begin{itemize}
    \item The proposed unregularized formulation achieves an initial condition NMSE of only $4.49 \times 10^{-3}$ and a full-field relative error of $3.09\%$, recovering $\gamma$ within $4.7\%$ and $\nu$ within $3.6\%$.
    \item In contrast, imposing Sobolev $H^1$ regularization degrades field error to $18.10\%$ and causes estimated physical viscosity $\nu$ to collapse to $7.66 \times 10^{-6}$ ($99.85\%$ relative error).
    \item We explain the physical mechanism in Section~3.5.3: penalizing the $H^1$ gradient semi-norm $\|\partial_x p_0\|_2^2 \propto \sum k_m^2 |U_{0,m}|^2$ introduces excessive smoothing on slope discontinuities and penalizes high spatial wavenumbers, systematically eliminating the spectral information required to estimate physical dissipation. Furthermore, when coupled to $\nu$ in the loss, the optimizer trades off physical viscosity to minimize the regularization penalty.
\end{itemize}
\end{response}

\begin{revcomment}{Comment R2C13: Toning Down ``Unprecedented Computational Efficiency''}
P107: ``achieve unprecedented computational efficiency and spatial accuracy.'' There is no solid proof of this in the thesis, this claim should be tuned down or retracted.
\end{revcomment}

\begin{response}
We completely retracted this claim in both Chapter~3 and Chapter~4. It was replaced with measured multi-threaded benchmarks: coarse ARD-PSTD is $1.63\times$ slower than coarse FDTD per step due to 3D DCT-II transforms, but executes $2.57\times$ faster than fine FDTD and $4.0\times$ faster than fine TD-FEM while completely eliminating numerical dispersion for resolved frequencies.
\end{response}

\begin{revcomment}{Comment R2C14: Retraction of Claim Extending ARD to FDTD}
P107: ``Subsequently, we extended the Adaptive Rectangular Decomposition (ARD) method to both FDTD and PSTD frameworks.'': No, this is incorrect, the FDTD framework is not treated. The author mentions it in a paragraph, but then explicitly writes that it won't be treated.
\end{revcomment}

\begin{response}
We corrected this statement in Chapter~4 (Summary of Main Contributions) and Section~3.4. The text now explicitly states that ARD was developed, formulated, and implemented exclusively for the PSTD framework.
\end{response}

\begin{revcomment}{Comment R2C15: Toning Down Superlative Claims on Sobolev Regularization}
P108: ``(...) Sobolev $H_1$ smoothness penalty gracefully suppressed high-frequency artifacts without smearing out sharp discontinuities'': This is a highly subjective and superlative statement. In fact, there is no real rigorous treatment of this, just a succinct qualitative experiment.
\end{revcomment}

\begin{response}
We completely retracted ``gracefully suppressed (...) without smearing out sharp discontinuities''. We replaced this subjective claim with the rigorous 5-configuration ablation study in Section~3.5.3 (Table~3.5) and the concluding assessment in Section~3.5.4:
\begin{itemize}
    \item We explicitly demonstrate that Sobolev $H^1$ regularization smears out sharp corners and does not preserve slope discontinuities.
    \item We show that $H^1$ penalization introduces excessive smoothing ($7.42\%$ error on Gaussian, $18.10\%$ on triangle) and creates an algorithmic trap where the optimizer drives physical viscosity $\nu \to 0$ to minimize the gradient penalty.
    \item We conclude that explicit regularization is not universally beneficial: for smooth profiles, an $L^1$ penalty provides modest benefits ($2.42\%$ error), whereas for triangular profiles, the unregularized formulation achieves the lowest error ($3.09\%$) and the most accurate recovery of physical attenuation parameters. Consequently, $H^1$ regularization was rejected from the proposed DPNN formulation.
\end{itemize}
\end{response}

\begin{revcomment}{Comment R2C16: Retraction of Unbacked Claim on Black-Box Inversions}
P108: ``(...), significantly outperforming purely data-driven black-box inversions'': this claim has no supporting evidence and should be retracted. No comparison is provided for this method, nor any quantitative results.
\end{revcomment}

\begin{response}
We retracted the unbacked phrase ``significantly outperforming purely data-driven black-box inversions''. Furthermore, we replaced it with the concrete quantitative comparison in Table~3.2 comparing our Differentiable Physics solvers against a high-capacity Physics-Informed Neural Network (PINN) baseline (4 layers, 128 units, $\sim 50{,}000$ parameters) \citep{raissi_pinn_2019,haghighat_sciann_2021}. The experimental data rigorously demonstrates that:
\begin{itemize}
    \item Structurally embedding the numerical wave propagator achieves lower reconstruction error on non-smooth initial conditions (PSTD-DPNN field error $3.09\%$ vs PINN $28.20\%$).
    \item PSTD-DPNN provides a $6.6\times$ to $7.7\times$ speedup per training epoch compared to soft-collocation PINNs, due to evaluating differential operators algebraically in modal space rather than through automatic differentiation of stiff third-order mixed derivatives.
    \item DPNN requires only $\sim 1{,}155$ parameters ($43\times$ fewer than the 50k-parameter PINN) by parameterizing only the initial spatial state $p_0(x)$ while leaving spatiotemporal propagation to the exact physical solver.
\end{itemize}
\end{response}

\begin{revcomment}{Comment R2C17: Expansion of Future Research Directions}
P108: 3 lines on future work is an order of magnitude less than what one would typically expect for a 3-year-long doctoral work.
\end{revcomment}

\begin{response}
Section~4.3 was expanded from 3 lines to over 4 pages of continuous technical prose, detailing 5 distinct research avenues:
\begin{enumerate}
    \item \emph{Section~4.3.1:} Geometric discretization and staircasing boundary error mitigation (Cartesian cut-cell, embedded boundary, curvilinear transformations).
    \item \emph{Section~4.3.2:} Energy stability and higher-order time integration for ARD (exact exponential Kick flow, midpoint Runge-Kutta, discrete energy estimates, high-order DtN Padé operators).
    \item \emph{Section~4.3.3:} Analytical min-max theory and adaptive transmission for OSWR.
    \item \emph{Section~4.3.4:} Scalability, identifiability, and architectural ablation for differentiable learning (3D scaling, sensor sparsity, ablation of SIREN vs Fourier features).
    \item \emph{Section~4.3.5:} Realistic acoustic modelling, advanced boundary conditions (CPML/CFS-PML with auxiliary memory variables, digital impedance filters), and GPU acceleration.
\end{enumerate}
\end{response}

\newpage
% ==============================================================================
\section{Part 3: Point-by-Point Response to Minor Comments and Typos}
% ==============================================================================

\begin{revcomment}{Comment R2m1: Timeline of Wave-Based Methods (P4)}
The author writes that ``The first time-dependent wave-based methods appeared around the same period of time [(1998)]'', but after that, papers from 1968, 1971 and 1974 are cited. Did time-dependent wave-based methods appear in the late 90s or early 70s? Please clarify.
\end{revcomment}

\begin{response}
We clarified Section~1.1: geometrical acoustics ray-tracing methods emerged in the late 1960s and early 1970s \citep{krokstad_calculating_1968}, whereas time-dependent wave-based numerical schemes (such as Yee's FDTD) were introduced in computational electromagnetics in 1966 and were later adopted in architectural acoustics in the late 1990s \citep{botteldooren_finitedifference_1995}.
\end{response}

\begin{revcomment}{Comment R2m2: Structural Placement of DPNN Paragraphs (P12)}
The three paragraphs on DPNNs talk about inverse problem while the 5 previous pages in the section are on numerical schemes for forward problems, so they seem a bit out of place. This should maybe be a different section.
\end{revcomment}

\begin{response}
We restructured Section~1.4 into dedicated, clearly demarcated subsections: Section~1.4.1 details forward numerical solvers, while Section~1.4.2 formally introduces inverse acoustic problems and motivates Differentiable Physics Neural Networks.
\end{response}

\begin{revcomment}{Comment R2m3: Equation Reference for $k(i\omega)$ (P18)}
``$k(i\omega)$ satisfies the dispersion relation (2.2).'': I think this is the wrong equation reference. This should probably point to the equation in section 2.1 (which doesn't have a number).
\end{revcomment}

\begin{response}
We assigned an equation label to the continuous dispersion relation in Section~2.1 and corrected the cross-reference to point directly to Equation~(2.1).
\end{response}

\begin{revcomment}{Comment R2m4: Clarification of the Three Models (P24)}
``An important qualitative difference between the three models considered...'' I am confused what are the three models referred to here. Section 2.1 only mentions two: multi-relaxation and telegrapher + viscoelastic damping. Is the undamped model considered as a third model?
\end{revcomment}

\begin{response}
We clarified the text: the three physical models compared in Section~2.1 are (1) the lossless undamped wave equation ($\gamma = \nu = 0$), (2) the telegrapher model ($\gamma > 0, \nu = 0$), and (3) the combined telegrapher and viscoelastic model ($\gamma > 0, \nu > 0$).
\end{response}

\begin{revcomment}{Comment R2m5: Role of Proposition 2.4 (P40)}
The usefulness of stating Proposition 2.4 in full detail is not clear, as this seems completely standard.
\end{revcomment}

\begin{response}
We added a concise remark in Section~2.2 explaining that Proposition~2.4 establishes the exact modal transition matrix under combined viscoelastic and telegrapher dissipation, providing the necessary mathematical foundation for the modal PSTD stability proofs in Proposition~2.5.
\end{response}

\begin{revcomment}{Comment R2m6: Clarification of ``Advanced in the Modal Domain'' (P54)}
``The numerical solution is advanced in the modal domain.'' I don't understand this sentence.
\end{revcomment}

\begin{response}
We replaced this sentence with a clear operational description: \emph{``The spatial acoustic field is transformed into discrete modal coefficients via the DCT-II, advanced in time analytically via decoupled modal transition operators, and synthesized back into physical space via the inverse DCT-II.''}
\end{response}

\begin{revcomment}{Comment R2m7: Grammatical Structure in Section 3.1 (P71)}
``For the classical Schwarz Waveform Relaxation method, for which Dirichlet transmission conditions are employed, whereas for Optimized Schwarz Waveform Relaxation...'' There appears to be a structure problem with this sentence (remove ``for which''?).
\end{revcomment}

\begin{response}
We corrected the sentence syntax by removing the redundant relative clause: \emph{``In the classical Schwarz Waveform Relaxation method, Dirichlet transmission conditions are employed, whereas Optimized Schwarz Waveform Relaxation utilizes generalized Robin transmission conditions.''}
\end{response}

\begin{revcomment}{Comment R2m8: Overlap Definition $a < b$ (P81)}
In the definition of $\Omega_1$ and $\Omega_2$: I assume that $a < b$ since the domains are overlapping? Please specify.
\end{revcomment}

\begin{response}
In Section~3.1.4, we explicitly stated that $0 < a < b < L$, defining the positive overlap length as $\delta = b - a > 0$.
\end{response}

\begin{revcomment}{Comment R2m9: Duplicate Citation Typo (P86)}
``(Strang Strang (1968))''
\end{revcomment}

\begin{response}
We fixed the citation command in Section~3.2 to use \texttt{\textbackslash citep\{strang\_construction\_1968\}}, eliminating the duplicate author name.
\end{response}

\begin{revcomment}{Comment R2m10: Citation Typo (P89)}
``Strang splitting Strang (1968)'': typo (one should use \texttt{\textbackslash citep} here).
\end{revcomment}

\begin{response}
We corrected the inline citation to standard parenthetical format: \emph{``Strang splitting \citep{strang_construction_1968}''}.
\end{response}

\begin{revcomment}{Comment R2m11: Definition of Drift and Kick Operators (P90)}
``$D(\Delta t)$ denotes the exact modal propagator and $K(\Delta t/2)$ represents...'' $D$ and $K$ are not adequately defined here.
\end{revcomment}

\begin{response}
In Section~3.2.2, we provided formal mathematical definitions: $D(\Delta t) = \mathcal{S}_{\text{local}}(\Delta t)$ denotes the exact modal time-evolution operator for the uncoupled subdomain interior, while $K(\tau) = I + \tau B$ denotes the discrete Kick operator integrating the singular interface transmission operator.
\end{response}

\vspace{0.8cm}
\noindent\rule{\textwidth}{0.4pt}\\[0.5em]
We believe that these extensive revisions, enhanced quantitative benchmarks, double comparison against the PINN baseline, and rigorous error analyses directly and comprehensively address every concern you raised. We thank you once again for your constructive guidance.

\vspace{0.5cm}
\noindent Sincerely,\\[0.3em]
\textbf{Gerardo Cicalese}\\
\emph{Ph.D. Candidate in Information Technology}\\
\emph{Politecnico di Milano}

\newpage
\bibliographystyle{unsrtnat-nourl}
\bibliography{Thesis_bibliography}

\end{document}
"""

with open('review/Response_to_Reviewer_2_Deleforge.tex', 'w', encoding='utf-8') as f:
    f.write(tex_content)

print("Generated clean review/Response_to_Reviewer_2_Deleforge.tex")
