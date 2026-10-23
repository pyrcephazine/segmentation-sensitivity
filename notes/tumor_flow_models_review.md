# Flow models for tumor-directed delivery

## 1. Models and strong forms

Start with the finite-radius lateral-average model of [Laurino–Zunino (2019)](#references), then change axial resistance, exchange location, or interface resistance separately. Compare models using **tumor delivery**, with healthy outlets retained because they compete for flow.

**Geometry and pressures.** Following Amare et al., let \(\Omega\) be the bulk tissue domain, \(x\in\Omega\) a spatial position, \(\Lambda\) the vessel centerline graph, and \(s\) arc length along a vessel. Write \(p_t(x)\) for bulk pressure and \(P(s)\) for vessel pressure. At position \(s\), let \(R(s)\) be vessel radius, \(\Theta(s)\) its cross section, \(|\Theta|=\pi R^2\) its area, and \(|\partial\Theta|=2\pi R\) its circumference.

**Flow coefficients.** Let \(\mu\) be dynamic viscosity, \(k_t\) bulk permeability, \(k_v(s)\) axial vessel permeability, \(\boldsymbol v_t\) bulk Darcy velocity, and \(q(s)\) axial volumetric flow. Define the axial coefficient \(a_v(s)=|\Theta(s)|k_v(s)/\mu\). A prime denotes differentiation with respect to \(s\). Then

$$
\boldsymbol v_t=-\frac{k_t}{\mu}\nabla p_t,\qquad
q=-a_vP',\qquad
k_v=\frac{R^2}{8},\qquad
 a_v=\frac{\pi R^4}{8\mu}
\quad\text{for a circular lumen}.
$$

Let \(P_{\rm CVP}\) be the downstream venous reference pressure and \(\alpha_v(x)\geq0\) a volumetric drainage coefficient, with units \(\mathrm{Pa}^{-1}\mathrm{s}^{-1}\). The term \(\alpha_v(p_t-P_{\rm CVP})\) is an optional physiological extension, not part of LZ’s original model. On the outer boundary \(\partial\Omega\), let \(\boldsymbol n\) denote the outward unit normal. Specify bulk pressure, no flux, or drainage there. All models use an arterial inlet pressure **or** inflow, continuous junction pressure, and conservation of junction flow. Artificial gap ends are closed, not new physiological outlets.

**Terminal exchange.** Let \(A_k\), \(k=1,\ldots,N_{\rm out}\), be the anatomical outlets; \(N_{\rm out}\) is their number. Let \(q(A_k)\) be the total flow leaving outlet \(A_k\), \(\epsilon>0\) the terminal source radius, and \(\eta_k^\epsilon(x)\) a nonnegative source kernel supported within distance \(\epsilon\) of \(A_k\), normalized on \(\Omega\). Define its pressure average \(\langle p_t\rangle_k\), the volumetric arterial source \(Q_{a,k}^\epsilon\), and the total terminal conductance \(G_{a,k}\), with units \(\mathrm{m}^3\mathrm{Pa}^{-1}\mathrm{s}^{-1}\), by

$$
\int_\Omega\eta_k^\epsilon\,dx=1,\qquad
\langle p_t\rangle_k=\int_\Omega\eta_k^\epsilon p_t\,dx,\qquad
Q_{a,k}^\epsilon(x)=q(A_k)\eta_k^\epsilon(x),
$$

$$
q(A_k)=G_{a,k}\bigl(P(A_k)-\langle p_t\rangle_k\bigr).
$$

Here \(G_{a,k}\) absorbs the area and viscosity factors accompanying Amare’s terminal coefficient \(\gamma_a\); it is a **total conductance**. In the coupled equations below, \(p_t\) describes an **unresolved blood bed** receiving terminal flow. If it instead describes **interstitium**, omit terminal blood sources from that compartment and use separate vascular downstream pressures \(P_{d,k}\) in \(q(A_k)=G_{a,k}(P(A_k)-P_{d,k})\).

**Tumor output.** Let \(q_\ell(s)\) be vessel-to-bulk exchange per unit vessel length, positive outward. Let \(\omega_k\in[0,1]\) be the tumor fraction supplied by outlet \(A_k\), and \(\omega_\ell(s)\in[0,1]\) the tumor fraction of lateral supply. The terminal tumor-flow output \(J_T\) is \(\sum_k\omega_kq(A_k)\). If lateral exchange also supplies the blood bed, use

$$
J_T=\sum_{k=1}^{N_{\rm out}}\omega_kq(A_k)
+\int_\Lambda\omega_\ell(s)q_\ell(s)\,ds.
$$

Keep these weights fixed across synthetic defects. Endpoint membership in a tumor mask only approximates supply territory. Do not count the same flow at an upstream feeder and its downstream outlets.

### N. Resistive graph baseline

Let \(i,j\) index connected graph vertices, \(P_i\) be nodal vessel pressure, \(q_{ij}\) flow from \(i\) to \(j\), \(g_{ij}\) their conductance, and \(f_i\) prescribed nodal inflow. Let \(\mathcal O\) be the set of outlet vertices; \(G_{a,i}\) and \(P_{d,i}\) are their downstream conductance and pressure. The indicator \(\mathbf1_{i\in\mathcal O}\) is one at an outlet and zero otherwise. For an edge \(e\), its conductance \(g_e\) satisfies

$$
q_{ij}=g_{ij}(P_i-P_j),\qquad
\sum_{j\sim i}g_{ij}(P_i-P_j)
+\mathbf1_{i\in\mathcal O}G_{a,i}(P_i-P_{d,i})=f_i,
\qquad
g_e=\left(\int_e\frac{ds}{a_v(s)}\right)^{-1}.
$$

This isolates geometry and connectivity from bulk coupling. Poiseuille gives \(g\propto R^4\); axial Darcy flow with **fixed** \(k_v\) gives \(a_v\propto R^2\). Setting \(k_v=R^2/8\) recovers Poiseuille, so the difference is the permeability assumption. [Derivation and discretization](#resistive-graph-derivation). [Numerical results](#radius-errors-and-repair).

### L. Finite-radius lateral average

Let \(\gamma(s)\) be wall hydraulic conductivity, with units \(\mathrm m\,\mathrm{Pa}^{-1}\mathrm{s}^{-1}\), and \(\overline{p_t}(s)\) the circumference-averaged bulk pressure:

$$
\overline{p_t}(s)=\frac1{|\partial\Theta(s)|}
\int_{\partial\Theta(s)}p_t\,d\ell,\qquad
q_\ell=\gamma|\partial\Theta|(P-\overline{p_t}).
$$

For any smooth scalar weighting function \(v_t\), let \(\overline{v_t}\) be its circumference average. Define the spreading operator \(S_R\) through its action under integration:
\(\langle S_Rq_\ell,v_t\rangle=\int_\Lambda q_\ell\overline{v_t}\,ds\). It deposits exchange on the finite lateral surface. The strong form is

$$
-\nabla\cdot\left(\frac{k_t}{\mu}\nabla p_t\right)
+\alpha_v(p_t-P_{\rm CVP})
=S_Rq_\ell+\sum_kQ_{a,k}^\epsilon,
\qquad -(a_vP')'+q_\ell=0.
$$

This retains LZ’s averaging while adding blood-flow coefficients, arterial boundaries, and the tumor output. Exchange means filtration for interstitial pressure, or effective unresolved vascular exchange for blood-bed pressure; these require different calibration. Finite-radius averaging avoids evaluating arbitrary bulk pressure on a line, but does not reconnect a broken artery. [Derivation and discretization](#lateral-average-derivation). [Numerical results](#gaps-and-tumor-delivery).

### T. Terminal-only exchange

$$
-\nabla\cdot\left(\frac{k_t}{\mu}\nabla p_t\right)
+\alpha_v(p_t-P_{\rm CVP})=\sum_kQ_{a,k}^\epsilon,
\qquad -(a_vP')'=0,
\qquad q(A_k)=G_{a,k}\bigl(P(A_k)-\langle p_t\rangle_k\bigr).
$$

This is L with \(\gamma=0\): resolved arteries feed unresolved microvasculature at anatomical outlets. A normalized sphere kernel gives a sphere-of-influence model; an anatomical territory can replace the sphere. Following [Amare et al. (2026)](#references), this tests **where exchange occurs**. Vary \(\epsilon\) separately from resistance \(1/G_{a,k}\). [Derivation and discretization](#terminal-exchange-derivation). [Numerical results](#gaps-and-tumor-delivery).

### C. Perfect-contact constraint

Let \(\lambda(s)\) be the unknown exchange flow per unit length needed to enforce equal averaged pressures:

$$
-\nabla\cdot\left(\frac{k_t}{\mu}\nabla p_t\right)
+\alpha_v(p_t-P_{\rm CVP})
=S_R\lambda+\sum_kQ_{a,k}^\epsilon,
\qquad -(a_vP')'+\lambda=0,
\qquad P-\overline{p_t}=0.
$$

This is the zero-interface-resistance limit, related to [Kuchta et al. (2021)](#references). It tests whether finite exchange resistance matters and avoids an arbitrarily large penalty coefficient. Perfect contact is an assumption, not automatic physiological improvement. [Derivation and discretization](#perfect-contact-derivation). [Numerical results](#gaps-and-tumor-delivery).

### R. Regularized distributed sources and the singular alternative

Let \(\epsilon_\ell>0\) be the lateral source width, \(\eta_{\epsilon_\ell}(x,s)\) a nonnegative normalized volume kernel, and \(\langle p_t\rangle_{\epsilon_\ell}\) its pressure average:

$$
\langle p_t\rangle_{\epsilon_\ell}(s)
=\int_\Omega\eta_{\epsilon_\ell}(x,s)p_t(x)\,dx,\qquad
\int_\Omega\eta_{\epsilon_\ell}(x,s)\,dx=1.
$$

Define \(S_{\epsilon_\ell}\) by
\(\langle S_{\epsilon_\ell}q_\ell,v_t\rangle=\int_\Lambda q_\ell\langle v_t\rangle_{\epsilon_\ell}\,ds\), and let \(\gamma_{\epsilon_\ell}\) be the effective wall conductivity at that averaging scale. Then

$$
\begin{aligned}
-\nabla\cdot\left(\frac{k_t}{\mu}\nabla p_t\right)
+\alpha_v(p_t-P_{\rm CVP})
&=S_{\epsilon_\ell}q_\ell+\sum_kQ_{a,k}^\epsilon,\\
-(a_vP')'+q_\ell&=0,\qquad
q_\ell=\gamma_{\epsilon_\ell}|\partial\Theta|
\bigl(P-\langle p_t\rangle_{\epsilon_\ell}\bigr).
\end{aligned}
$$

Volume exchange reduces sensitivity to mesh alignment but introduces a support scale and a pressure correction/calibration question; see [Koch et al. (2020)](#references). For the distinct singular alternative of [D’Angelo–Quarteroni (2008)](#references), let \(\delta_\Lambda\) be a source concentrated on the centerline and \(\Pi p_t\) a specified finite-radius pressure average. Replace the bulk lateral source by \(q_\ell\delta_\Lambda\), retaining vessel conservation and \(q_\ell=\gamma|\partial\Theta|(P-\Pi p_t)\). Its logarithmic pressure singularity requires weighted spaces or singularity treatment. Neither formulation identifies missing lumen connectivity. [Derivation and discretization](#distributed-source-derivation). [Numerical results](#mesh-and-exchange-assumptions).

### G. Explicit lumen bridge

Let \(L_b\) be candidate bridge length, \(R_b(s)\) its radius, \(\widehat P(s)\) bridge pressure, \(q_b\) bridge flow, and \(P_-,P_+\) its endpoint pressures. Define \(a_{v,b}=\pi R_b^4/(8\mu)\). Then

$$
q_b'=0,\qquad q_b=-a_{v,b}\widehat P',\qquad
\widehat P(0)=P_-,\quad\widehat P(L_b)=P_+.
$$

This restores a hypothesized arterial route. Choose candidates using distance, tangent alignment, radius compatibility, and image evidence. Nearby endpoints need not belong to one vessel. A no-leak bridge is an approximation; an oracle repair in L must also restore the removed lateral coupling. [Derivation and discretization](#lumen-bridge-derivation). [Numerical results](#radius-errors-and-repair).

## 2. Weak forms and discretization

A **test function** weights the conservation residual. Multiplying by every admissible test and integrating gives the weak form. Integration by parts transfers one derivative from pressure to the test function. Use \(v_t\) for the bulk test and \(V\) for the vessel test. Tests vanish at prescribed-pressure boundaries; graph tests are continuous at junctions, so conservation cancels internal endpoint terms.

Let \(\partial\Omega_{\rm out}\) be the bulk drainage boundary, \(\gamma_v\geq0\) its hydraulic conductivity, and \(dS\) surface area measure. Define the bulk bilinear form \(\mathfrak a_t\) by

$$
\mathfrak a_t(p_t,v_t)=
\int_\Omega\left(\frac{k_t}{\mu}\nabla p_t\cdot\nabla v_t
+\alpha_vp_tv_t\right)dx
+\int_{\partial\Omega_{\rm out}}\gamma_vp_tv_t\,dS.
$$

This boundary term corresponds to \(\boldsymbol v_t\cdot\boldsymbol n=\gamma_v(p_t-P_{\rm CVP})\). Let \(F(v_t,V)\) collect prescribed inflow and reference-pressure loads, including the \(P_{\rm CVP}\) contributions. A finite-element discretization expands pressure in piecewise-linear basis functions and tests against each basis function, giving one equation per free unknown.

### Resistive graph derivation

Let \(\varrho\) be distance from the tube axis and \(v_s(\varrho)\) axial fluid speed. In this equation only, derivatives of \(v_s\) are with respect to \(\varrho\); \(P'\) remains the axial pressure derivative. Steady axial Stokes balance, regularity at the axis, and no slip give

$$
\frac{\mu}{\varrho}\frac{d}{d\varrho}
\left(\varrho\frac{dv_s}{d\varrho}\right)=P',\qquad
v_s(R)=0
\quad\Longrightarrow\quad
v_s(\varrho)=-\frac{P'}{4\mu}(R^2-\varrho^2),\quad
q=-\frac{\pi R^4}{8\mu}P'.
$$

Conservation makes \(q\) constant on a no-leak edge; integrating \(P'=-q/a_v\) gives its resistance. Let \(V_i\) be the nodal vessel tests. Multiply nodal balances by \(V_i\) and collect each edge once:

$$
\sum_{(i,j)}g_{ij}(P_i-P_j)(V_i-V_j)
+\sum_{i\in\mathcal O}G_{a,i}P_iV_i
=\sum_i f_iV_i+\sum_{i\in\mathcal O}G_{a,i}P_{d,i}V_i.
$$

Let \(B_g\) be oriented node–edge incidence, \(E\) select outlet pressures, \(\Gamma_a=\operatorname{diag}(G_{a,k})\), and \(L\) be the assembled matrix. Then
\(L=B_g\operatorname{diag}(g_e)B_g^T+E^T\Gamma_aE\).
Impose inlet pressure by elimination. For variable radius, use the resistance integral; conductance at an arithmetic mean radius can miss a constriction.

### Lateral average derivation

A slice satisfies \(q'+q_\ell=0\). Integrating the wall law \(v_n=\gamma(P-p_t)\), where \(v_n\) is outward wall-normal fluid speed, around its circumference gives \(q_\ell=\gamma|\partial\Theta|(P-\overline{p_t})\). This assumes circumferentially constant \(\gamma\). LZ also assume approximately uniform cross-sectional vessel pressure, neglect circumferential fluctuation products in bulk coupling, and extend the bulk into the excised thin-vessel volume. Poiseuille supplies the axial blood-flow closure.

Multiply bulk conservation by \(v_t\) and vessel conservation by \(V\), integrate by parts, and add. Bulk exchange contributes \(-q_\ell\overline{v_t}\), vessel exchange \(q_\ell V\), and outlet flow \(q(A_k)V(A_k)\). Hence

$$
\begin{aligned}
\mathfrak a_t(p_t,v_t)+\int_\Lambda a_vP'V'\,ds
&+\int_\Lambda\gamma|\partial\Theta|
(P-\overline{p_t})(V-\overline{v_t})\,ds\\
&+\sum_kG_{a,k}\bigl(P(A_k)-\langle p_t\rangle_k\bigr)
\bigl(V(A_k)-\langle v_t\rangle_k\bigr)=F(v_t,V).
\end{aligned}
$$

For the interstitial alternative, replace the terminal product by \(\sum_kG_{a,k}P(A_k)V(A_k)\) and add \(\sum_kG_{a,k}P_{d,k}V(A_k)\) to \(F\). Terminal blood then does not enter interstitium.

Let \(\phi_i\) and \(\psi_j\) be bulk and vessel basis functions, respectively, and \(K_t,K_v\) their stiffness matrices:

$$
(K_t)_{ik}=\mathfrak a_t(\phi_k,\phi_i),\qquad
(K_v)_{jl}=\int_\Lambda a_v\psi_l'\psi_j'\,ds.
$$

Let \(s_m\) and \(w_m\) be line quadrature points and integration weights. Define the circumference-average matrix \(T_{mi}=\overline{\phi_i}(s_m)\), vessel evaluation matrix \(N_{mj}=\psi_j(s_m)\), and diagonal weight matrix \(W_{mm}=w_m\gamma(s_m)|\partial\Theta(s_m)|\). Let \(A_{ki}=\langle\phi_i\rangle_k\) be terminal averaging and \(D=[-A\ E]\) the terminal pressure-difference matrix. The coupled matrix is

$$
\begin{pmatrix}
K_t+T^TWT&-T^TWN\\
-N^TWT&K_v+N^TWN
\end{pmatrix}+D^T\Gamma_aD.
$$

Using the same averaging and spreading weights gives equal-and-opposite exchange and nonnegative coupling energy. For element length \(h\) and constant \(a_v\), the vessel element contributes \((a_v/h)\begin{pmatrix}1&-1\\-1&1\end{pmatrix}\).

### Terminal exchange derivation

Set \(\gamma=0\) in the preceding weak form. Let \(L_0=B_g\operatorname{diag}(g_e)B_g^T\) be the axial graph matrix without outlet loads. The coupled matrix becomes

$$
\operatorname{diag}(K_t,L_0)+D^T\Gamma_aD.
$$

For a bulk cell \(\Omega_i\), use \(A_{ki}=\int_{\Omega_i}\eta_k^\epsilon\,dx\) both to average pressure and deposit source \(q(A_k)A_{ki}\). Normalize each outlet’s weights to sum to one, including boundary-clipped spheres. On orthogonal isotropic grids, a cell face \(F_{ij}\) of area \(|F_{ij}|\), with cell-center separation \(\ell_{ij}\) and face permeability \(k_{t,ij}\), has conductance \((k_{t,ij}/\mu)|F_{ij}|/\ell_{ij}\). Full anisotropy needs a compatible discretization.

### Perfect-contact derivation

A **Lagrange multiplier** is an extra unknown enforcing a constraint. Here it is \(\lambda\), the exchange needed to match pressures. Let \(X=(p_t,P)\), \(Y=(v_t,V)\), and \(\mathfrak a_0(X,Y)\) contain the uncoupled bulk/axial terms plus terminal coupling. Let \(\zeta\) be an independent scalar test for the pressure constraint. Then

$$
\mathfrak a_0(X,Y)+\int_\Lambda\lambda(V-\overline{v_t})\,ds=F(Y),
\qquad
\int_\Lambda\zeta(P-\overline{p_t})\,ds=0.
$$

Let \(K_0\) be the matrix of \(\mathfrak a_0\), \(B\) the tested pressure-difference matrix, \(\boldsymbol x\) the pressure coefficients, \(\boldsymbol\lambda_h\) multiplier coefficients, and \(\boldsymbol f\) the load vector. The system is

$$
\begin{pmatrix}K_0&B^T\\B&0\end{pmatrix}
\binom{\boldsymbol x}{\boldsymbol\lambda_h}=\binom{\boldsymbol f}{0}.
$$

For finite wall resistance \(R_w=1/(\gamma|\partial\Theta|)\), replace the constraint by \(P-\overline{p_t}-R_w\lambda=0\), adding a negative resistance mass block. Eliminating \(\lambda\) recovers Robin coupling. Choose compatible multiplier/pressure spaces or justified stabilization; arbitrary multiplier discretizations can be unstable.

### Distributed source derivation

Interchanging integrals gives

$$
\int_\Omega v_tS_{\epsilon_\ell}q_\ell\,dx
=\int_\Lambda q_\ell\langle v_t\rangle_{\epsilon_\ell}\,ds.
$$

Thus R has L’s weak form and matrix, replacing circumference averages by volume averages and \(\gamma\) by \(\gamma_{\epsilon_\ell}\). Hold physical support fixed during mesh refinement. For a straight vessel in isotropic bulk, let \(p_t(\varrho)\) denote radial pressure. Radial Darcy flow outside the source gives
\(p_t(R)-p_t(\epsilon_\ell)=\mu q_\ell\log(\epsilon_\ell/R)/(2\pi k_t)\).
This illustrates why changing the sampling scale requires a pressure correction.

For the singular alternative, deposition tests the **centerline** while feedback uses \(\Pi p_t\), giving exchange

$$
\int_\Lambda\gamma|\partial\Theta|(P-\Pi p_t)
(V-v_t|_\Lambda)\,ds.
$$

This differs from L’s symmetric coupling. Since \(|\nabla p_t|\sim1/\varrho\), ordinary gradient energy diverges as \(\int_0^R d\varrho/\varrho\). Weighted trial/test spaces permit this singularity while controlling line testing. Use a justified weighted or singularity-removal discretization, not casual centerline pressure sampling.

### Lumen bridge derivation

Let \(g_b\) be the bridge conductance. Integrating its equations gives

$$
q_b=g_b(P_--P_+),\qquad
g_b=\left(\int_0^{L_b}\frac{8\mu}{\pi R_b(s)^4}\,ds\right)^{-1}.
$$

Let \(V_-,V_+\) be endpoint test values and \(\boldsymbol b_b\) the bridge incidence vector. The weak contribution is \(g_b(P_--P_+)(V_--V_+)\), giving matrix contribution \(g_b\boldsymbol b_b\boldsymbol b_b^T\). This is exact no-leak continuum reduction; a coarse variable-coefficient linear element need not reproduce the resistance integral.

For a local 3D bridge, let \(K_{pp},K_{pc},K_{cp},K_{cc}\) be stiffness blocks partitioned into port (\(p\)) and interior (\(c\)) pressure coefficients. Eliminating interior pressures gives the port matrix \(S=K_{pp}-K_{pc}K_{cc}^{-1}K_{cp}\). It reduces to one conductance only for two lumped ports with no other sources, ports, or drainage. Shared tissue generally produces additional couplings, not one missing arterial edge.

### Conservative flux discretization

This changes numerical representation, not physiology. For T or volume-kernel R, retain \(\boldsymbol v_t\) explicitly and let \(\boldsymbol w\) be a vector test. Multiply Darcy’s law by \(\boldsymbol w\) and integrate by parts:

$$
\int_\Omega\frac{\mu}{k_t}\boldsymbol v_t\cdot\boldsymbol w\,dx
-\int_\Omega p_t\nabla\cdot\boldsymbol w\,dx
=-\int_{\partial\Omega}p_t\boldsymbol w\cdot\boldsymbol n\,dS.
$$

Test conservation directly with cell pressures. Shared face-flux unknowns enforce **outgoing minus incoming flow plus drainage equals sources** cell by cell. Exact L’s surface source instead needs broken/interface flux spaces or a conservative volume projection; ordinary global finite-divergence elements do not directly accommodate the surface measure.

## 3. Sensitivity and the next experiment

Let \(\boldsymbol H\boldsymbol y=\boldsymbol f\) be the assembled system after prescribed pressures are eliminated, with matrix \(\boldsymbol H\), unknown coefficients \(\boldsymbol y\), and loads \(\boldsymbol f\). Let \(\boldsymbol j_T\) be output coefficients and \(J_{T,0}\) the fixed-boundary contribution, so \(J_T=\boldsymbol j_T^T\boldsymbol y+J_{T,0}\). Define the adjoint \(\boldsymbol z\) by \(\boldsymbol H^T\boldsymbol z=\boldsymbol j_T\). For a small parameter change, denoted by \(\delta\), differentiating and substituting gives

$$
\delta J_T=(\delta\boldsymbol j_T)^T\boldsymbol y+\delta J_{T,0}
+\boldsymbol z^T(\delta\boldsymbol f-\delta\boldsymbol H\,\boldsymbol y).
$$

For a constant-radius Poiseuille edge \(e\), let \(R_e\) be its radius, \(\Delta P_e\) its full forward pressure difference, and \(\Delta z_e\) its adjoint difference. With fixed output coefficients and only axial conductance varied,

$$
\frac{\partial J_T}{\partial\log R_e}
=-4g_e\Delta P_e\Delta z_e.
$$

Use zero adjoint values at prescribed-pressure nodes. Geometry-dependent averaging and exchange require additional derivatives. **A small radius derivative does not imply tolerance to deletion:** a low-resistance proximal voxel can disconnect the entire tumor supply.

1. Implement N, L, T, and explicit bridges first. Use C and R as targeted coupling/numerical comparisons.
2. Verify a Y-tree with tumor and healthy outlets under both pressure-driven and flow-driven inflow. Compare intact, gapped, oracle-repaired, estimated-radius-repaired, and tissue-mediated cases, then repeat on one intact in-vivo tree.
3. Produce separate maps of radius sensitivity and finite-gap tumor-flow loss. Evaluate held-out repair error, false connections, conservation, and uncertainty in outlet loads/supply territories. Synthetic intact-tree recovery validates robustness relative to that model, not physiological or dose accuracy.

## 4. Numerical results

The models were run on a 24 mm straight vessel and a Y-shaped vessel with an 11 mm parent and two 13.9 mm daughters. A sphere of radius 4.5 mm around the straight vessel tip or the upper daughter tip defines tumor. Blood enters at 80 mmHg and drains toward 10 mmHg. There is no independent tissue source: arterial exchange supplies the unresolved blood bed. The outer box has zero normal flow. [Geometry](../results/figures/toy_geometry.pdf). The [complete numerical study](../results/sensitivity_results.md) includes all figures and parameters.

### Gaps and tumor delivery

Gaps were removed at nine positions in each vessel, with lengths from 0.25 to 4 mm. A separate experiment removed the center of the bifurcation. Original outlets were retained; newly cut ends were closed. Each comparison used an intact vessel with the same mesh partitions around the affected interval.

A 1 mm gap halfway along the Y tumor branch left essentially zero tumor delivery under vessel-only or terminal-only flow, **0.71%** under wall exchange, **0.63%** under spread exchange, and **77.1%** under equal mean pressures. Thus nonzero flow after a gap can come from the surrounding blood bed rather than a restored lumen. Cutting the healthy branch instead increased terminal-model tumor delivery by **22.2%**: vessels outside the target can matter through competition for inflow.

[Gap position and tumor flow](../results/figures/gap_position_y.pdf) · [Gap length](../results/figures/gap_length_y.pdf) · [Spatial map](../results/figures/gap_effect_map.pdf) · [Bifurcation gap](../results/figures/bifurcation_gap.pdf)

### Radius errors and repair

Radius changes were applied to a 3 mm interval. Halving the radius halfway along the Y tumor branch left **46.1%** of the original tumor flow under Poiseuille resistance, compared with **81.0%** when lumen permeability was held fixed. This isolates the different radius dependence; both laws were matched on the intact vessel.

A bridge restored a 2 mm missing interval, and its radius was varied from half to one and a half times the original radius. Using the original radius recovered the matched intact solution; narrower bridges restored connectivity but could still restrict delivery. Coupled bridges restored their exchange terms as well as axial conductance.

[Radius perturbations](../results/figures/radius_change_y.pdf) · [Bridge radius and recovered flow](../results/figures/bridge_radius_y.pdf) · [Pressure before and after repair](../results/figures/pressure_along_straight_vessel.pdf)

### Mesh and exchange assumptions

Tissue spacing was varied from 5 to 1.5 mm, and vessel resolution was varied separately while physical gap boundaries and source widths stayed fixed. Intact Y tumor flow under wall exchange changed from **2.576 to 2.552 mL min⁻¹** over that tissue refinement. The figures show absolute flow alongside the fraction remaining after a gap, so a stable fraction cannot conceal a changing intact prediction. Equal-pressure constraints used fixed 4 mm wall intervals, separated across a gap. The conservative volume formulation used shared face fluxes for terminal and spread exchange.

Terminal source radius, wall conductivity, and spread radius were also varied. These change different assumptions: where blood is assigned, how readily it exchanges, and how pressure is averaged. Wider sources can cross the tumor boundary. A false connection between daughters with unequal outlet loads could redistribute tumor flow even though no vessel was deleted.

[Tissue resolution](../results/figures/tissue_resolution.pdf) · [Vessel resolution](../results/figures/vessel_resolution.pdf) · [Conservative flux formulation](../results/figures/conservative_flux.pdf) · [Source radius](../results/figures/terminal_source_radius.pdf) · [Wall conductivity](../results/figures/wall_conductivity.pdf) · [Spread radius](../results/figures/spread_radius.pdf) · [False connection](../results/figures/extra_connection.pdf)

These results support separate measures of small radius errors, finite connectivity loss, and recovery after repair. They quantify sensitivity of the toy models; they do not establish physiological accuracy or reactive bolus delivery.

## References

- Laurino & Zunino (2019), [doi:10.1051/m2an/2019042](https://doi.org/10.1051/m2an/2019042). Finite-radius reduction: equations (2.6), (2.8), (3.2); surface-source strong form: (3.18). [Local paper](../papers/laurino-zunino-2019-derivation-analysis-coupled-pdes.pdf).
- D’Angelo & Quarteroni (2008), [doi:10.1142/S0218202508003108](https://doi.org/10.1142/S0218202508003108). Singular line-source coupling and weighted analysis. [Local paper](../papers/dangelo-quateroni-2008-on-the-coupling-1d-3d.pdf).
- Kuchta, Laurino, Mardal & Zunino (2021), [doi:10.1137/20M1329664](https://doi.org/10.1137/20M1329664). Multiplier coupling and discretization stability. [Local paper](../papers/miroslav-2021-analysis-approximation-mixed-dimensional.pdf).
- Amare et al. (2026), [doi:10.64898/2026.07.12.738114](https://doi.org/10.64898/2026.07.12.738114), bioRxiv preprint. Lateral versus terminal/SOI exchange; its LAM already includes terminal coupling. [Local paper](../papers/amare-2026-quantitative-comparison-3d-1d.pdf).
- Koch, Schneider, Helmig & Jenny (2020), [doi:10.1016/j.jcp.2020.109370](https://doi.org/10.1016/j.jcp.2020.109370). Distributed-source tissue perfusion and exchange-flux correction.

The notation follows Amare et al. where practical; the equations retain consistent flow signs, total terminal conductances, and normalized exchange operators rather than copying every printed coefficient convention.
