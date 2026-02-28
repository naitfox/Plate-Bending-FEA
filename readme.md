# Plate Bending (FEA)

## Problem Statement

A simply supported square plate of side **$a = 2\text{ m}$** and thickness **$\Delta z = 0.01\text{ m}$** is subjected to a uniform load **$q = 33.6\text{ kN/m}^2$**.

**Given Parameters:**

* Elastic Modulus: $E = 2 \times 10^{11}\text{ Pa}$
* Poisson's Ratio: $\sigma = 0.3$
* Boundary Conditions: $w = 0$ and $\frac{\partial^2 w}{\partial n^2} = 0$ on all edges (simply supported).
* Discretization for test run: $\Delta x = \Delta y = 0.5$

## Governing Equation

The deflection $w(x, y)$ is given by:

$$\frac{\partial^4 w}{\partial x^4} + 2\frac{\partial^4 w}{\partial x^2 \partial y^2} + \frac{\partial^4 w}{\partial y^4} = q/D$$

Where **$D$** is the flexural rigidity:


$$D = \frac{E \Delta z^3}{12(1 - \sigma^2)}$$

We multiply by a weighting function '$v$' and integrate over the plate area. Using partial integration twice, we arrive at the **weak form**:

$$\int_A D(\nabla^2 w)(\nabla^2 v) \, dA = \int_A pv \, dA$$

---

## Finite Element Formulation

According to **Kirchhoff's plate theory**, we need $C^1$ continuity because the weak form has 2nd derivatives of $w$.

### Element Type

We will use **rectangular 4-noded elements** as they conform with the plate geometry. For each element, we approximate the deflection as:

$$w(x, y) = N_1(x, y)w_1 + N_2(x, y)w_2 + N_3(x, y)w_3 + N_4(x, y)w_4$$

Where each node has **3 Degrees of Freedom (DOFs)**: $\{w, \theta_x, \theta_y\}$.

* $\theta_x = \partial w / \partial y$
* $\theta_y = -\partial w / \partial x$

The displacement matrix is defined as:


$$[w_i, \theta_{xi}, \theta_{yi}]^T = \{d_i\} \quad \text{for } i=1, 2, 3, 4$$

### Polynomial Basis

Chosen from Pascal's Triangle for this application:


$$w(x, y) = a_1 + a_2x + a_3y + a_4x^2 + a_5xy + a_6y^2 + a_7x^3 + a_8x^2y + a_9xy^2 + a_{10}y^3 + a_{11}x^3y + a_{12}xy^3$$

Also, $\{d_i, d_j, d_k, d_l\}^T = [A][a]$
Where $[a]$ is the coefficient matrix: $[a] = [a_1, a_2, a_3, \dots, a_{12}]^T$. From this, we can calculate $[A]$.

---

## Matrices and Integration

### Strain Matrix

The strain matrix $\{\varepsilon\}$ is given by:

\left\{ -\frac{\partial^2 w}{\partial x^2}, -\frac{\partial^2 w}{\partial y^2}, -2\frac{\partial^2 w}{\partial x \partial y} \right\}

### Element Stiffness Matrix ($K_e$)

$$K_e = \iint_A B^T D B \, dx dy$$

Where **$D$** is the rigidity matrix:


$$D = \frac{E \Delta z^3}{12(1 - \sigma^2)} \begin{bmatrix} 1 & \sigma & 0 \\ \sigma & 1 & 0 \\ 0 & 0 & \frac{1-\sigma}{2} \end{bmatrix}$$

### Force Matrix ($F_e$)

For uniform load conditions:


$$F_e = \iint_A N^T p(x, y) \, dx dy$$


*(Note: $N$ is the interpolation matrix found for each element)*

---

## Global Assembly and Verification

By summing up all element stiffness matrices with the direct stiffness method and all elemental force matrices, we get the global system:


$$\mathbf{K}d = \mathbf{F}$$

### Discretization

For the test run, taking $\Delta x = \Delta y = 0.5$ results in **16 elements** with a total of **25 nodes**.

### Theoretical Verification

For a simply supported square plate under uniform load, the maximum central deflection is:


$$w_{max} = \frac{\alpha q a^4}{D}$$

Given $\alpha = 0.00406$ (for $\sigma = 0.3$):


$$w_{max} = \frac{0.00406 \times 33.6 \times 10^3 \times 2^4}{1.83 \times 10^4} = 0.119 \text{ m}$$

**Conclusion:** The theoretical maximum central deflection is approximately **$0.12 \text{ m}$**, which can be used to verify the code results.

---

## Authors

Jena