# MATLAB Codes for Numerical Examples in “Barycentric Rational Collocation for Volterra Integral Equations with Diagonal and Boundary Singularities”

This repository contains MATLAB implementations of the three numerical examples presented in the paper:

> Yuee Zhong, Kelong Zhao, Shuhuang Xiang. *Barycentric rational collocation method for second-kind Volterra integral equations with product kernels having diagonal and boundary singularities*. (2026)

The codes reproduce the computed solutions for the Volterra integral equations of the second kind with weakly singular kernels. The numerical method is the **barycentric rational collocation scheme** based on the **2nd‑tapered exponentially clustered poles** (BIEP) and **Gauss–Laguerre quadrature**, which achieves the optimal root‑exponential convergence rate \(O(e^{-2\pi\sqrt{\alpha n}})\).

## Requirements

- **MATLAB** (R2018b or later recommended)
- **Chebfun** toolbox – required for computing Gauss–Laguerre nodes and weights with high accuracy.  
  Download and install Chebfun from [https://www.chebfun.org](https://www.chebfun.org).

## Usage

1. **Install Chebfun** and make sure it is on your MATLAB path.
2. Run any example script directly in the MATLAB command window:
   ```matlab
   example1
   example2
   example3
