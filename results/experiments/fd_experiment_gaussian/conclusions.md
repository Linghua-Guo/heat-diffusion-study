# Gaussian Source Experiment Conclusions

This experiment uses a centered Gaussian heat source:

$$
f(x,y)=A\exp\left(-\frac{(x-0.5)^2+(y-0.5)^2}{2\sigma^2}\right),
\quad \sigma=0.06.
$$

## Solver Accuracy

The relative residuals are near machine precision:

```text
mesh study:     about 4e-15 to 9e-13
strength study: about 5e-14
```

This confirms that the computed solution satisfies the discrete linear system $AU=b$ very accurately.

## Mesh Convergence

The Gaussian source gives a much smoother mesh convergence trend than the square source:

```text
n=16  -> t_max 0.7464
n=24  -> t_max 0.7592
n=32  -> t_max 0.7638
n=48  -> t_max 0.7670
n=64  -> t_max 0.7682
n=96  -> t_max 0.7690
n=128 -> t_max 0.7693
n=160 -> t_max 0.7695
n=192 -> t_max 0.7695
n=256 -> t_max 0.7696
```

The result is effectively converged by the largest meshes:

$$
T_{\max}\approx 0.7696.
$$

The change from $n=192$ to $n=256$ is only about $7.5\times10^{-5}$, or roughly $0.01\%$ relative change.

This behavior is smoother because the Gaussian source is continuous and does not depend on a sharp geometric cutoff.

## Source Integral Stability

The discrete source integral is essentially constant across mesh sizes:

```text
n=16  -> source_integral 2.2619466997
n=24  -> source_integral 2.2619467106
n=32  -> source_integral 2.2619467106
n=64  -> source_integral 2.2619467106
n=128 -> source_integral 2.2619467106
n=256 -> source_integral 2.2619467106
```

This is an important difference from the square source case, where the active source area jumps on coarse grids.

## Source Strength Linearity

At fixed mesh size $n=64$, the maximum temperature scales linearly with source strength:

```text
strength 25  -> t_max 0.1920
strength 50  -> t_max 0.3841
strength 100 -> t_max 0.7682
strength 200 -> t_max 1.5364
strength 500 -> t_max 3.8409
```

This confirms the linearity of the Poisson equation: multiplying $f$ by $\alpha$ multiplies $T$ by $\alpha$.

## Main Takeaways

- The solver residual is very small, so the linear solves are reliable.
- The Gaussian source gives smooth mesh convergence.
- The maximum temperature is essentially converged at the largest tested meshes.
- The source integral is stable even on coarse grids.
- Gaussian sources are better than discontinuous square sources for clean convergence studies.
