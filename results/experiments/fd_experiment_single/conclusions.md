# Single Square Source Experiment Conclusions

This experiment uses the square heat source:

$$
f=100 \quad \text{for } 0.4<x<0.6,\ 0.4<y<0.6.
$$

## Solver Accuracy

The relative residuals are near machine precision:

```text
mesh study:     about 3e-15 to 2e-13
strength study: about 5e-14
```

This confirms that the computed solution satisfies the discrete system $AU=b$ very accurately. The observed temperature trends are not caused by linear-solver error.

## Source Strength Linearity

At fixed mesh size $n=64$, the maximum temperature scales linearly with source strength:

```text
strength 25  -> t_max 0.3273
strength 50  -> t_max 0.6546
strength 100 -> t_max 1.3092
strength 200 -> t_max 2.6183
strength 500 -> t_max 6.5458
```

This matches the linearity of the Poisson equation. If $f$ is multiplied by $\alpha$, then $T$ is also multiplied by $\alpha$.

## Mesh Study

The mesh convergence is reasonable but not monotone:

```text
n=16  -> t_max 1.656
n=24  -> t_max 0.923
n=32  -> t_max 1.128
n=48  -> t_max 1.346
n=64  -> t_max 1.309
n=96  -> t_max 1.368
n=128 -> t_max 1.322
```

The non-monotone behavior is expected because the square source has a sharp cutoff. On coarse grids, the number of grid points inside the source region changes unevenly. This can be seen from the discrete source integral:

```text
n=16  -> source_integral 5.536
n=24  -> source_integral 2.560
n=32  -> source_integral 3.306
n=48  -> source_integral 4.165
n=64  -> source_integral 4.000
n=96  -> source_integral 4.251
n=128 -> source_integral 4.062
```

The expected continuous source integral is:

$$
100\times0.2\times0.2=4.
$$

By $n=64$ and above, the discrete source integral is close to 4, and $T_{\max}$ settles around order 1.3.

## Main Takeaways

- The solver residual is very small, so the linear solves are reliable.
- The temperature response is linear in source strength.
- The square-source mesh study is bumpy because the heat source is discontinuous.
- Gaussian sources are better for smooth mesh-convergence tests.

