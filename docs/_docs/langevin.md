<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
});
</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

<!-- # Deterministic and Stochastic Dynamics -->
## Langevin Dynamics

Conducting hybrid molecular dynamics – Monte Carlo (MDMC) schemes is possible by addition of the Langevin (pseudo) move in the `moves` section. Example:

~~~ yaml
moves:
    - langevin_dynamics:
        nsteps: 25
        integrator: {time_step: 0.005, friction: 5}
    - ...
~~~

`langevin`      | Description
--------------- | --------------------------------------------
`nsteps`        | Number of time iterations for each Langevin dynamics event
`integrator`    | Object with the following two keywords:
`time_step`     | Time step for integration <!--∆𝑡--> $\Delta t$ (ps)
`friction`      | Friction coefficient <!--𝛾--> $\gamma$ (1/ps)

This move will solve the Langevin equation for the particles in the system on the form

$$
d\begin{bmatrix} \mathbf{q} \\ \mathbf{p} \end{bmatrix} =
\underbrace{\begin{bmatrix} M^{-1}\mathbf{q} \\ 0 \end{bmatrix} \mathrm{d}t}\_{\text{A}} +
\underbrace{\begin{bmatrix} 0 \\ -\nabla U(\mathbf{q}) \end{bmatrix}\mathrm{d}t}\_{\text{B}} +
\underbrace{\begin{bmatrix} 0 \\ -\gamma \mathbf{p} + \sqrt{2\text{k}\_{\mathrm{b}}T\gamma} \sqrt{M} ~\mathrm{d}\mathbf{W} \end{bmatrix}}\_{\text{O}}
$$

Where A, B, and O makes up the terms for solving the Langevin equation, which can be individually solved to obtain a trajectory given by
\begin{aligned}
\varphi^{\mathrm{A}}(\mathbf{q}, \mathbf{p}) &= \left(\mathbf{q} + \Delta t \sqrt{M}~ \mathbf{p}, \mathbf{p}\right) \\
\varphi^{\mathrm{B}}(\mathbf{q}, \mathbf{p}) &= \left(\mathbf{q}, \mathbf{p} - \Delta t \nabla U(\mathbf{q})\right) \\
\varphi^{\mathrm{O}}(\mathbf{q}, \mathbf{p}) &= \left(\mathbf{q}, e^{-\gamma \Delta t}\mathbf{p} + \sqrt{\mathrm{k}\_{\mathrm{B}}T (1 - e^{-2\gamma \Delta t})} \sqrt{M} \mathbf{R} \right)
\end{aligned}

We currently use the splitting scheme "BAOAB" ([Symmetric Langevin Velocity-Verlet](http://doi.org/ggrnfs)) since it is less errorprone with increasing timestep [Leimkuhler & Matthews, pp. 279-281](http://doi.org/dx7v).


<!--

The keyword `splitting` thus refers to the string constructed from the labels A, B and O. In the current implementation a selection of splitting schemes are available

`splitting`     | Description
--------------- | --------------------------------------------
ABAO            | [Geometric Langevin Position-Verlet](http://doi.org/fh8jpp)
BABO            | [Geometric Langevin Velocity-Verlet](http://doi.org/fh8jpp)
ABOBA           | [Symmetric Langevin Position-Verlet](http://doi.org/ggrnfs)
BAOAB           | [Symmetric Langevin Velocity-Verlet](http://doi.org/ggrnfs)
OABAO           | [Bussi-Parrinello Langevin Position-Verlet](http://doi.org/cjpdjx)
OBABO           | [Bussi-Parrinello Langevin Velocity-Verlet](http://doi.org/cjpdjx)

We recommend the usage of the splitting scheme "BAOAB" (Symmetric Langevin Velocity-Verlet) due to it being less errorprone with increasing timestep [Leimkuhler & Matthews, pp. 279-281](http://doi.org/dx7v).

## Molecular dynamics
To be implemented. It is currently possible to sample the NVE ensemble using the Langevin move by setting `friction = 0`.

For examples of the usage of the Langevin move and combined usage with Monte Carlo moves see the [molecular dynamics](https://github.com/mlund/faunus/blob/master/examples/moleculardynamics) directory.
-->
