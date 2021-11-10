# Automatically find spring-mass-damper behavior to match observations using adjoint-based sensitivity: Part I Linear system
matlab code which solves the adjoint equations for spring mass damper with generic time dependent mass, spring, damper

Executing the matlab code will generate 
1) Observational data v. initial guess and its adjoint solution
2) Adjoint Error versus stepsize (alpha) compared to finite different approximations
3) Cost during optimization
4) y position of spring mass damper systems during optimization
5) Comparison between observation, initial guess, optimized guess
6) The optimized k,m,c versus the `unknown k-m-c'

The directory also contains a copy of data-expected which provides you `what you should see' if .m is run correctly
