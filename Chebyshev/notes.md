Chebyshev
========================================

Goal
----------------------------------------

- Quasinormal modes!
- this time not with any package, but from scratch!

Aufpassen
----------------------------------------

- Indices fangen in Julia mit 1 an, aber in den Formeln mit 0


ToDo
----------------------------------------


- Chebyshev grid ✅
- Chebyshev polynomials ✅
- Cardinal Functions ✅
- check if cardinal functions are correct
- get coefficients for a function ✅
- plot function to some order ✅
- arbitrary intervals
- get derivative matrices

Documentation
----------------------------------------
- grid erstellen mit `chebGrid(N)`
    * Chebyshev grid mit N+1 Punkten
    * zwischen -1 und 1
- Chebyshev polynomials `chebT(i)(x)`
- Chebyshev cardinal functions for that grid with `chebCardinal(j, grid)`
- coefficients of a function foo: `coeffs = foo.(grid)`
- to see Chebyshev interpolation: `coeffsToFunction(coeffs, grid`


Notes
----------------------------------------

- `@show i`
