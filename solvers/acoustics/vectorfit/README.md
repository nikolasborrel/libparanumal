# Vector fitting for locally-reacting boundary conditions

`vectorfitDriverLR.m` produces the coefficient files consumed by the acoustics
solver setting `[LR VECTORFIT FILE]`, i.e. boundary flag 3, the locally-reacting
(LR) frequency-dependent impedance boundary condition. It comes from
[dtu-act/libparanumal](https://github.com/dtu-act/libparanumal), MIT licensed.

`../data/LRDATA14.dat` was generated with this script, and the parameters it
currently carries are the ones recorded in that file's footer: a 5 cm porous
layer of flow resistivity 47700 Pa·s/m², rigidly backed, fitted with 14 poles
over 50–2000 Hz in air.

## Theory

A locally-reacting surface is one where the wall responds only to the local
pressure, so the wall-normal particle velocity follows the surface admittance

    v_n(ω) = Y(ω) p(ω),        Y = 1/Z

The solver needs this relation in the time domain, but a convolution against
`Y(t)` at every boundary node is not affordable. Vector fitting replaces `Y`
with a rational approximation

    Y(s) ≈ Y_inf + Σ_i  A_i/(s + λ_i)  +  Σ_k  2·Re[ (B_k + i·C_k)/(s + α_k + i·β_k) ]

with real poles `λ_i` and complex-conjugate pole pairs `α_k ± i·β_k`. Each term
is a first-order ODE, so the boundary condition becomes a small set of
accumulator states `φ` per boundary node, integrated alongside the field:

    φ'_i = -λ_i φ_i + p                              (real pole)
    φ'_k = -α_k φ_k - β_k ψ_k + p,  ψ'_k = -α_k ψ_k + β_k φ_k   (complex pair)
    v_n  = Y_inf p + Σ_i A_i φ_i + Σ_k 2(B_k φ_k + C_k ψ_k)

`Y_inf` is the instantaneous part of the response: with all residues zero the
condition collapses to the frequency-independent impedance BC with `Z = 1/Y_inf`.

The material model in the script is the Miki/Delany-Bazley empirical fit for a
porous absorber of flow resistivity `sigma` and thickness `dmat`, optionally
backed by an air cavity of depth `d0`, evaluated at incidence angle `thetai`.

Because the poles are decay rates in 1/s, the fit is tied to the medium it was
made for (`rho`, `c` in the script) and to the frequency range `f_range`. Outside
that range the rational function is an extrapolation, and the accumulator ODEs
are integrated explicitly, so `max|pole|·dt` must stay inside the LSERK4
stability region — stiff fits need a finer mesh or a smaller CFL number.

## Usage

1. Install the [Matrix Fitting Toolbox](https://www.sintef.no/projectweb/vectorfitting/downloads/)
   and put it on the MATLAB path (`VFdriver` and `RPdriver`).
2. Set the material (`sigma`, `dmat`, `d0`), the medium (`rho`, `c`), the fitting
   range `f_range`, the pole count `Npoles` (must be even) and `fileName`.
3. Run the script. It fits `Y`, enforces passivity with `RPdriver`, and writes
   the coefficient file.
4. Point the solver at it:

       [LR VECTORFIT FILE]
       data/LRDATA14.dat

   and give the LR faces boundary flag 3, as in
   `../setups/setupRoomLRTet3D.rc`.

Check the plotted fit against the reference admittance before using it. The wall
is passive only where `Re{Y(ω)} >= 0`; `RPdriver` enforces that, and a fit that
loses it feeds energy into the room instead of absorbing it.

## File format

    Npoles  NRealPoles  NImagPoles
    A       (NRealPoles values, real-pole residues)
    B       (NImagPoles values, real part of the complex-pair residues)
    C       (NImagPoles values, imaginary part)
    lambda  (NRealPoles values, real poles)
    alpha   (NImagPoles values, real part of the complex poles)
    beta    (NImagPoles values, imaginary part)
    Yinf    (1 value)
    -----
    material parameters, ignored by the solver

`Npoles = NRealPoles + 2*NImagPoles`, and the poles are stored negated, so all
values of `lambda` and `alpha` are positive for a stable fit.
