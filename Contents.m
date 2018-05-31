% TEMINV1D
%
% Files
%   correctRampTime        - add ramp time distortion to transient
%   driverInversion        - example driver for synthetic data
%   driverPROTEM           - example driver for PROTEM data
%   getHankelFC            - provides Fast Hankel transform coefficients
%   getJ                   - calculate Jacobian for current model
%   getVMDLayeredDownward  - VMD evaluation of solution in z > 0
%   getVMDLayeredHarmonic  - VMD in the frequency domain
%   getVMDLayeredTransient - VMD in the time domain
%   getVMDLayeredUpward    - VMD evaluation of layer admittances
%   interp_transient       - interpolate PROTEM onto log-equispaced times
%   minmax                 - evaluate smallest and largest elements
%   plot_model             - plot model
%   plotTransient          - plot transient
%   PROTEM_1D_Inversion    - PROTEM_1D_Inversion inversion of PROTEM data
%   simulate_PROTEM        - simulatePROTEM(t, r, rho, thk, t0)
%   taylortest             - Taylor test for jacobian
%   TEM_1D_Inversion       - TEM_1D_Inversion inversion of synthetic TEM data
