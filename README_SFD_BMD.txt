VM Mechanical Engineering Virtual Laboratory – SFD and BMD

FILES
1. vm_sfd_bmd_virtual_lab.html – browser-based step-by-step virtual laboratory.
2. vm_sfd_bmd_solver.m – MATLAB Euler-Bernoulli beam finite-element solver and plotting script.

RUN THE WEB LAB
- Keep the HTML and MATLAB files in the same folder.
- Open vm_sfd_bmd_virtual_lab.html in Chrome, Edge or Firefox.
- Select beam type, enter span and section properties.
- Add one or more loads. Magnitudes are positive; vertical loads act downward and applied moments act clockwise.
- Press Prepare experiment, then Next Step. Previous and Reset are available.

BEAM TYPES
- Simply supported
- Cantilever
- Overhanging (enter both support coordinates)
- Fixed at both ends

LOAD TYPES
- Point load (kN)
- UDL (kN/m)
- UVL (linearly increasing from zero at its start to the entered intensity at its end)
- Applied moment (kN.m, clockwise)

MATLAB
- Open vm_sfd_bmd_solver.m in MATLAB.
- Edit USER INPUTS near the top of the script.
- Run the script.
- Units: m, kN, kN.m, GPa and m^4.

MODELLING NOTES / LIMITATIONS
- The web and MATLAB solvers use a linear Euler-Bernoulli beam finite-element model.
- The diagrams are reconstructed from calculated support reactions and applied loads.
- The web plot samples the beam at finite intervals, so concentrated-load jumps may appear as short vertical transitions.
- UVL is assumed to vary linearly from zero at its start to the specified maximum at its end.
- The fixed and overhanging cases are statically indeterminate and are solved using the beam stiffness method.
- The script is an educational tool. Verify unusual loading/support combinations independently before using results for design.
