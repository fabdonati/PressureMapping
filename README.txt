Input files
•	bVTK:		Binary mask in .vtk format
•	Velocity: 		Velocity (Ntimeframes, Ncomponents, Nx, Ny, Nz) in mm/s
•	opts.dt:		Time spacing dt
•	opts.rho:		Fluid density in kg/mm3
•	opts.mu:		Fluid dynamic viscosity in kg/mm/s

User defined options 
•	opts.interp2d:	‘nn’ / ‘linear’ to define interpolation type when generating 2D planes to compute fluxes and advective energy component
1.	Use ‘nn’ to have nearest point interpolation
2.	Use ‘linear’ to have linear interpolation
•	opts.stencil:	‘standard’ / ‘filtered’ to define stencil used for FDM scheme.
1.	Use ‘standard’ to apply II order centered scheme standard stencil
2.	Use ‘filtered’ to cull noise effect
•	opts.timescheme:	‘cdt’ / ‘cdt2’ to defined central time derivative scheme.
1.	Use ‘cdt‘ to have standard central differences with error O(dt2)
2.	Use ‘cdt2’ to have standard central differences with error O((dt/2)2)
•	opts.MeshName:	‘filename’ to define name of output quadratic hexahedral mesh for PPE computation on the masked region
•	opts.compliance: ‘0’/ ‘1’
1.	Use ‘0’ if you want to exclude walls contribution in the boundary flux term
2.	Use ‘1’ if you want to include walls contribution in the boundary flux term
•	opts.AdveFilter:	any number
1.	Use ‘0’ if you want to compute the advective energy without filtering
2.	Use ‘1’ to define error and advective component of the same order of magnitude,  ‘10’ for 1 order of magnitude, ‘100’ for 2 orders of magnitude etc.
•	opts.AdveAVG:	‘0’ / ‘1’
1.	Use ‘0’ if you want to get the advective component as it is
2.	Use ‘1’ if you want to compute averaged value over consecutive planes for the advective component

When asked, type number of points (planes to compute advective component and boundary inlet and outlet fluxes): e.g. ‘3’ if you want to get 3 points at the inlet and 3 points at the outlet. Then, define 1 inlet point and 1 outlet point and the code select the neighbouring points automatically.

Output files
•	werp.pdo:		pressure drop obtained using outlet flux as Lambda 
•	werp.pdi:		pressure drop obtained using inlet flux as Lambda
•	werp.kpdo:	kinetic drop obtained using outlet flux as Lambda
•	werp.kpdi:		kinetic drop obtained using inlet flux as Lambda
•	werp.apdo:	advective drop obtained using outlet flux as Lambda
•	werp.apdi:		advective drop obtained using inlet flux as Lambda
•	werp.vpdo:	viscous drop obtained using outlet flux as Lambda
•	werp.vpdi:		viscous drop obtained using inlet flux as Lambda

