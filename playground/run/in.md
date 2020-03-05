units		metal
atom_style	atomic

neighbor	0.3 bin
neigh_modify	 every 2 delay 10 check yes

variable	a equal 4.062
lattice		fcc $a 
region		box block 0 5 0 5 0 5 
create_box	1 box
create_atoms	1 box

mass		1 26.98

group		Al type 1

velocity	all create 600.03 8728
timestep	0.0005

#pair_style	adp 
#pair_coeff	* * ./Si_Au_Al.gs.mod.5.adp.txt Al
pair_style 	pace
pair_coeff 	* * Al.pbe.ace Al

dump		snap all cfg 10 ./dump.fcc.*.cfg mass type xs ys zs fx fy fz
dump_modify	snap element Al

min_style	cg
minimize	1.0e-10 1.0e-10 1000 1000

fix		1 all nve


thermo 		1
thermo_style 	custom step etotal pe ke pxx pyy pzz press temp vol

run		100

