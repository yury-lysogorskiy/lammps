LAMMPS (4 Feb 2020)
OMP_NUM_THREADS environment is not set. Defaulting to 1 thread. (src/comm.cpp:94)
  using 1 OpenMP thread(s) per MPI task
echo both
units		metal
atom_style	atomic

neighbor	0.3 bin
neigh_modify	 every 2 delay 10 check yes

variable	a equal 3.597
lattice		fcc $a
lattice		fcc 3.597
Lattice spacing in x,y,z = 3.597 3.597 3.597
region		box block 0 2 0 2 0 2
create_box	1 box
Created orthogonal box = (0 0 0) to (7.194 7.194 7.194)
  1 by 1 by 1 MPI processor grid
create_atoms	1 box
Created 32 atoms
  create_atoms CPU = 0.000699997 secs

mass		1 26.98

group		Al type 1
32 atoms in group Al

pair_style 	pace/al
ACE/AL version: 2021.10.25
pair_coeff  * * Cu.yaml Cu.asi Cu
Loading Cu.yaml
Total number of basis functions
	Cu: 15 (r=1) 717 (r>1)
Mapping LAMMPS atom type #1(Cu) -> ACE species type #0
Loading ASI Cu.asi

fix 1 all box/relax iso 0.0 vmax 0.001
min_style	cg
minimize	1.0e-10 1.0e-10 1000 1000
WARNING: Using 'neigh_modify every 1 delay 0 check yes' setting during minimization (src/min.cpp:190)
Neighbor list info ...
  update every 1 steps, delay 0 steps, check yes
  max neighbors/atom: 2000, page size: 100000
  master list distance cutoff = 7.3
  ghost atom cutoff = 7.3
  binsize = 3.65, bins = 2 2 2
  1 neighbor lists, perpetual/occasional/extra = 1 0 0
  (1) pair pace/al, perpetual
      attributes: full, newton on
      pair build: full/bin/atomonly
      stencil: full/bin/3d
      bin: standard
Per MPI rank memory allocation (min/avg/max) = 4.12 | 4.12 | 4.12 Mbytes
Step Temp E_pair E_mol TotEng Press Volume 
       0            0   -118.28129            0   -118.28129    42382.243    372.31566 
      11            0   -118.41313            0   -118.41313 -0.017250551    382.57568 
Loop time of 0.622325 on 1 procs for 11 steps with 32 atoms

100.0% CPU use with 1 MPI tasks x 1 OpenMP threads

Minimization stats:
  Stopping criterion = energy tolerance
  Energy initial, next-to-last, final = 
        -118.281285202     -118.413125729     -118.413125729
  Force two-norm initial, final = 29.5465 1.2246e-05
  Force max component initial, final = 29.5465 1.2246e-05
  Final line search alpha, max atom move = 0.53673 6.57282e-06
  Iterations, force evaluations = 11 13

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0.62188    | 0.62188    | 0.62188    |   0.0 | 99.93
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 8.75e-05   | 8.75e-05   | 8.75e-05   |   0.0 |  0.01
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 0.0003603  |            |       |  0.06

Nlocal:    32 ave 32 max 32 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    1067 ave 1067 max 1067 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
FullNghs:  4480 ave 4480 max 4480 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 4480
Ave neighs/atom = 140
Neighbor list builds = 0
Dangerous builds = 0

velocity        all create 300 8728
timestep        0.0005
fix		2 all nve

thermo 		100
thermo_style 	custom step etotal pe ke pxx pyy pzz press temp vol
print "Run NVE for 1000 steps"
Run NVE for 1000 steps
run		500
Per MPI rank memory allocation (min/avg/max) = 2.995 | 2.995 | 2.995 Mbytes
Step TotEng PotEng KinEng Pxx Pyy Pzz Press Temp Volume 
      11   -117.21101   -118.41313    1.2021193    2891.1072     3876.156    3301.3198    3356.1943          300    382.57568 
     100   -117.21093   -117.71615   0.50522436    10042.501     10407.66    13864.396    11438.186    126.08341    382.57568 
     200   -117.21089   -117.76035   0.54946074    11610.616    11078.137    12379.165    11689.306    137.12301    382.57568 
     300   -117.21089   -117.65097   0.44007917    11120.708    10367.736    15099.204    12195.882    109.82583    382.57568 
     400   -117.21094   -118.08796   0.87701808    7258.3759    7866.4075    7921.5507    7682.1114    218.86797    382.57568 
     500   -117.21088   -117.60905   0.39816243    11200.959    11621.361    15969.034    12930.452    99.365116    382.57568 
     511   -117.21092   -117.81995   0.60903359    8540.4502    9190.1514     12884.01    10204.871    151.98996    382.57568 
Loop time of 21.6831 on 1 procs for 500 steps with 32 atoms

Performance: 0.996 ns/day, 24.092 hours/ns, 23.059 timesteps/s
100.0% CPU use with 1 MPI tasks x 1 OpenMP threads

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 21.675     | 21.675     | 21.675     |   0.0 | 99.96
Neigh   | 0.001698   | 0.001698   | 0.001698   |   0.0 |  0.01
Comm    | 0.0035813  | 0.0035813  | 0.0035813  |   0.0 |  0.02
Output  | 0.00021911 | 0.00021911 | 0.00021911 |   0.0 |  0.00
Modify  | 0.00086093 | 0.00086093 | 0.00086093 |   0.0 |  0.00
Other   |            | 0.001673   |            |       |  0.01

Nlocal:    32 ave 32 max 32 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    910 ave 910 max 910 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
FullNghs:  4480 ave 4480 max 4480 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 4480
Ave neighs/atom = 140
Neighbor list builds = 9
Dangerous builds = 0
print "DONE"
DONE

Total wall time: 0:00:22
