Measuring performance
=====================

Before trying to make your simulation run faster, you should understand
how it currently performs and where the bottlenecks are.  We generally
distinguish between serial performance (how fast can a single process do
the calculations?) and parallel efficiency (how much faster does a
calculation get by using more processes?).  There are many factors
affecting either and below are some lists discussing some commonly
known but also some less known factors.

Factors affecting serial performance (in no specific order):

* CPU hardware: clock rate, cache sizes, CPU architecture (instructions
  per clock, vectorization support, fused multiply-add support and more)
* RAM speed and number of channels that the CPU can use to access RAM
* Cooling: CPUs can change the CPU clock based on thermal load, thus the
  degree of cooling can affect the speed of a CPU.  Sometimes even the
  temperature of neighboring compute nodes in a cluster can make a
  difference.
* Compiler optimization: most of LAMMPS is written to be easy to modify
  and thus compiler optimization can speed up calculations. However, too
  aggressive compiler optimization can produce incorrect results or
  crashes (during compilation or at runtime).
* Source code improvements: styles in the OPT, OPENMP, and INTEL package
  can be faster than their base implementation due to improved data
  access patterns, cache efficiency, or vectorization. Compiler optimization
  is required to take full advantage of these.
* Number and kind of fixes, computes, or variables used during a simulation,
  especially if they result in collective communication operations
* Pair style cutoffs and system density: calculations get slower the more
  neighbors are in the neighbor list and thus for which interactions need
  to be computed.  Force fields with pair styles that compute interactions
  between triples or quadruples of atoms or that use embedding energies or
  charge equilibration will need to walk the neighbor lists multiple times.
* Neighbor list settings: tradeoff between neighbor list skin (larger
  skin = more neighbors, more distances to compute before applying the
  cutoff) and frequency of neighbor list builds (larger skin = fewer
  neighbor list builds).
* Proximity of per-atom data in physical memory that for atoms that are
  close in space improves cache efficiency (thus LAMMPS will by default
  sort atoms in local storage accordingly)
* Using r-RESPA multi-timestepping or a SHAKE or RATTLE fix to constrain
  bonds with higher-frequency vibrations may allow a larger (outer) timestep
  and thus fewer force evaluations (usually the most time consuming step in
  MD) for the same simulated time (with some tradeoff in accuracy).

Factors affecting parallel efficiency (in no specific order):

* Bandwidth and latency of communication between processes. This can vary a
  lot between processes on the same CPU or physical node and processes
  on different physical nodes and there vary between different
  communication technologies (like Ethernet or InfiniBand or other
  high-speed interconnects)
* Frequency and complexity of communication patterns required
* Number of "work units" (usually correlated with the number of atoms
  and choice of force field) per MPI-process required for one time step
  (if this number becomes too small, the cost of communication becomes
  dominant).
* Choice of parallelization method (MPI-only, OpenMP-only, MPI+OpenMP,
  MPI+GPU, MPI+GPU+OpenMP)
* Algorithmic complexity of the chosen force field (pair-wise vs. many-body
  potential, Ewald vs. PPPM vs. (compensated or smoothed) cutoff-Coulomb)
* Communication cutoff: a larger cutoff results in more ghost atoms and
  thus more data that needs to be communicated
* Frequency of neighbor list builds: during a neighbor list build the
  domain decomposition is updated and the list of ghost atoms rebuilt
  which requires multiple global communication steps
* FFT-grid settings and number of MPI processes for kspace style PPPM:
  PPPM uses parallel 3d FFTs which will drop much faster in parallel
  efficiency with respect to the number of MPI processes than other
  parts of the force computation.  Thus using MPI+OpenMP parallelization
  or :doc:`run style verlet/split <run_style>` can improve parallel
  efficiency by limiting the number of MPI processes used for the FFTs.
* Load (im-)balance: LAMMPS' domain decomposition assumes that atoms are
  evenly distributed across the entire simulation box. If there are
  areas of vacuum, this may lead to different amounts of work for
  different MPI processes. Using the :doc:`processors command
  <processors>` to change the spatial decomposition, or MPI+OpenMP
  parallelization instead of only-MPI to have larger sub-domains, or the
  (fix) balance command (without or with switching to communication style
  tiled) to change the sub-domain volumes are all methods that
  can help to avoid load imbalances.

The best way to do this is run the your system (actual number of
atoms) for a modest number of timesteps (say 100 steps) on several
different processor counts, including a single processor if possible.
Do this for an equilibrium version of your system, so that the
100-step timings are representative of a much longer run.  There is
typically no need to run for 1000s of timesteps to get accurate
timings; you can simply extrapolate from short runs.

For the set of runs, look at the timing data printed to the screen and
log file at the end of each LAMMPS run.  The
:doc:`screen and logfile output <Run_output>` page gives an overview.

Running on one (or a few processors) should give a good estimate of
the serial performance and what portions of the timestep are taking
the most time.  Running the same problem on a few different processor
counts should give an estimate of parallel scalability.  I.e. if the
simulation runs 16x faster on 16 processors, its 100% parallel
efficient; if it runs 8x faster on 16 processors, it's 50% efficient.

The most important data to look at in the timing info is the timing
breakdown and relative percentages.  For example, trying different
options for speeding up the long-range solvers will have little impact
if they only consume 10% of the run time.  If the pairwise time is
dominating, you may want to look at GPU or OMP versions of the pair
style, as discussed below.  Comparing how the percentages change as
you increase the processor count gives you a sense of how different
operations within the timestep are scaling.  Note that if you are
running with a Kspace solver, there is additional output on the
breakdown of the Kspace time.  For PPPM, this includes the fraction
spent on FFTs, which can be communication intensive.

Another important detail in the timing info are the histograms of
atoms counts and neighbor counts.  If these vary widely across
processors, you have a load-imbalance issue.  This often results in
inaccurate relative timing data, because processors have to wait when
communication occurs for other processors to catch up.  Thus the
reported times for "Communication" or "Other" may be higher than they
really are, due to load-imbalance.  If this is an issue, you can
use the :doc:`timer sync <timer>` command to obtain synchronized timings.
