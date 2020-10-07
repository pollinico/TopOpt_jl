# Topology Optimization with Julia

Code for 2-D topology optimization with Julia (https://julialang.org/).
The code translates in Julia the 88 line Matlab code (<a href="https://www.topopt.mek.dtu.dk/Apps-and-software/Efficient-topology-optimization-in-MATLAB">top88.m</a>) for 2-D topology optimization.
Details of the Matlab implementation are discussed in the paper:
Efficient topology optimization in MATLAB using 88 lines of code, E. Andreassen, A. Clausen, M. Schevenels, B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, Volume 43, Issue 1, p.1 - 16, (2011) .   


To use the code, type in Julia:
```
include("TopOpt_Julia.jl")
TopOptOC()
```   
In this way the code runs with the default input values.  
The main function is defined as: TopOptOC(nelx, nely, rmin, volfrac, penal, ft), where:
<ul>
  <li><i>nelx</i> is the number of finite elements in the x direction</li> 
  <li><i>nely</i> is the number of finite elements in the y direction</li> 
  <li><i>rmin</i> is the filter radius</li>
  <li><i>volfrac</i> is the solid volume fraction allowed</li> 
  <li><i>penal</i> is the parameter the defines the SIMP interpolation scheme</li>
  <li><i>ft</i> is the type of filtering, 1 activates the sensitivity filter, 2 activates the density filter.</li> 
</ul>

If you have used the topology optimization Julia code provided in this repository in your research work, or find it useful in any way, please consider to cite:
```
@misc{nipol2020tojl,
  author = {Pollini, Nicol{\`o}},
  title = {A 2-D topology optimization code written in Julia},
  year = {2020},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/pollinico/topopt_jl}},
  }
  ```
