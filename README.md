# Topology Optimization with Julia

Code for 2-D topology optimization with Julia (https://julialang.org/).
The code translates in Julia the 88 line MATLAB code (<a href="https://www.topopt.mek.dtu.dk/Apps-and-software/Efficient-topology-optimization-in-MATLAB">top88.m</a>) for 2-D topology optimization.
Details of the MATLAB implementation are discussed in the paper:  
"Efficient topology optimization in MATLAB using 88 lines of code, E. Andreassen, A. Clausen, M. Schevenels, B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, Volume 43, Issue 1, p.1 - 16, (2011)."  


To use the code, run in Julia the file ```TopOpt_Julia.jl```.   

The main function is defined as: TopOpt(nelx, nely, rmin, volfrac, penal, ft, alg), where:
<ul>
  <li><i>nelx</i> is the number of finite elements in the x direction;</li> 
  <li><i>nely</i> is the number of finite elements in the y direction;</li> 
  <li><i>rmin</i> is the filter radius;</li>
  <li><i>volfrac</i> is the solid volume fraction allowed;</li> 
  <li><i>penal</i> is the parameter that defines the SIMP interpolation scheme;</li>
  <li><i>ft</i> is the type of filtering, ft=1 activates the sensitivity filter, ft=2 activates the density filter.</li> 
  <li><i>alg</i> is the optimization algorithm chosen, alg="OC" uses the optimality criteria method, alg="MMA" uses the method of moving asymptotes.</li> 
</ul>

The Julia code provided in this respository uses an implementation of the Method of Moving Asymptotes (Svanberg, 1987), which is based on the GCMMA-MMA-code originally written by Prof. Krister Svanberg for MATLAB. The original code can be obtained from http://www.smoptit.se/ and is distributed under the GNU General Public License.   

For instructions on using and citing the GCMMA-MMA-code originally written by Prof. Krister Svanberg, please visit: https://github.com/pollinico/GCMMA-MMA-Julia.   

Extensions and improvements are welcome and strongly encouraged :smiley:   
  
  
