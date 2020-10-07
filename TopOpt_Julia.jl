# Tested and developed in Windows with Julia 1.5.1, October 2020.
# 2-D Topoplogy Optimization in Julia
# IMPORT PACKAGES
using Printf
using SparseArrays
using LinearAlgebra
using Plots: heatmap
# FUNCTIONS
function optimize!(F::Array{Float64,1}, U::Array{Float64,1}, 
                  freedofs::Array{Int64,1}, 
                  penal::Float64, E0::Float64, Emin::Float64, ft::Int64,
                  nelx::Int64, nely::Int64, 
                  iK::Array{Int64,1}, 
                  jK::Array{Int64,1}, 
                  edofMat::Array{Int64,2}, 
                  KE::Array{Float64,2}, 
                  H::SparseMatrixCSC{Float64,Int64}, 
                  Hs::Array{Float64,2}, volfrac::Float64)::Array{Float64,2}
    x = repeat([volfrac],nely,nelx)
    xPhys = copy(x)
    xnew = copy(x)
    loop = 0
    change = 1.0
    sK = zeros(64*nelx*nely)
    K = spzeros(2*(nely+1)*(nelx+1),2*(nely+1)*(nelx+1))
    dc = zeros(nely, nelx)
    dv = zeros(nely, nelx)
    ce = zeros(nely, nelx)
    while change > 0.01
      loop = loop +1
      # FE-ANALYSIS
      sK = reshape(vec(KE)*(Emin.+vec(xPhys)'.^penal*(E0-Emin)),64*nelx*nely)
      K = sparse(iK,jK,sK)
      U[freedofs] .= Symmetric(K[freedofs,freedofs])\F
      # OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
      ce = reshape(sum((U[edofMat]*KE).*U[edofMat], dims=2),nely,nelx)
      c = sum((Emin.+xPhys.^penal.*(E0-Emin)).*ce)
      dc = -penal*(E0-Emin).*xPhys.^(penal-1).*ce
      dv = ones(nely,nelx)
      # FILTERING/MODIFICATION OF SENSITIVITIES
      if ft == 1 # sensitivity filtering
        dc = reshape(H*(vec(x).*vec(dc))./Hs./max.(1e-3,vec(x)),nely,nelx)
      elseif ft == 2 # density filtering
        dc = reshape(H*(vec(dc)./Hs),nely,nelx)
        dv = reshape(H*(vec(dv)./Hs),nely,nelx)
      end
      # OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
      xPhys, xnew = optimalityCriteria!(x, xPhys, xnew, dc, dv, H, Hs, nelx, nely, ft, volfrac)
      change = maximum(abs.(vec(xnew) - vec(x)))
      x = copy(xnew)
      s = @sprintf "It: %i  Obj.: %.2f  Ch: %.3f" loop c change;
      println(s)
    end
    return x
end
function optimalityCriteria!(x::Array{Float64,2}, xPhys::Array{Float64,2}, xnew::Array{Float64,2}, 
                            dc::Array{Float64,2}, dv::Array{Float64,2}, H::SparseMatrixCSC{Float64,Int64}, 
                            Hs::Array{Float64,2}, nelx::Int64, nely::Int64, ft::Int64,
                            volfrac::Float64)::Tuple{Array{Float64,2},Array{Float64,2}}
    l1 = 0.0; l2 = 1.0e9; move = 0.2
    while (l2-l1)/(l1+l2) > 1.0e-3
      lmid = 0.5*(l2+l1)
      xnew = max.(0.0,max.(x.-move,min.(1.0,min.(x.+move,x.*sqrt.(-dc./dv./lmid)))))
      if ft == 1
        xPhys = copy(xnew)
      elseif ft == 2
        xPhys = reshape((H*vec(xnew))./Hs,nely,nelx)
      end
      if (sum(xPhys) > volfrac*nelx*nely) l1 = lmid else l2 = lmid end
    end
    return xPhys, xnew
end
function prepFilter(nelx::Int64, nely::Int64, rmin::Int64)::Tuple{SparseMatrixCSC{Float64,Int64}, Array{Float64,2}}
  iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2)
  jH = ones(size(iH))
  sH = zeros(size(iH))
  k = 0
  for i1 in 1:nelx
    for j1 in 1:nely
      e1 = (i1-1)*nely+j1
      for i2 in max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
        for j2 in max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
          e2 = (i2-1)*nely+j2
          k = k+1
          iH[k] = e1
          jH[k] = e2
          sH[k] = max(0.0,rmin-sqrt((i1-i2)^2+(j1-j2)^2))
        end
      end
    end
  end
  iH = Int.(iH)
  jH = Int.(jH)
  H = sparse(iH,jH,sH)
  Hs = sum(H, dims=2)
  return (H,Hs)
end
################################################################################
# MAIN CODE
function TopOptOC(nelx = 150::Int64, nely = 50::Int64, 
                  rmin = 6::Int64, volfrac = 0.5::Float64, 
                  penal = 3.0::Float64, ft = 1::Int64)
  # TOPOPT SETTINGS
  # nelx, elements in x
  # nely, elements in y
  # rmin, filter radius
  # volfrac, volume fraction
  # penal, SIMP penalty
  # ft, filtering option: 1 = sensitivity filter, 2 density filter
  nele = nelx * nely
  # MATERIAL PROPERTIES
  E0 = 1.0
  Emin = 1.0e-9
  nu = 0.3
  # ELEMENT STIFFNESS MATRIX
  A11 = Float64[12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12.0]
  A12 = Float64[-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6.0]
  B11 = Float64[-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4.0]
  B12 = Float64[2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2.0]
  KE = 1/(1-nu^2)/24 * ([A11 A12; A12' A11] + nu*[B11 B12; B12' B11])
  # CONNECTIVITY
  nodenrs = reshape(1:(1+nelx)*(1+nely), 1+nely, 1+nelx)
  edofVec = reshape(2*nodenrs[1:end-1,1:end-1]+ones(nely,nelx), nelx*nely, 1)
  edofMat = Int.(repeat(edofVec, 1, 8) + repeat([0 1 2*nely .+ [2 3 0 1] -2 -1], nelx*nely, 1))
  iK = reshape(kron(edofMat, ones(8,1))', 64*nelx*nely)
  jK = reshape(kron(edofMat, ones(1,8))', 64*nelx*nely)
  iK = convert(Array{Int64}, iK)
  jK = convert(Array{Int64}, jK)
  # DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
  F = zeros(2*(nely+1)*(nelx+1))
  F[2,1] = -1.0
  U = zeros(2*(nely+1)*(nelx+1))
  fixeddofs = union(collect(1:2:2*(nely+1)), [2*(nelx+1)*(nely+1)])
  alldofs = collect(1:2*(nely+1)*(nelx+1))
  freedofs = setdiff(alldofs,fixeddofs)
  # PREPARE FILTER
  H,Hs = prepFilter(nelx,nely,rmin)
  # INITIALIZE ITERATION
  x = Array{Float64}(undef, nely, nelx)
  x = optimize!(F[freedofs], U, freedofs, penal, E0, Emin, ft, nelx, nely, iK, jK, edofMat, KE, H, Hs, volfrac)
  # PLOT FINAL DESIGN
  heatmap(1.0.-x[end:-1:1,:],  yaxis=false, xaxis=false, legend = :none, 
  color = :greys, grid=false, border=nothing, aspect_ratio=:equal)
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This code was written by Nicolo Pollini,                                %
# Department of Wind Energy,                                              %  
# Technical University of Denmark.                                        %
#                                                                         %
# Contact: nicolo@alumni.technion.ac.il                                   %
#                                                                         %
# Code repository: https://github.com/pollinico/topopt_jl                 %
#                                                                         %
# Disclaimer:                                                             %
# The author reserves all rights but does not guarantee that the code is  %
# free from errors. Furthermore, the author shall not be liable in any    %
# event caused by the use of the program.                                 %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%