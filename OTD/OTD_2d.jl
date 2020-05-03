#!/usr/bin/julia
# Source code for Exercise 2 of Homework 1


println("Visualizing the changing approximated eigenvector/eigenvalue")
println("from the OTD method in a 2x2 system, for a single OTD mode")

# import Pkg
# #Pkg.add("PyPlot")
# #Pkg.add("Blink")
# 
# #using Blink
# using PyPlot,Colors,PyCall
# #using Plots
# using LinearAlgebra

lafs = 16;


# Temporal discretization
bdf1  = [ 1.  -1.  0.  0.]/1.;
bdf2  = [ 3.  -4.  1.  0.]/2.;
bdf3  = [11. -18.  9. -2.]/6.;
ex0   = [0. 1.  0. 0.];
ex1   = [0. 2. -1. 0.];
ex2   = [0. 3. -3. 1.];

close("all")

dt     = 0.01;
Nstep  = 10000;
egvupd = 50;
ifrot  = true;    # Rotate Matrix


# Create Matrix A
v0 = [-0.01 -0.02];

n = length(v0);
A = zeros(Float64,n,n);

for i in 1:n
  A[i,i] = copy(v0[i]);
end

if (ifrot)
  global A    
  U = rand(Float64,n,n);
  U = U/norm(U);

  v = U[:,1];
  U[:,1] = v/norm(v);    
  for i=2:n
    v = U[:,i];
    α = U[:,1:i-1]'*v;
    v = v - U[:,1:i-1]*α;
    v = v/norm(v);
    U[:,i] = v;
  end
  
  A  = U'*A*U;
end  


nmodes = 1;
Vinit  = rand(Float64,n,nmodes);
V      = zeros(Float64,n,nmodes);
Vlag   = zeros(Float64,n,3,nmodes);
R1lag  = zeros(Float64,n,2,nmodes);
R2lag  = zeros(Float64,n,2,nmodes);

v = copy(Vinit[:,1]);
#v[1] = 1.e-1;
#v[2] = 1.;
v = v/norm(v);
Vinit[:,1] = v;

# Orthogonalize initial field
for i=2:nmodes
  v = Vinit[:,i];
  alpha = (Vinit[:,1:i-1])'*v;
  w = v - Vinit[:,1:i-1]*alpha;
  w = w/norm(w);
  Vinit[:,i] = w;
end

V   = copy(Vinit);
Rhs = 0*copy(V);

Evals = zeros(Complex,Nstep,nmodes);
Ermax = zeros(Complex,Nstep);

t     = 0.

v1    = [0, 1.];
v2    = [1., 0.];

#h1  = figure(num=1,figsize=[18.,6.]);
#ax1 = subplot(121);
#ax2 = subplot(122);
 

time = range(dt,step=dt,length=Nstep);

F   = eigen(A);  
ee0 = F.values;
ee  = sort(ee0,rev=true);
emax1 = ee[1]*ones(Float64,Nstep);
emax2 = ee[2]*ones(Float64,Nstep);

pl11 = plot(time,emax1,lw=1.);
pl12 = plot(time,emax2,lw=1.);


xlabel(L"time",fontsize=lafs)
ylabel(L"\lambda",fontsize=lafs)
title("Approximated Eigenvalue")

#println("Plotted the first two")

# Plot the eigenvectors
egvs  = F.vectors
cm    = get_cmap("tab10");
rgba0 = cm(0); 
rgba1 = cm(1); 

#ax2.arrow(0.,0.,egvs[1,1],egvs[2,1],width=0.02,length_includes_head=true,color=rgba1);
#ax2.arrow(0.,0.,egvs[1,2],egvs[2,2],width=0.02,length_includes_head=true,color=rgba0);

ax2.set_xlabel(L"x_{1}",fontsize=lafs)
ax2.set_ylabel(L"x_{2}",fontsize=lafs)
ax2.set_xlim([-1.1,1.1])
ax2.set_ylim([-1.1,1.1])

sleep(0.001)

for i in 1:Nstep
  global V, Vlag,Rlag,Rhs,Evals,Ermax
  global t, ax2
     
  t = t + dt;

  bdf = bdf3;
  ext = ex2;
     
  if i==1
    bdf = bdf1;
    ext = ex0;
  elseif i==2  
    bdf = bdf2;
    ext = ex1;
#  elseif i==3 
#    bdf = bdf3;
#    ext = ex2;
  end

  A1     = copy(A);

  Ar = V'*A1*V;
  ee = eigvals(Ar);
  Evals[i,:] = ee;

#  Ar = Ar + Ar';
  ee = eigvals(Ar);
  Ermax[i] = maximum(ee);

  for j in 1:nmodes
    global V,Vlag,R1lag,R2lag,Rhs

#   Extrapolate A*x
    rhs     = A1*V[:,j];
    rhs1    = ext[2]*rhs + ext[3]*R1lag[:,1,j] + ext[4]*R1lag[:,2,j];

    R1lag[:,2,j] = R1lag[:,1,j];
    R1lag[:,1,j] = rhs;


    Ax      = A1*V[:,j]; 
    alpha   = V'*Ax; 
    w       = V*alpha;
    rhs2    = 1*(ext[2]*w + ext[3]*R2lag[:,1,j] + ext[4]*R2lag[:,2,j]);


    R2lag[:,2,j] = R2lag[:,1,j];
    R2lag[:,1,j] = w;

    bdlag    = 1. /dt*(bdf[2]*V[:,j] + bdf[3]*Vlag[:,1,j] + bdf[4]*Vlag[:,2,j]); 
    Rhs[:,j] = rhs1 - rhs2 - bdlag;

    Vlag[:,3,j] = Vlag[:,2,j];
    Vlag[:,2,j] = Vlag[:,1,j];
    Vlag[:,1,j] = V[:,j];

    V[:,j]      = copy(Rhs[:,j])*dt/bdf[1];

    vdiff       = V[:,j] - Vlag[:,1,j];  

    if (mod(i,egvupd)==0)
      global vdiff    
      global pl2,pl3,pl4

      if i>egvupd    
        pl2.remove();
#        pl3.remove();
        pl4[1].remove();        
      end  

#      ax2.plot([0., V[1,j]],[0., V[2,j]]);
      pl2 = ax2.arrow(0.,0.,V[1,j],V[2,j],width=0.01,color="black",length_includes_head=true);

#      scale = 5000.;
#      pl3 = ax2.arrow(Vlag[1,j],V[2,j],scale*vdiff[1],scale*vdiff[2],width=0.01,color="red",length_includes_head=true);     

#      ax2.set_xlabel(L"x_{1}")
#      ax2.set_ylabel(L"x_{2}")
#      ax2.set_xlim([-1.25,1.25])
#      ax2.set_ylim([0.,1.25])
#      ax2.set_title("Approximated Eigenvctor")
 
      
      pl4 = ax1.plot(time[1:i],real(Evals[1:i,j]),color="black");

#      IJulia.clear_output(true)

      sleep(0.001)
    end  

  end

# Orthogonalize field
# Do we also need to rescale?
  v = V[:,1];
  v = v/norm(v);
  V[:,1] = v;

  for j=2:nmodes
    v = V[:,j];
    alpha = (V[:,1:j-1])'*v;
    w = v - V[:,1:j-1]*alpha;
    w = w/norm(w);
    V[:,j] = w;
  end


end


println("Done.")





