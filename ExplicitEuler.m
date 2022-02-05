%%% Euler Explicit in time with CDS scheme in space %%%
clear all 
clc

%%% Setting up variables & conditions %%%
T=0.16;           %Total Time 
dt=0.0001;        %Stepsize
t=T/dt;           %Total Number of steps 
dx=1/40;          %Node size in x-direction 
dy=1/40;          %Node size in y-direction
o=0:dx:1;         %used for plotting(Plate dimension)
p=0:dy:1;         %used for plotting(Plate dimension)
Lx=1;             %Plate size in x-direction
Ly=1;             %Plate size in y-direction
Nx=(Lx/dx)+1;     %Nb of nodes in x-direction
Ny=(Ly/dy)+1;     %Nb of nodes in y-direction  
w=zeros(Nx,Ny,t); %initializing a matrix of zeroes
y=0:dx:1;
for a=1:1:Ny
  w(a,1,:)=1-y(a)^3;             %BC Left Side
endfor
for b=1:1:Nx
  w(b,Nx,:)=1-sin((pi/2)*y(b));  %BC Right Side
endfor
w(1,:,:)=1;                      %BC Top Side
w(Nx,:,:)=0;                     %BC Bottom Side

%%% Solution for a 3D matrix %%%
for k=1:1:t       
  for i=2:1:Nx-1
    for j=2:1:Ny-1
      w(i,j,k+1)=w(i,j,k)+dt*(((w(i-1,j,k)-2*w(i,j,k)+w(i+1,j,k))/dx^2)+((w(i,j-1,k)-2*w(i,j,k)+w(i,j+1,k))/dy^2));
    endfor
  endfor
  k
endfor 

%%%Storing final temperature values%%%
vecsize=(Nx-2)*(Ny-2);                 %Total number of elements
Tr=transpose(w(2:1:Nx-1,2:1:Ny-1,t));  
M=reshape(Tr,vecsize,1);               %From matrix to a vector
%These 3 above are use to export Temperature values 
%into an excel sheet for example by using a writing function

%%% Plotting the data from the solution matrix %%%
figure 1
contourf(o,p,w(:,:,t),10);
colorbar;
ylabel('Y(dimensionless)')
xlabel('X(dimensionless)')
title('Temperature profile at t=0.16s')
saveas(1,'ContourPlot.png')

%%% Plot for stable dt at x=y=0.4 for T=0.16s %%%
for i=1:1:t
  v(i)=w(0.4/dx+1,0.4/dx+1,i);  %Temperature value
endfor 
timevector=[0:dt:T-dt];
figure 2
plot (timevector,v)
xlabel('Time(s)')
ylabel('Temperature(dimensionless)')
title('Temperature variation at x=y=0.4 for t=0.16s')
saveas(2,'stable x=y=0.4.png')

%%% Plot for stable dt at x=0.4 at t=0.16 %%%
for i=1:1:Ny
  B(i)=w(i,0.4/dx+1,t);        %Temperature value
endfor  

figure 3
plot (o,B)
xlabel('Y(dimensionless)')
ylabel('Temperature(dimensionless)')
title('Vertical T profile for x=0.4 at t=16s')
saveas(3,'stable x=4.png')

 