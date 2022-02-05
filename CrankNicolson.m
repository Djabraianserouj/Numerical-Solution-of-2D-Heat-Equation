%initializing dimensions for discretization in space
n = 41;                 %number of nodes
L = 1;                  %length
m =(n-2)*(n-2);         %dimension of the vector containing inner node temperature values
sm = sqrt(m);           %dimension of matrix containing inner node temperature values
h = L/(n-1);            %height, which is in this case 1/40

%initializing temperatuer vector with all the inner node temperature vectors
% and the corresponding temperature solution matrix
T_array  = zeros(m,1);    
T_solution = zeros(n,n);
T_contour = zeros(5,m);  

%% Setting up the boundary conditions into arrays
T_solution(1,1:n) = 0;         %Top Wall
T_solution(n,1:n) = 1;         %Bottom Wall      
for i_Tl = 1 : 1 : n           %Left wall, applying the boundary condition
  y = ( n- i_Tl )/( n - 1 );
  T_solution(i_Tl, 1 ) = 1 - y^3;
end
for i_Tr = 1 : 1 : n           %Right wall, applying the boundary condition
  y = ( n - i_Tr )/( n - 1 );
  T_solution(i_Tr, n ) = 1 - sin(y*pi/2); 
end
T_contour_matrix = permute(repmat(T_solution, 1,1,5), [3 1 2]);

%initializing coefficient matrices ( A and C) and boundary conditions vector C
A  = zeros(m,m);  
C = A;     
B  = zeros(sm,sm);

%initiliazing parameters for time discretization
dt = 0.0001;                
t = 0 : dt : 0.16;
r = dt /(h^2);  

% vectors for figures and assignment questions
T_time_evolution = zeros(length(t),1);
T_vertical_profile = zeros(sm+2, 1);

%% Setup coefficient matrix A, coefficient matrix for temperature at n+1
for ix = 1 : 1 : m
    for jx =1 : 1 : m
        if (ix == jx)
                A(ix,jx) = (1 + 2*r); 
        elseif ( (ix == jx + 1) && ( (ix - 1) ~= sm * round( (ix-1)/sm) ) ) %RHS
                A(ix,jx) = -r/2;
        elseif ( (ix == jx - 1) && ( ix ~= sm*round(ix/sm) ) )
                A(ix,jx) = -r/2;
        elseif (ix == jx + sm)
                A(ix,jx) = -r/2;
        elseif (jx == ix + sm)
                A(ix,jx) = -r/2;
        else
                A(ix,jx) = 0;
        end
    end
end
 
%Setup boundary condition array
for iy = 1 : 1 : sm
    for jy = 1 : 1 : sm
        if (iy == 1) && (jy == 1)
            B(iy,jy) = r *(T_solution(1,2)+ T_solution(2,1)); %Tl + Tt
        elseif (iy == 1) && (jy == sm)
            B(iy,jy) = r *( T_solution(1,sm +1) + T_solution(2,sm+2));   %Tt + Tr                         %RHS
        elseif (iy == sm) && (jy == sm)
            B(iy,jy) = r * (T_solution(sm+2,sm +1) + T_solution(sm+1,sm+2)); % Tb + Tr
        elseif (iy == sm) && (jy == 1)
            B(iy,jy) = r * (T_solution(sm+1,1) + T_solution(sm+2,2)); % Tb + Tl
        elseif (iy == 1)&&(jy > 1 || jy < sm)
            B(iy,jy) = 0;
        elseif (jy == sm) && (iy > 1 || iy < sm)
            B(iy,jy) = r * T_solution(iy + 1, sm + 2);
        elseif (iy == sm) && ( jy > 1 || jy < sm)
            B(iy,jy) = r * 1;
        elseif (jy == 1) && ( iy > 1 || jy < sm)
            B(iy,jy) = r * T_solution(iy+1, 1);
        else
            B(iy,jy) = 0;
        end
    end   
end
Bx = reshape(B,[],1);  %Convert matrix to array
 
% Setup coefficient matrix C for temperature array at n
for iz = 1 : 1 : m
    for jz =1 : 1 : m
        if (iz == jz)
                C(iz,jz) = (1 - 2*r); 
        elseif ( (iz == jz + 1) && ( (iz - 1) ~= sm * round( (iz-1)/sm) ) ) %RHS
                C(iz,jz) = r/2;
        elseif ( (iz == jz - 1) && ( iz ~= sm*round(iz/sm) ) )
                C(iz,jz) = r/2;
        elseif (iz == jz + sm)
                C(iz,jz) = r/2;
        elseif (jz == iz + sm)
                C(iz,jz) = r/2;
        else
                C(iz,jz) = 0;
        end
    end
end
%% Solution
selected_time = [ 0.01, 0.02, 0.04, 0.08, 0.16 ];
counter = 1;
for l = 1 : length(t)   %time steps
     Xx = ( C*T_array + Bx ); %RHS of the equation where all the values are known
     T_array = A \ Xx; % LHS of the equation to get value of T at n+1
     T_time_evolution(l) = T_array(15*39+24);
     
     %storing values for contour plots at different t
     if any( selected_time == t(l))
       T_contour(counter,:) = T_array;
       T_contour_matrix(counter, 2:n-1, 2:n-1 ) = reshape( T_array , sm , sm);
       counter = counter + 1;
     endif
       
     
    
end
%Convert inner node temperature solution to matrix
T_matrix = reshape( T_array , sm , sm); %convert array to matrix


%Appending the inner node temperature solution matrix component to the whole temperature solution matrix

T_solution(2:n-1, 2:n-1 ) = T_matrix;
T_vertical_profile = flip(T_solution( :, 17 ));



%initializing meshgrid with 41 evenly spaced points in length and width
x = linspace(0,L,n);    
y = x;    
[X,Y]=meshgrid(x,y);   
%% Plot
figure( 'Position', [0, 0, 1000, 1200]);
for i = 1 : 1 :5
  subplot(3,2,i)
  contourf(X,flipud(Y),squeeze(T_contour_matrix(i,:,:)),10);
  ylabel('Y(dimensionless)','FontSize',10);
  xlabel('X(dimensionless)','FontSize',10);
  title(sprintf('Crank Nicolson scheme t = %.2f s',  selected_time(i)), 'FontSize',10);
  colorbar
  
  
end
print -dpng -color "-S1000,1200" ContourSubplots.png

figure
plot(t,T_time_evolution)
xlabel('Time(s)')
ylabel('Temperature(dimensionless)')
title('Temperature variation at x=y=0.4 against time')


figure
plot(y,T_vertical_profile)
xlabel('Y(dimensionless)')
ylabel('Temperature(dimensionless)')
title('Vertical T profile for x=0.4 at t=16s')
