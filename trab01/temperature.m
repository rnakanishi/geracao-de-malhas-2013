% initial conditions
Tair    =   -30.0;               % Temperature of air
Tin     =    21; 

% setting initial values for grid
for i=1:(nodes)
    for j=1:(nodes)
Told(i,j)     =   Tin;
Tnew(i,j)     =   Tin;
frozen(i)   =   0;                      
latent(i)   =   Qs*mass(i)*Water/dt;    
k(i)        =   0.5;                    
cp(i)       =   cw;                     
W(i)        =   Water;                  
l(i)        =   0;                      
S(i)        =   1-Water;                
     end
end


%Simulation conditions
J       =    9;              % No. of space steps   
nodes   =   J+1;             % Number of nodes along radius or theta direction
dt      =0.1;


t       =   0;               % time index on start
tmax    =   7000;            % Time simmulation is supposed to run


R       =  d/2;

dr      =  (d/2)/J;         %  space steps in r direction

 y   =        pi/2;        % (theta) for hemisphere
dy      =  (pi/2)/J;       % space steps in Theta direction

    % Top surface condition for hemisphere
      i=nodes; 
        for j=1:1:(nodes-1) 

       Qcd_ot(i,j) = ((k(i)+ k(i-1))/2)*A(i-1)*(( Told(i,j)-Told(i-1,j))/(dr)); % heat conduction out of nod 

        Qcv(i,j) = h*(Tair-Told(i,j))*A(i); % heat transfer through convectioin on surface 

        Tnew(i,j) = ((Qcv(i,j)-Qcd_ot(i,j))/(mass(i)*cp(i))/2)*dt + Told(i,j); 
          end             %end of for loop 

  % Temperature profile for inner nodes

  for i=2:1:(nodes-1)     
    for j=2:1:(nodes-1)  

Qcd_in(i,j)=   ((k(i)+ k(i+1))/2)*A(i) *((2/R)*(( Told(i+1,j)-Told(i,j))/(2*dr)) + ((Told(i+1,j)-2*Told(i,j)+Told(i-1,j))/(dr^2)) + ((cot(y)/(R^2))*((Told(i,j+1)-Told(i,j))/(2*dy))) + (1/(R^2))*(Told(i,j+1)-2*Told(i,j)+ Told(i,j-1))/(dy^2));
Qcd_out(i,j)=  ((k(i)+ k(i-1))/2)*A(i-1)*((2/R)*(( Told(i,j)-Told(i-1,j))/(2*dr)) +((Told(i+1,j)-2*Told(i,j)+Told(i-1,j))/(dr^2)) + ((cot(y)/(R^2))*((Told(i,j)-Told(i,j-1))/(2*dy))) + (1/(R^2))*(Told(i,j+1)-2*Told(i,j)+ Told(i,j-1))/(dy^2));

Tnew(i,j)     =   (Qcd_in(i,j)-Qcd_out(i,j))/(mass(i)*cp(i))*dt + Told(i,j);
    end           
end               

   %bottom of the hemisphere solid
  Tnew(:,nodes)=-30;



 Told=Tnew;
 t=t+dt;