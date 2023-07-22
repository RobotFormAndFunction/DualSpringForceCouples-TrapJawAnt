%%% dual spring force couple mandible simulation
%%% R. St. Pierre - University at Buffalo

function [poweroutput, poweroutput2] =...
    mandiblelever_energypartion_fxn(Lratio,partition,ktval)

% mandible model, direct newton approach, using assumptions from Greg
Etot=24.2e-6; %J mean is 24.2 and std is 6.0

m=54.7e-9; %kg mean is 54.7 and std is 10.0
L=1.4e-3; %m mean is 1.37 and std is 0.06
I=1/12 * m * L^2; %kg-m^2 - mass moment of inertia at the geometric 
%center of the mandible


a=partition; %energy split E_tensile/Etotal
kt=ktval; %N/m tensile spring

%some algebra relating the two springs given the energy partitions
%Et/Etot * Etot/Ec = Et/Ec = kc/kt   this comes from springs in series
%partition formulas
%(a)/(1-a)=kc/kt -> kc=kt*(a)/(1-a);

kc=kt*(a/(1-a)); %N/m setting kc based upon the energy partition

%some algebra relating energy and displacements
%Et/Etot * Etot/Ec = Et/Ec = xt/xc   this comes from springs in series
%partition formulas
%Etot=1/2 kt xt^2 + 1/2 kc xc^2
%2*Etot = kt xt^2 + kc xc^2
%2*Etot/xc^2 = kt xt^2/xc^2 + kc
%2*Etot/xc^2 = kt (a/(1-a))^2 + kc

xc=sqrt((2*Etot)/(kt*(a/(1-a))^2+kc)); %solving the initial displacement
%of the compression spring given the energy partion factor, total energy,
%and spring constant
xt=(a/(1-a))*xc; %related through the partition formula

%lever arm lengths
Lt=L; %m - distance from mandible tip to tensile spring attachment
Lc=Lratio*L; %m - distance from mandible tip to compression spring 
%attachment

%find point between two springs, consider that the moment of inertia 
d=(Lc-L/2)+(Lt-Lc)/2; %m half the distance between the two springs plus 
% half the mandible length. This is for the mass moment of inertia at the
% spring midpoint
Ip=I+m*d^2; %kg-m^2 - mass moment of inertia at the 



% simulation setup
dt=1e-9; %s - time step
t_end=1000e-6; %how long the simulation should take
t_vec=[0:dt:t_end]'; %initial a time vector

int=[0;0;0;0]; % initial conditions, these are the intiail conditions of
%theta, y, and their derivatives
dy(1,1)=int(1);
y(1,1)=int(2);
dtheta(1,1)=int(3);
theta(1,1)=int(4);



%this is the forward newton-euler method for simulating the mandible
%closure
for i=1:length(t_vec)
    
    Ft(i,1)=kt*xt+kt*y(i,1)+kt*(Lt-Lc)*0.5*sin(theta(i,1));
    if Ft(i,1) < 0 %the tendon spring can only pull, not push so this
        Ft(i,1)=0; %if stattement catches that
        %'slack'
    end
    Fc(i,1)=kc*xc-kc*y(i,1)+kc*(Lt-Lc)*0.5*sin(theta(i,1));
    
    %to pin the mandible at the point between the two springs, set all of
    %these values to zero
    ddy(i,1)=1/m*(Fc(i,1)-Ft(i,1));
    dy(i+1,1)=ddy(i,1)*dt+dy(i,1);
    y(i+1,1)=dy(i,1)*dt+y(i,1);
    
    ddtheta(i,1)=1/Ip*(-Ft(i,1)*(Lt-Lc)*0.5*cos(theta(i,1))-...
        Fc(i,1)*(Lt-Lc)*0.5*cos(theta(i,1)));
    dtheta(i+1,1)=ddtheta(i,1)*dt+dtheta(i,1);
    theta(i+1,1)=dtheta(i,1)*dt+theta(i,1);
    
    r(i,1)=-ddy(i,1)/ddtheta(i,1);
    
    %calculating various instantaneous powers and energies
    Power(i,1)=(abs(Ip*ddtheta(i,1)*dtheta(i,1))+abs(m*ddy(i,1)*dy(i,1)));
    Powerrot(i,1)=abs(Ip*ddtheta(i,1)*dtheta(i,1));
    Powerlin(i,1)=abs(m*ddy(i,1)*dy(i,1));
    KE(i,1)=0.5*Ip*dtheta(i,1)^2+0.5*m*dy(i,1)^2;
    KErot(i,1)=0.5*Ip*dtheta(i,1)^2;
    KElin(i,1)=0.5*m*dy(i,1)^2;
    
    %set break at 65 deg instead
    %65/360 * 2 * pi
    
    if theta(i,1) <= -65/360*2*pi || theta(i,1) >= 65/360*2*pi
        additionaltime=abs((25/360*2*pi)/dtheta(i,1));
        %the spring actuation only powers the first 65 degrees of rotation,
        %this additional time captures how much additional time is needed
        %to rotate to 90 degrees given the rotational velocity when the
        %mandibles are at 65 degrees
       
        break
    else
        additionaltime=0;
    end 
    
end
%%
    if t_vec(i)==t_end 
        KErot=nan(i,1);        
    end    



poweroutput=KErot(end)/(t_vec(length(KErot))+additionaltime);
poweroutput2=KErot(end)/(t_vec(length(KErot))+additionaltime+(1/300000));
%power output 2 has an additional 'frame' of data because it is hard to
%figure out when the mandibles start moving in real life. 

end
