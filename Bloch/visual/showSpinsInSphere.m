function Showdata=showSpinsInSphere(spins,Mtype,Mindex)
%%
% Mtype: 'Mxy','Mz','Mx','My'
% Mindex:{':',':',':',':'},{1,1,':',1}....
%%
x=unique(spins.rcs(:,1));
y=unique(spins.rcs(:,2));
z=unique(spins.rcs(:,3));
cs=unique(spins.rcs(:,4));

Nx=length(x);
Ny=length(y);
Nz=length(z);
Ncs=length(cs);

Mx=reshape(spins.M(:,1),Nx,Ny,Nz,Ncs);
My=reshape(spins.M(:,2),Nx,Ny,Nz,Ncs);
Mxy=reshape(spins.M(:,1)+1i*spins.M(:,2),Nx,Ny,Nz,Ncs);
Mz=reshape(spins.M(:,3),Nx,Ny,Nz,Ncs);

switch Mtype
    case 'Mx'
        Showdata=squeeze(Mx(Mindex{:}));   
    case 'My'
        Showdata=squeeze(My(Mindex{:}));
    case 'Mz'
        Showdata=squeeze(Mz(Mindex{:}));
    case 'Mxy'
        Showdata=squeeze(Mxy(Mindex{:}));      
end
%% Create a sphere background 
r = 1;                  %radius
d = [-r,r]';
t=0:0.001:(2*pi);       % circle
t=t';
figure();
hold on;
x = [zeros(size(t)), r*cos(t), r*sin(t)];
y = [r*cos(t), zeros(size(t)), r*cos(t)];
z = [r*sin(t), r*sin(t), zeros(size(t))];
x0=[zeros(size(d)), d, zeros(size(d))];
y0=[ zeros(size(d)), zeros(size(d)), d];
z0=[d, zeros(size(d)), zeros(size(d))];
plot3(x(:,1),y(:,1),z(:,1),'r--',x(:,2),y(:,2),z(:,2),'g--',x(:,3),y(:,3),z(:,3),'c--');
plot3(x0(:,1),y0(:,1),z0(:,1),'r',x0(:,2),y0(:,2),z0(:,2),'g',x0(:,3),y0(:,3),z0(:,3),'c');
backColor = [0 0 0];
set(gcf, 'color', backColor);
set(gca, 'color', backColor);
view(3)
axis([-1 1 -1.5 1.5 -1 1]*r);axis off;
theta=-250;
phi=10;
view(theta,phi);
%%
hold on;
% a = farrow(0, 0, 0, 0, 0, 1, 'w', 2)

a = farrow(0, 0, 0, 1, 0, 0, 'w', 2)
% a = farrow(0, 0, 0, 0, 1, 0, 'w', 2)
% a = farrow(0, 0, 0, 0, 0, 1, 'w', 2)