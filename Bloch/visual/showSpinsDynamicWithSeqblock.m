function showSpinsDynamicWithSeqblock(Mt,block,gifname)
%%
% Mt: nspin*3*Nt

%% Create a sphere background 
r = 1;                  %radius
d = [-r,r]';
t=0:0.001:(2*pi);       % circle
t=t';
axes1=subplot(121);showseqblock({block})
axes2=subplot(122);hold on;
x = [zeros(size(t)), r*cos(t), r*sin(t)];
y = [r*cos(t), zeros(size(t)), r*cos(t)];
z = [r*sin(t), r*sin(t), zeros(size(t))];
x0=[zeros(size(d)), d, zeros(size(d))];
y0=[ zeros(size(d)), zeros(size(d)), d];
z0=[d, zeros(size(d)), zeros(size(d))];
plot3(x(:,1),y(:,1),z(:,1),'b:',x(:,2),y(:,2),z(:,2),'g:',x(:,3),y(:,3),z(:,3),'c:','linewidth',2);
plot3(x0(:,1),y0(:,1),z0(:,1),'b',x0(:,2),y0(:,2),z0(:,2),'g',x0(:,3),y0(:,3),z0(:,3),'c','linewidth',2);
set(gcf, 'color', 'w');
set(gca, 'color', 'w');
view(3)
axis([-1 1 -1 1 -1 1]*r);axis off;
theta=-255;
phi=10;
view(theta,phi);
% set(axes2,'position',[0.5 0.05,0.5,0.9]);
%%
Nt=size(Mt,3);
Nspins=size(Mt,1);
hold on;
linetrace=cell(1,Nspins);
for spin_idx=1:Nspins
    linetrace{spin_idx} = animatedline(axes2,'linewidth', 1,'color',[.6 .2 .5]);           % Init Animated Line
end

if nargin>2 && ~ischar('gifname')
    gifname='SpinDynamic.gif';
end
h=zeros(1,Nspins);
cspins = spring(Nspins);
for t_idx=1:Nt
    if t_idx>1
       delete(h);
       delete(htline);
    end
    
    subplot(121);
    hold on;htline=line([block.t(t_idx),block.t(t_idx)],[0 40], 'LineWidth',2,  'Color','r','LineStyle','--' ); 
    subplot(122);
    hold on;
    for spin_idx=1:Nspins
        x_cur=Mt(spin_idx,1,t_idx);
        y_cur=Mt(spin_idx,2,t_idx);
        z_cur=Mt(spin_idx,3,t_idx);
%         addpoints(linetrace{spin_idx},x_cur,y_cur,z_cur);                  % - add to animated line
        h(spin_idx)=farrow(0,0,0,x_cur,y_cur,z_cur,cspins(spin_idx,:),2);
    end
    if (nargin>2)
        frame=getframe(gcf);
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
          % Write to the GIF File 
        if t_idx == 1 
           imwrite(imind,cm,gifname,'gif','DelayTime',0.1, 'Loopcount',inf); 
        else
           imwrite(imind,cm,gifname,'gif','DelayTime',0.1,'WriteMode','append'); 
        end 
    else
        pause(0.01);
    end
end

