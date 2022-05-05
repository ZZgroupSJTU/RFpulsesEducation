function params=showMagnetizationWithSeqblock_Dynamic_B1RotateFrame(Mt,seqblock,params)
%%
%   Mt: 3D matrix with a size of number_of_spins * 3 (Mx,My,Mz) * Nt (Dynamic)
%   seqblock: sequence seqblock generated by function Generate_seqseqblock.m or combineseqblocks.m
%   params:Animation setting

%   Code made available for the ISMRM 2022 Educational Course
% 
%   Zhiyong Zhang (zhiyong.zhang@sjtu.edu.cn)
%% Default Animation setting and updated with input params
paramsDefault.SaveFormat=[];
paramsDefault.viewSeqblock=true;
paramsDefault.viewTrace=true;
paramsDefault.viewBfield=0;   %0 no Bfield, 1 Beff, 2 B0, B1, 3 B0 B1 Beff 4 Beff, M
paramsDefault.LastKeep=2;     %[Keep the last moment for serveral seconds]
paramsDefault.filepath=pwd;
paramsDefault.filename='SeqM_Dynamic';
paramsDefault.FrameRate=30;
paramsDefault.TracePath= true ;
paramsDefault.view3D=[-255,10];    %Default 3D view
paramsDefault.BGColor='k';
paramsDefault.B0t=zeros(1,length(seqblock.t));
paramsDefault.Ntjump=1;

if nargin>2
    params = UpdateParams(paramsDefault, params, true);
else
    params=paramsDefault;
end

graycolor=gray(64);
switch params.BGColor
    case 'k'
        params.Seqblock_baseLineColor=[0.5 0.5 0.5];
        params.Seqblock_LineColor='w';
        params.Seqblock_BGColor=graycolor(4:7,:);
    case 'w'
        params.Seqblock_baseLineColor=[0.5 0.5 0.5];
        params.Seqblock_LineColor='k';
        params.Seqblock_BGColor=graycolor(57:60,:);
    otherwise
        warning('please make sure set the Seqblock_baseLineColor,Seqblock_LineColor,Seqblock_BGColor')
end
%% show the sequence seqblock plot as a fixed background for dynamic plot
set(gcf, 'color', params.BGColor);
if params.viewSeqblock
    t=seqblock.t;
    B1amp=3*seqblock.B1amp/max(abs(seqblock.B1amp(:)+eps));
    B1phase=3*seqblock.B1phase/max(abs(seqblock.B1phase(:)+eps));
    Gx=3*seqblock.Gx/max(abs(seqblock.Gx(:)+eps));
    Gy=3*seqblock.Gy/max(abs(seqblock.Gy(:)+eps));
    Gz=3*seqblock.Gz/max(abs(seqblock.Gz(:)+eps));

    axes1=subplot(121);
    cmap=params.Seqblock_BGColor;
    rectangle(axes1,'Position',[0,0,20,8],'FaceColor',cmap(1,:),'EdgeColor','none', 'Tag','Gzbox');
    rectangle(axes1,'Position',[0,8,20,8],'FaceColor',cmap(2,:),'EdgeColor','none', 'Tag','Gybox');
    rectangle(axes1,'Position',[0,16,20,8],'FaceColor',cmap(3,:),'EdgeColor','none','Tag','Gxbox');
    rectangle(axes1,'Position',[0,24,20,16],'FaceColor',cmap(4,:),'EdgeColor','none','Tag','RFbox');

    N=length(t);
    line(axes1,t,4*ones(1,N), 'LineWidth',2,  'Color',params.Seqblock_baseLineColor,'LineStyle',':' ,'Tag','LineGz0');
    line(axes1,t,12*ones(1,N), 'LineWidth',2,  'Color',params.Seqblock_baseLineColor,'LineStyle',':' ,'Tag','LineGy0');
    line(axes1,t,20*ones(1,N), 'LineWidth',2,  'Color',params.Seqblock_baseLineColor,'LineStyle',':' ,'Tag','LineGx0');
    line(axes1,t,28*ones(1,N), 'LineWidth',2,  'Color',params.Seqblock_baseLineColor,'LineStyle',':' , 'Tag','RFphase0');
    line(axes1,t,32*ones(1,N), 'LineWidth',2,  'Color',params.Seqblock_baseLineColor,'LineStyle',':' , 'Tag','RFamp0');

    line(axes1,t,Gz+4,'LineWidth',2, 'Color',params.Seqblock_LineColor,'Tag','LineGz');
    line(axes1,t,12+Gy,'LineWidth',2,  'Color',params.Seqblock_LineColor, 'Tag','LineGy');
    line(axes1,t,20+Gx,'LineWidth',2, 'Color',params.Seqblock_LineColor, 'Tag','LineGx');
    line(axes1,t,28+B1phase,'LineWidth',2,  'Color',params.Seqblock_LineColor,'Tag','LineB1phase');
    line(axes1,t,32+2*B1amp,'LineWidth',2, 'Color',params.Seqblock_LineColor, 'Tag','LineB1amp');
    axis([0,t(end),0,40]);
    set(axes1,'XTick',linspace(0,t(end),3),'FontSize',12,'FontWeight','bold','XColor',params.Seqblock_LineColor);
    set(axes1,'xticklabels',{num2str(0),[num2str(t(end)/2*1e3),' ms'], [num2str(t(end)*1e3),' ms']});
    set(axes1,'YTick',[4 12 20 28 32],'FontSize',12,'FontWeight','bold','YColor',params.Seqblock_LineColor);
    set(axes1,'yticklabels',{'G_z','G_y','G_x','RFphase','RFamp'});
    axes1.XAxis.LineWidth=0.001;
    axes1.YAxis.LineWidth=0.001;
end
%% Create a sphere background and show as a fixed background for dynamic plot
if params.viewSeqblock
    axes2=subplot(122);
    set(axes2,'position',[0.5 0.05,0.5,0.9]);
else
    axes2=gca;
    set(axes2,'position',[0.05 0.05,0.9,0.9]);
end
set(axes2, 'color', params.BGColor);

hold on;
r = 1;                  %radius
d = [-r,r]';
t=0:0.001:(2*pi);       % circle
t=t';
x = [zeros(size(t)), r*cos(t), r*sin(t)];  
y = [r*cos(t), zeros(size(t)), r*cos(t)]; 
z = [r*sin(t), r*sin(t), zeros(size(t))];  
x0=[zeros(size(d)), d, zeros(size(d))];    
y0=[ zeros(size(d)), zeros(size(d)), d];
z0=[d, zeros(size(d)), zeros(size(d))];
plot3(x(:,1),y(:,1),z(:,1),':',x(:,2),y(:,2),z(:,2),':',x(:,3),y(:,3),z(:,3),':','color',[0.5 0.5 0.5],'linewidth',2);
plot3(x0(:,1),y0(:,1),z0(:,1),'b',x0(:,2),y0(:,2),z0(:,2),'g',x0(:,3),y0(:,3),z0(:,3),'c','linewidth',2);

view(3)
axis([-1 1 -1 1 -1 1]*r);
axis off;
view(params.view3D);

%%
Nt=size(Mt,3);
Nspins=size(Mt,1);

phase4freq=phase(exp(1i*seqblock.B1phase));
w1freqHz=(phase4freq(2:end)-phase4freq(1:end-1))./(seqblock.t(2:end)-seqblock.t(1:end-1))/2/pi;
w1freqHz=[w1freqHz,w1freqHz(end)];
if params.viewBfield==3
   B1field=cat(3,repmat(seqblock.B1amp,[Nspins,1]),zeros(Nspins,Nt),zeros(Nspins,Nt));  
else
   B1field=cat(3,repmat(seqblock.B1amp.*cos(seqblock.B1phase),[Nspins,1]),repmat(seqblock.B1amp.*sin(seqblock.B1phase),[Nspins,1]),zeros(Nspins,Nt));
end
B0field=cat(3,zeros(Nspins,Nt),zeros(Nspins,Nt),1e-3*(params.B0t-repmat(w1freqHz,[Nspins,1])));
Beff=cat(3,repmat(seqblock.B1amp,[Nspins,1]),zeros(Nspins,Nt),1e-3*(params.B0t-repmat(w1freqHz,[Nspins,1])));
Beffnorm=sqrt(Beff(:,:,1).^2+Beff(:,:,2).^2+Beff(:,:,3).^2);
Beffphase=atan(Beff(:,:,1)./abs(Beff(:,:,3)));
BeffphaseRate=(Beffphase(:,2:end)-Beffphase(:,1:end-1))./(seqblock.t(2:end)-seqblock.t(1:end-1))/2/pi;
BeffphaseRate=abs(cat(2,BeffphaseRate,BeffphaseRate(:,end))./1e3);
Beff=Beff./max(Beffnorm(:));
% Beff=Beff./repmat(Beffnorm,[1 1 3]);
B1field=B1field./max(Beffnorm(:));
B0field=B0field./max(Beffnorm(:));
%%
hold on;
cspins = spring(Nspins);   % create a colormap for many spins' dynamics
if params.viewTrace
    h_M=zeros(1,Nspins);         % For the animatedline handle indicators of the tracing of the spin pathways.
    h_B1eff=zeros(1,Nspins);         % For the animatedline handle indicators of the tracing of the spin pathways.
    h_B1=zeros(1,Nspins);         % For the animatedline handle indicators of the tracing of the spin pathways.
    h_B0=zeros(1,Nspins);         % For the animatedline handle indicators of the tracing of the spin pathways.

    
    linetrace=cell(1,Nspins);  % Init the animatedline handles for all the spins
    linetrace_B1eff=cell(1,Nspins);  % Init the animatedline handles for all the spins
    linetrace_B0=cell(1,Nspins);  % Init the animatedline handles for all the spins
    linetrace_B1=cell(1,Nspins);  % Init the animatedline handles for all the spins      
    for spin_idx=1:Nspins
        switch params.viewBfield
            case 0
                 linetrace{spin_idx} = animatedline(axes2,'linewidth', 0.8,'color',cspins(spin_idx,:),'linestyle',':');
            case 1
                 linetrace_B1eff{spin_idx} = animatedline(axes2,'linewidth', 0.8,'color',params.Seqblock_LineColor,'linestyle',':');           % Init Animated Line
            case 2
                 linetrace_B0{spin_idx} = animatedline(axes2,'linewidth', 0.8,'color','y','linestyle',':');           % Init Animated Line
                 linetrace_B1{spin_idx} = animatedline(axes2,'linewidth', 0.8,'color','r','linestyle',':');           % Init Animated Line
            case 3
                 linetrace_B1eff{spin_idx} = animatedline(axes2,'linewidth', 0.8,'color',params.Seqblock_LineColor,'linestyle',':');           % Init Animated Line
                 linetrace_B0{spin_idx} = animatedline(axes2,'linewidth', 0.8,'color','y','linestyle',':');           % Init Animated Line
                 linetrace_B1{spin_idx} = animatedline(axes2,'linewidth', 0.8,'color','r','linestyle',':');           % Init Animated Line
            case 4
                 linetrace_B1eff{spin_idx} = animatedline(axes2,'linewidth', 0.8,'color',params.Seqblock_LineColor,'linestyle',':');           % Init Animated Line
                 linetrace{spin_idx} = animatedline(axes2,'linewidth', 0.8,'color',cspins(spin_idx,:),'linestyle',':');
        end
    end
end

if strcmp(params.SaveFormat,'mp4')
    mp4name=fullfile(params.filepath,[params.filename,'.mp4']);
    writerObj=VideoWriter(mp4name,'MPEG-4');
    writerObj.FrameRate=params.FrameRate;
    open(writerObj);
end
            
for t_idx=1:params.Ntjump:Nt             % Start the dynamic plot loop
    if t_idx>1
       switch params.viewBfield
           case 0
               delete(h_M);          % remove the previous Animated Lines in the spin dynamic plot
           case 1
               delete(h_B1eff);          % remove the previous Animated Lines in the spin dynamic plot
           case 2
               delete(h_B1);          % remove the previous Animated Lines in the spin dynamic plot
               delete(h_B0);          % remove the previous Animated Lines in the spin dynamic plot
           case 3
               delete(h_B1eff);          % remove the previous Animated Lines in the spin dynamic plot
               delete(h_B1);          % remove the previous Animated Lines in the spin dynamic plot
               delete(h_B0);          % remove the previous Animated Lines in the spin dynamic plot      
           case 4
               delete(h_B1eff);          % remove the previous Animated Lines in the spin dynamic plot
               delete(h_M);          % remove the previous Animated Lines in the spin dynamic plot
       end
       if params.viewSeqblock
          delete(htline);     % remove the previous moment line in the sequence seqblock plot
       end
    end
    
    if params.viewSeqblock
        set(gcf,'CurrentAxes',axes1);         
        hold on;               % moment line plot (t line) in the sequence seqblock plot
        htline=line([seqblock.t(t_idx),seqblock.t(t_idx)],[0 42], 'LineWidth',1,  'Color','r' ); 
    end
    set(gcf,'CurrentAxes',axes2);
    hold on;
    for spin_idx=1:Nspins
        x_cur=Mt(spin_idx,1,t_idx);
        y_cur=Mt(spin_idx,2,t_idx);
        z_cur=Mt(spin_idx,3,t_idx);
        phi=seqblock.B1phase(t_idx);
        x_cur_B1frame=x_cur*cos(phi)-y_cur*sin(phi);
        y_cur_B1frame=x_cur*sin(phi)+y_cur*cos(phi);
        if params.viewTrace
           switch params.viewBfield
                case 0
                    addpoints(linetrace{spin_idx},x_cur_B1frame,y_cur_B1frame,z_cur);                  % - add to animated line
                case 1
                    addpoints(linetrace_B1eff{spin_idx},Beff(spin_idx,t_idx,1),Beff(spin_idx,t_idx,2),Beff(spin_idx,t_idx,3));                  % - add to animated line
                case 2
                    addpoints(linetrace_B1{spin_idx},B1field(spin_idx,t_idx,1),B1field(spin_idx,t_idx,2),B1field(spin_idx,t_idx,3));                  % - add to animated line
                    addpoints(linetrace_B0{spin_idx},B0field(spin_idx,t_idx,1),B0field(spin_idx,t_idx,2),B0field(spin_idx,t_idx,3));                  % - add to animated line
                case 3
                    addpoints(linetrace_B1eff{spin_idx},Beff(spin_idx,t_idx,1),Beff(spin_idx,t_idx,2),Beff(spin_idx,t_idx,3));                  % - add to animated line
                    addpoints(linetrace_B1{spin_idx},B1field(spin_idx,t_idx,1),B1field(spin_idx,t_idx,2),B1field(spin_idx,t_idx,3));                  % - add to animated line
                    addpoints(linetrace_B0{spin_idx},B0field(spin_idx,t_idx,1),B0field(spin_idx,t_idx,2),B0field(spin_idx,t_idx,3));                  % - add to animated line
               case 4
                    addpoints(linetrace{spin_idx},x_cur_B1frame,y_cur_B1frame,z_cur);                  % - add to animated line
                    addpoints(linetrace_B1eff{spin_idx},Beff(spin_idx,t_idx,1),Beff(spin_idx,t_idx,2),Beff(spin_idx,t_idx,3));                  % - add to animated line
           end
        end

        switch params.viewBfield
            case 0
                h_M(spin_idx)=farrow(0,0,0,x_cur_B1frame,y_cur_B1frame,z_cur,cspins(spin_idx,:),2);        
            case 1
                h_B1eff(spin_idx)=farrow(0,0,0,Beff(spin_idx,t_idx,1),Beff(spin_idx,t_idx,2),Beff(spin_idx,t_idx,3),params.Seqblock_LineColor,2);
            case 2
                h_B1(spin_idx)=farrow(0,0,0,B1field(spin_idx,t_idx,1),B1field(spin_idx,t_idx,2),B1field(spin_idx,t_idx,3),'r',2);
                h_B0(spin_idx)=farrow(0,0,0,B0field(spin_idx,t_idx,1),B0field(spin_idx,t_idx,2),B0field(spin_idx,t_idx,3),'y',2);
            case 3
                h_B1eff(spin_idx)=farrow(0,0,0,Beff(spin_idx,t_idx,1),Beff(spin_idx,t_idx,2),Beff(spin_idx,t_idx,3),params.Seqblock_LineColor,2);
                h_B1(spin_idx)=farrow(0,0,0,B1field(spin_idx,t_idx,1),B1field(spin_idx,t_idx,2),B1field(spin_idx,t_idx,3),'r',2);
                h_B0(spin_idx)=farrow(0,0,0,B0field(spin_idx,t_idx,1),B0field(spin_idx,t_idx,2),B0field(spin_idx,t_idx,3),'y',2);
            case 4
                h_M(spin_idx)=farrow(0,0,0,x_cur_B1frame,y_cur_B1frame,z_cur,cspins(spin_idx,:),2);        
                h_B1eff(spin_idx)=farrow(0,0,0,Beff(spin_idx,t_idx,1),Beff(spin_idx,t_idx,2),Beff(spin_idx,t_idx,3),params.Seqblock_LineColor,2);
       	end
     end
    
    switch params.SaveFormat
        case 'none'
            pause(0.01);
        case 'gif'
            gifname=fullfile(params.filepath,[params.filename,'.gif']);
            frame=getframe(gcf);
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 
            % Write to the GIF File 
            if t_idx == 1 
               imwrite(imind,cm,gifname,'gif','DelayTime',1/params.FrameRate, 'Loopcount',inf); 
            elseif t_idx==Nt
               for m=1:max(1,round(params.LastKeep/params.FrameRate))
                   imwrite(imind,cm,gifname,'gif','DelayTime',1/params.FrameRate,'WriteMode','append');
               end
            else
               imwrite(imind,cm,gifname,'gif','DelayTime',1/params.FrameRate,'WriteMode','append'); 
            end 
        case 'mp4'
            if t_idx==Nt
                for m=1:max(1,round(params.LastKeep/params.FrameRate))
                    frame=getframe(gcf);
                    writeVideo(writerObj,frame);
                end
            else
               frame=getframe(gcf);
               writeVideo(writerObj,frame);
            end
    end
end
if strcmp(params.SaveFormat,'mp4')
    close(writerObj);
end
