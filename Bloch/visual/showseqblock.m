function showseqblock(block,showvalue)

t=0;B1amp=[];B1phase=[];Gx=[];Gy=[];Gz=[];
for k=1:length(block)
    SAR(k)=sum(power(block{k}.B1amp,2)*(block{1}.t(2)-block{1}.t(1)))*1000;
    SARpostion(k)=t(end)+block{k}.t(round(length(block{k}.t)*0.45));
    t=[t,t(end)+block{k}.t];
    if (max(abs(block{k}.B1amp))==0)  || ((max(abs(block{k}.B1amp))<3) && (max(abs(block{k}.B1amp))>1.5))
        B1amp=[B1amp,block{k}.B1amp];
    else
        B1amp=[B1amp,2*block{k}.B1amp/max(block{k}.B1amp)];
    end
    B1phase=[B1phase,block{k}.B1phase];
    %%
%     if (max(abs(block{k}.Gx))==0)  || ((max(abs(block{k}.Gx))<3.5) && (max(abs(block{k}.Gx))>1.5)) 
%         Gx=[Gx,block{k}.Gx];
%     else
%         Gx=[Gx,2*block{k}.Gx/max(abs(block{k}.Gx))];
%     end
%     if (max(abs(block{k}.Gy))==0) || ((max(abs(block{k}.Gy))<3.5) && (max(abs(block{k}.Gy))>1.5)) 
%         Gy=[Gy,block{k}.Gy];
%     else
%         Gy=[Gy,2*block{k}.Gy/max(abs(block{k}.Gy))];
%     end
%     if (max(abs(block{k}.Gz))==0)  || ((max(abs(block{k}.Gz))<3.5) && (max(abs(block{k}.Gz))>1.5)) 
%         Gz=[Gz,block{k}.Gz];
%     else
%         Gz=[Gz,2*block{k}.Gz/max(abs(block{k}.Gz))];
%     end
    %%
    Gx=[Gx,block{k}.Gx];
    Gy=[Gy,block{k}.Gy];
    Gz=[Gz,block{k}.Gz];   
end

if (max(abs(Gx))>0)  
    Gx=3*Gx/max(abs(Gx(:)));
end
if (max(abs(Gy))>0)  
    Gy=3*Gy/max(abs(Gy(:)));
end
if (max(abs(Gz))>0)  
    Gz=3*Gz/max(abs(Gz(:)));
end

t=t(2:end);

N=length(t);

rectangle('Position',[0,0,0.2,8],'FaceColor',[.9 .9 .9],'EdgeColor','none', 'Tag','Gzbox');
rectangle('Position',[0,8,0.2,8],'FaceColor',[.8 .8 .8],'EdgeColor','none', 'Tag','Gybox');
rectangle('Position',[0,16,0.2,8],'FaceColor',[.9 .9 .9],'EdgeColor','none','Tag','Gxbox');
rectangle('Position',[0,24,0.2,16],'FaceColor','w','EdgeColor','none','Tag','RFbox');

line(t,4*ones(1,N), 'LineWidth',2,  'Color','k','LineStyle','--' ,'Tag','LineGz0');
line(t,12*ones(1,N), 'LineWidth',2,  'Color','k','LineStyle','--' ,'Tag','LineGy0');
line(t,20*ones(1,N), 'LineWidth',2,  'Color','k','LineStyle','--' ,'Tag','LineGx0');
line(t,28*ones(1,N), 'LineWidth',2,  'Color','k','LineStyle','--' , 'Tag','RFphase0');
line(t,32*ones(1,N), 'LineWidth',2,  'Color','k','LineStyle','--' , 'Tag','RFamp0');

line(t,Gz+4,'LineWidth',2, 'Color','b','Tag','LineGz');
line(t,12+Gy,'LineWidth',2,  'Color','b', 'Tag','LineGy');
line(t,20+Gx,'LineWidth',2, 'Color','b', 'Tag','LineGx');
line(t,28+B1phase,'LineWidth',2,  'Color','b','Tag','LineB1phase');
line(t,32+B1amp*4,'LineWidth',2, 'Color','b', 'Tag','LineB1amp');

if nargin>1 && showvalue==1
    for k=1:length(block)
        text(SARpostion(k),35,num2str(SAR(k)),'FontSize',8);
        text(SARpostion(k),21,[sprintf('%.1f',max(abs(block{k}.Gx))),' G/cm'],'FontSize',8);
        text(SARpostion(k),13,[sprintf('%.1f',max(abs(block{k}.Gy))),' G/cm'],'FontSize',8);
        text(SARpostion(k),5,[sprintf('%.1f',max(abs(block{k}.Gz))),' G/cm'],'FontSize',8);
    end
end

axis([0,t(end),0,40]);
ax=gca;
set(ax,'XTick',linspace(0,t(end),3),'FontSize',12,'FontWeight','bold');
set(ax,'xticklabels',{num2str(0),[num2str(t(end)/2*1e3),' ms'], [num2str(t(end)*1e3),' ms']});
set(ax,'YTick',[4 12 20 28 32],'FontSize',12,'FontWeight','bold');
set(ax,'yticklabels',{'G_z','G_y','G_x','RFphase','RFamp'});
ax.XAxis.LineWidth=0.001;
ax.YAxis.LineWidth=0.001;
set(gcf, 'color', 'w');

