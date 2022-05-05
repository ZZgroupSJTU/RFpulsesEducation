function Showdata=showspins(spins,Mtype,Mindex,subplotmatIn)
%%
% Mtype: 'Mxy','Mz','Mx','My'
% Mindex:{':',':',':',':'},{1,1,':',1}....
% subplotmat: [ncor nrow namp] or [ncor nrow namp nph]


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

subplotmat=[2 1 1 2];
if nargin==4
    subplotmat=subplotmatIn;
end
    
if length(subplotmat)<3
    subplotmat=[2 1 1 2];
end

switch nnz(size(Showdata)>1)
    case 1
        for k=1:4
            if Mindex{k}==':'
                xdim=unique(spins.rcs(:,k));
                ydim=reshape(Showdata,1,[]);
                if strcmp(Mtype,'Mz')                    
                    hold on; subplot(subplotmat(1),subplotmat(2),subplotmat(3));plot(xdim,abs(ydim),'.-'); title('Amplitute plot'); ylim([-1,1]);
                    if length(subplotmat)==4
                        hold on;  subplot(subplotmat(1),subplotmat(2),subplotmat(4));plot(xdim,phase(ydim),'.-');title('Phase plot');
                    end
                else 
                    hold on; subplot(subplotmat(1),subplotmat(2),subplotmat(3));plot(xdim,abs(ydim),'.-'); title('Amplitute plot'); ylim([0,1]);
                    if length(subplotmat)==4
                        hold on;  subplot(subplotmat(1),subplotmat(2),subplotmat(4));plot(xdim,phase(ydim),'.-');title('Phase plot');
                    end
                end
            end
        end
    case 2   
       if strcmp(Mtype,'Mz')
           subplot(subplotmat(1),subplotmat(2),subplotmat(3));
           h2=imagesc(abs(Showdata));caxis([-1,1]);
           set(h2,'Xdata',z,'Ydata',y);axis([z(1) z(end) y(1) y(end)]);  %Need update
           if length(subplotmat)==4
               subplot(subplotmat(1),subplotmat(2),subplotmat(4));
               h3=imagesc(angle(Showdata));caxis([-pi,pi]);
               set(h3,'Xdata',z,'Ydata',y);axis([z(1) z(end) y(1) y(end)]); %Need update
           end
       else
           subplot(subplotmat(1),subplotmat(2),subplotmat(3));
           h2=imagesc(abs(Showdata));caxis([0,1]);
           set(h2,'Xdata',z,'Ydata',y);axis([z(1) z(end) y(1) y(end)]); %Need update
           if length(subplotmat)==4
               subplot(subplotmat(1),subplotmat(2),subplotmat(4));
               h3=imagesc(angle(Showdata));caxis([-pi,pi]);
               set(h3,'Xdata',z,'Ydata',y);axis([z(1) z(end) y(1) y(end)]); %Need update
           end
       end
    case 3
        imshow3Dfull(abs(Showdata));
        figure();imshow3Dfull(abs(angle(Showdata)));
%         disp('support in the future');
end

    


