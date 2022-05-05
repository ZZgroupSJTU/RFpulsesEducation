function rfstruct=calc_rf(rfstruct,p1pat,p1,flip1,varargin)

fulfilename=['/home/Zhiyongi/vnmrsys/shapelib/',p1pat,'.RF'];
if length(varargin)>=2
    fulfilename=[varargin{2},'\',p1pat,'.RF'];
end

fid = fopen(fulfilename,'r'); 
if fid==-1
    error('Error opening the aglient file');
end

rfstructmember=fieldnames(rfstruct);
for k=1:length(rfstructmember)
    membervalue=rfstruct.(rfstructmember{k});
    if isstruct(membervalue)
        rfstructsubmember=fieldnames(rfstruct.(rfstructmember{k}));
        for m=1:length(rfstructsubmember)
            submembervalue=rfstruct.(rfstructmember{k}).(rfstructsubmember{m});
            parvalue=freadpar(fid,rfstructsubmember{m},submembervalue);
            rfstruct.(rfstructmember{k}).(rfstructsubmember{m})=parvalue;
        end
    elseif isnumeric(membervalue) || ischar(membervalue)
        parvalue=freadpar(fid,rfstructmember{k},membervalue);
        rfstruct.(rfstructmember{k})=parvalue;
    end      
end
fclose(fid);
rfstruct.pulseName=p1pat;
rfstruct.duration=p1;
rfstruct.flip=flip1;

switch rfstruct.pulseName
    case 'sinc'
        rfstruct.bandwidth=5/rfstruct.duration;
    case 'gauss'
        rfstruct.bandwidth=2/rfstruct.duration;
end

if ~isempty(varargin)
    rfstruct.bandwidth=varargin{1};
end
