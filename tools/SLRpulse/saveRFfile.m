function saveRFfile(rfstruct,path)
%%
fulfilename=[path,'/',rfstruct.pulseName,'.RF'];
fid = fopen(fulfilename,'w'); 
%---------------write head file--------------------%

headname=fieldnames(rfstruct.head);
fprintf(fid,'# ********************************************************\n');
for k = 1:length(headname);
    headvalue=rfstruct.head.(headname{k});
    if isnumeric(headvalue)
        fprintf(fid,'# %s  %f \n', headname{k},headvalue);
    else 
        fprintf(fid,'# %s  %s \n',headname{k},headvalue);
    end
end
fprintf(fid,'# ********************************************************\n');

fprintf(fid,'\n');

headname=fieldnames(rfstruct);
fprintf(fid,'# ******************RF structure@Zhiyong*******************\n');
for k = 1:length(headname);
    headvalue=rfstruct.(headname{k});
    if isnumeric(headvalue)
        if length(headvalue) ==1
            fprintf(fid,'# %s  %f \n',headname{k},headvalue);
        end
    elseif ischar(headvalue) 
        fprintf(fid,'# %s  %s \n',headname{k},headvalue);        
    end
end
fprintf(fid,'# ********************************************************\n');

%---------------wirte rfdata -----------------------%
rfdata=[rfstruct.phase;rfstruct.amp;ones(1,rfstruct.pts);ones(1,rfstruct.pts)];
fprintf(fid,'%3.3f    %3.3f    %1.1f    %0.0f \n',rfdata);
fclose(fid);