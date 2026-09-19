function configure_model()
% Reproducible, idempotent repairs to the migrated model via Stateflow API.
here=fileparts(mfilename('fullpath'));
model='main_lowpath';
load_system(fullfile(here,[model '.slx']));
root=sfroot;
charts=root.find('-isa','Stateflow.EMChart');
patched=0;
for i=1:numel(charts)
    ch=charts(i);
    if strcmp(ch.Path,[model '/MATLAB Function5'])
        % Restore Iwai's pack-current calculation. The cell current is
        % obtained downstream by dividing by the parallel string count.
        ch.Script=strrep(ch.Script,'Ib_bar = min(max(Pb_bar/Vb,-36),36);', ...
            'Ib_bar = Pb_bar/Vb;');
        patched=patched+1;
        continue
    end
    if ~strcmp(ch.Path,[model '/NMPC']), continue; end
    s=ch.Script;
    if ~contains(s,'forwardpass_runge'), continue; end
    s=strrep(s,'x_bar(1,:)/0.1173','x_bar(1,ii)/0.1173');
    s=strrep(s,'du(:,ii));','du(:,ii),Ref);');
    s=strrep(s,'du(:,ii-1));','du(:,ii-1),Ref);');
    ch.Script=s;
    patched=patched+1;
    % The controller is fixed: four air states, one voltage, ten steps.
    % Make the fixed contract explicit instead of relying on size inference.
    data=ch.find('-isa','Stateflow.Data');
    sizes=struct('Ifc','1','z_ref','1','x','[4 1]','u_bar','[1 10]', ...
        'Pd','[10 1]','du1','1','K0','[1 4]','Smax','1', ...
        'u_b','[1 10]','ub1','1','rammda','[7 11]','sannso','1');
    for j=1:numel(data)
        d=data(j);
        if isfield(sizes,d.Name)
            d.Props.Array.Size=sizes.(d.Name);
            d.Props.Array.IsDynamic=false;
        end
    end
end
assert(patched==2,'Expected the battery-current and active NMPC charts');
save_system(model);
close_system(model,0);
end
