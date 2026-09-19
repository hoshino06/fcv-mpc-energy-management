function configure_models()
% Align lowpath battery with the generated model used by integrated NMPC.
expdir=fileparts(fileparts(mfilename('fullpath')));
mdir=fullfile(expdir,'methods','lowpath');
load_system(fullfile(mdir,'main_lowpath.slx'));
root=sfroot; charts=root.find('-isa','Stateflow.EMChart');
fixed_battery=0; fixed_voltage=0;
for i=1:numel(charts)
    ch=charts(i);
    if ~startsWith(ch.Path,'main_lowpath/'), continue; end
    s=ch.Script;
    if contains(s,'dBatt = Battery(Batt,Ib,0);') && contains(s,'Ib = Ib/n_parallel;')
        s=strrep(s,'Ib = Ib/n_parallel;', '% Battery() already converts pack current to cell current.');
        fixed_battery=fixed_battery+1;
    end
    if contains(s,'Vb = Voc-R0*Ib;')
        s=strrep(s,'Vb = Voc-R0*Ib;  %-Vs-Vf;', 'Vb = Voc-R0*Ib-Vs-Vf;');
        fixed_voltage=fixed_voltage+1;
    end
    ch.Script=s;
end
save_system('main_lowpath'); close_system('main_lowpath',0);
fprintf('Battery alignment: %d current scaling, %d terminal voltage charts repaired.\n', ...
    fixed_battery,fixed_voltage);
source=fileread(fullfile(expdir,'methods','ddp','solve_snapshot.m'));
for N=[5 10]
    model=sprintf('main_ddp_N%d',N);
    load_system(fullfile(expdir,'methods','ddp',[model '.slx']));
    charts=root.find('-isa','Stateflow.EMChart'); count=0;
    for i=1:numel(charts)
        ch=charts(i);
        if startsWith(ch.Path,[model '/']) && contains(ch.Script,'dBatt = Battery(Batt,Ib,0);')
            ch.Script=strrep(ch.Script,'Ib = Ib/params.battery.parallel;', ...
                '% Battery() already converts pack current to cell current.');
        end
        if startsWith(ch.Path,[model '/']) && ...
                (contains(ch.Script,'time_ILQR = toc') || contains(ch.Script,'solve_snapshot('))
            s=strrep(source,'N = params.MPC.N;',sprintf('N = %d;',N));
            % Preserve existing chart port names.
            s=strrep(s,'x_next','x_'); s=strrep(s,'lambda','rammda');
            s=strrep(s,'solve_time','time_ILQR');
            ch.Script=s; count=count+1;
            parent=get_param(ch.Path,'Parent'); ports=get_param(ch.Path,'PortHandles');
            names={'snap_x','snap_warm','snap_pd','snap_nominal'}; input_index=[1 2 3 5];
            for q=1:numel(names)
                dest=[parent '/' names{q}];
                if getSimulinkBlockHandle(dest)<0
                    add_block('simulink/Sinks/To Workspace',dest,'VariableName',names{q}, ...
                        'SaveFormat','Timeseries','SampleTime','params.MPC.dt', ...
                        'MaxDataPoints','inf','Position',[1000 100+q*50 1120 125+q*50]);
                    line=get_param(ports.Inport(input_index(q)),'Line');
                    src=get_param(line,'SrcPortHandle'); target=get_param(dest,'PortHandles');
                    add_line(parent,src,target.Inport,'autorouting','on');
                end
            end
        end
    end
    assert(count==1,'Expected exactly one active timed DDP chart');
    save_system(model); close_system(model,0);
end
end
