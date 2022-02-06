clear all
close all
clc


chls = [ "R" , "T", "Z"];   % Define channels 
r = 0.1;                    % Threshold  
r_str = num2str(r); 
delta_v = 'p-0.2'

rawdata_dir = strcat('../data/', delta_v, '/processed/')

%m = ["1.5", "8"; "2", "4"; "2", "5"; "2", "6"; "2", "7"; "2", "8"; "2.5", "6"; "2.5", "7"; "2.5", "8"];
m = ["1", "2"; 
     "2", "4"; 
     "2", "6"; 
     "2", "8"]

loop_len = size(m);

for mfp_i = 1:loop_len(1);
    
    rad      = m(mfp_i,1);
    mfp      = m(mfp_i,2);
    sim_name = strcat(delta_v, '_2hz_', mfp, '_mfp_', rad,'_rad')

    
    
    
    % Create output file file:
    folder_path = strcat('./MSE/', delta_v ,'/', mfp,'_', rad);
    mkdir(folder_path)
    
    for ch = 1:3;
      
        channel = chls(ch);

        input_path = strcat();
        output_name = strcat();

        
        disp('Loading data from')
        disp(input_path)
        t = textread(input_path); 
        disp('Load completed.') 


        out = zeros(50, 11);

        for i=1:11;
            t1 = t(:,i);
            t2 = t1/(max(abs(t1)));
            c1 = coarsegrain_notime(t2, 30);


              %tolerance threshold 35
            scale = 50; %scale factor for decimation


            e1 = movavg_mse(c1, r, scale);

            out(:,i) = e1;


        end

        writematrix(out,output_name,'Delimiter',',')
        disp('Completed:')
        disp(output_name)
        disp('_______________________________________')
        
    end
end 

disp('************************************************')
disp('**************  FINISHED RUN   ***************')
disp('************************************************')
