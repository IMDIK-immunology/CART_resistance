
% Creation of a directory with date and time for copied and created files
    
    data_string_tmp = string(datetime);
    data_string = strrep(data_string_tmp,':','-');
    output_data_folder_string = "Histograms " + data_string;
    mkdir(output_data_folder_string);

% Copying the currently executed file to this directory
    executed_file_name = string(mfilename);
    copyfile(executed_file_name + ".m", output_data_folder_string);


% name of the directory with the file with input data
    input_data_folder_string = "Simulation 1M Nonselective 30-May-2025 22-02-46";
    % Simulation 1M Selective 30-May-2025 22-16-38
    % Simulation 1M Nonselective 30-May-2025 22-02-46
    % Simulation 10M Selective 30-May-2025 20-01-23
    % Simulation 10M Nonselective 30-May-2025 20-19-34

% For a population of 1 million set PF=1 , for 10 million set PF=10
PF = 1;
% PF = 10;

for iteration = 1:50

    % Reading data from the file
    
        iteration_file_fullpath = input_data_folder_string + "/iteration_"+int2str(iteration)+"_current_population_antigen.mat";
        load(iteration_file_fullpath,'current_population_antigen');
        
        iteration_file_fullpath = input_data_folder_string + "/iteration_"+int2str(iteration)+"_scalars.mat";
        load(iteration_file_fullpath,'iteration', ...
            'RESISTANT_ANTIGEN_MAX_VALUE', 'LOW_ANTIGEN_MAX_VALUE', ...
            'LOW_ANTIGEN_DEATH_PROBABILITY', 'HIGH_ANTIGEN_DEATH_PROBABILITY', ...
            'LOW_ANTIGEN_TRANSITION_PROBABILITY', 'HIGH_ANTIGEN_TRANSITION_PROBABILITY', ...
            'resistant_cells_cardinality', 'low_antigen_cell_cardinality', 'high_antigen_cells_cardinality', ...
            'POPULATION_REDUCTION_GOAL', 'population_reduction_real');
    
    % Displaying of the histogram
    
    % Parameters related to the presentation of histograms
    
       bins = 10.^((-100:500)/100);
    
    % Data preprocessing

        current_population_size = size(current_population_antigen,1);
        clear pokaz;
        pokaz=[];
        
        % Sigmoid parameters - antigen expression correction
        %    a = 1.0; %amplitude (maximum value)
        %    b = -1; % stepness ("heat" parameter)
        %    c = 12; % x-axis shift (srodkowa iteracja)
            
         %   x = iteration;
            
         %   sigmoid = a / (1+exp(-b*(x-c)));
            
         %   LOWEST_SIGMOID = 1/30;
            
         %   if (sigmoid > (LOWEST_SIGMOID))
         %       mult_factor = sigmoid;
         %   else
         %       mult_factor = LOWEST_SIGMOID;
         %  end


         % if sigmoid is inactive uncomment this line
            mult_factor = 1.0;
        
       for c = 1:current_population_size
             if (current_population_antigen(c) <= RESISTANT_ANTIGEN_MAX_VALUE)
                 mu = 0;
                   sigma = 5;
                  pokaz(end+1) = current_population_antigen(c)+random('Normal',mu,sigma);
             elseif (current_population_antigen(c) <= LOW_ANTIGEN_MAX_VALUE)  
                   sigma = current_population_antigen(c)/1.0;
                   pokaz(end+1) = current_population_antigen(c)+random('Normal',mu,sigma);
             else
                  mu = 0;
                   darken_pupulation_antigen = current_population_antigen(c)*mult_factor;
                   sigma = darken_pupulation_antigen/1.0;
                   pokaz(end+1) = darken_pupulation_antigen+random('Normal',mu,sigma);
              end
         end
    
    
    % Actual presentation of charts
    
        h = histogram(pokaz, bins);
        h.FaceColor = 'black';
        h.EdgeColor = 'black';
        h.LineWidth = 2.0;
        ylim([0 14000*PF]);
        xlim([0.1 200000]);
        set(gca, "XScale", "log");
        set(get(gca,'YLabel'),'Rotation',0);
        
        str = sprintf(['Resistant-A Max Value: %d\n Low-A Max Value: %d\n\n' ...
            'Low-A Death P %f\n High-A Death P %f\n\n' ...
            'Low-A Transition P %f\n High-A Transition P %f\n\n' ...
            'resistant-a cells cardinality %d\n low-a cells cardinality %d\n high-a cells cardinality %d\n\n' ...
            'Popul Reduct Goal %f\n popul reduct real %f\n\n'], ...
            RESISTANT_ANTIGEN_MAX_VALUE, LOW_ANTIGEN_MAX_VALUE, ...
            LOW_ANTIGEN_DEATH_PROBABILITY, HIGH_ANTIGEN_DEATH_PROBABILITY, ...
            LOW_ANTIGEN_TRANSITION_PROBABILITY, HIGH_ANTIGEN_TRANSITION_PROBABILITY, ...
            resistant_cells_cardinality, low_antigen_cell_cardinality, high_antigen_cells_cardinality, ...
            POPULATION_REDUCTION_GOAL, population_reduction_real);
        ylabel(str);
        str = sprintf('iteration %d, population %d', iteration, current_population_size);
        xlabel(str);
        set(gca,'FontSize',15);
        tickXVector=[1 10 100 1000 10000 100000];
        set(gca,'XTick',tickXVector);
        tickYVector=[0 2000*PF  4000*PF 6000*PF 8000*PF 10000*PF 12000*PF 14000*PF];
        set(gca,'YTick',tickYVector);
        drawnow;
        
        iteration_file_fullpath = output_data_folder_string + "/iteration_" + int2str(iteration) + "_histogram_with_legend " + data_string +".png";
        exportgraphics(gca, iteration_file_fullpath);
        
        
        str = sprintf('iteration %d', iteration);
        title(gca,str);
        xlabel("MFI");
        ylabel("Count");
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',90,'VerticalAlignment','bottom')
        set(gca,'FontSize',30);
        tickXVector=[1 10 100 1000 10000 100000];
        set(gca,'XTick',tickXVector);
        tickYVector=[0 2000*PF 4000*PF 6000*PF 8000*PF 10000*PF 12000*PF 14000*PF];
        set(gca,'YTick',tickYVector);
        
        drawnow;
        
        iteration_file_fullpath = output_data_folder_string + "/iteration_" + int2str(iteration) + "_histogram_without_legend " + data_string +".png";
        exportgraphics(gca, iteration_file_fullpath);
    
    iteration

end
