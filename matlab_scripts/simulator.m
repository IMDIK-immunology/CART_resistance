
% Script control parameters

    PAUSE_ON_ITERATION = false;
    LOAD_POPULATION = true;

% Additional functions

    function y = linearProbability(left_probability, right_probability, min_value, max_value, sample)
        y = left_probability + (right_probability-left_probability)*((sample-min_value)/(max_value-min_value));
    end

% Parameters related to the presentation of histograms

    bins = 10.^((-100:500)/100);

% Other parameters

    BASE_POPULATION_SIZE = 1000000;
    
    RESISTANT_ANTIGEN_MAX_VALUE = 10;
    LOW_ANTIGEN_MAX_VALUE = 2000;

    % Choose one of the two therapies by commenting accordingly
    LOW_ANTIGEN_DEATH_PROBABILITY = 0.02; % Selective
    % LOW_ANTIGEN_DEATH_PROBABILITY = 0.1; % Non-selective
    HIGH_ANTIGEN_DEATH_PROBABILITY = 0.1;
    
    LOW_ANTIGEN_MUTATION_PROBABILITY = 0.00005;
    HIGH_ANTIGEN_MUTATION_PROBABILITY = 0.00005;
    
    POPULATION_REDUCTION_GOAL = 0.6;

% Creation of a directory with date and time for copied and created files

    data_folder_string_tmp = string(datetime);
    data_folder_string = "Simulation "+strrep(data_folder_string_tmp,':','-');
    mkdir(data_folder_string);

% Copying the currently executed file to this directory

    executed_file_name = string(mfilename);
    copyfile(executed_file_name + ".m",data_folder_string);
    
    copyfile("histogram_generator.m",data_folder_string);

% clear initial_population_antigen;

    initial_population_antigen=[];

% Preparation of the base population

if (LOAD_POPULATION)
    disp("Loading base population from file");

    base_population_folder_string = "Base population 10M 01-Dec-2024 17-33-38";

    base_population_file_fullpath = base_population_folder_string + "/base_population_antigen.mat";
    load(base_population_file_fullpath,'initial_population_antigen');
    copyfile(base_population_file_fullpath, data_folder_string);
    
    base_population_generator_fullpath = base_population_folder_string + "/base_population_generator.m";
    copyfile(base_population_generator_fullpath, data_folder_string);
else
    disp("Generating base population");
    z = randn(BASE_POPULATION_SIZE,1);
    
    % Dodajemy tylko pospólstwo i dwór
    for c = 1:size(z,1)
        new_element = z(c)*3000+9000;
        if (new_element > RESISTANT_ANTIGEN_MAX_VALUE) 
            initial_population_antigen(end+1) = new_element;
        end
    end
    
    initial_population_antigen=initial_population_antigen.';

end

% Displaying of the base population

    h = histogram(initial_population_antigen, bins);
    set(gca, "XScale", "log");
    xlabel('Base population');


% For the purpose of the first iteration, when it is still undefined
population_reduction_real = 0;

for iteration = 1:50

    current_population_antigen = initial_population_antigen;
    
    current_population_size = size(current_population_antigen,1);
    % Code for the proper operation of the while loop below
    previous_population_size = current_population_size +100;
    
    
    % The first iteration is without mutation, killing etc.
    
    if (iteration>1)
    
        % Killing
        
           while ((current_population_size > BASE_POPULATION_SIZE*(1-POPULATION_REDUCTION_GOAL)) ...
                && (current_population_size < (previous_population_size - 10)))
        
                previous_population_size = current_population_size;
            
                clear current_population_antigen;
                current_population_antigen=[];
                
                for c = 1:size(initial_population_antigen,1)
                    if (initial_population_antigen(c) <= RESISTANT_ANTIGEN_MAX_VALUE)
                        if (rand > 0.0)
                            current_population_antigen(end+1) = initial_population_antigen(c);
                         end
                    elseif (initial_population_antigen(c) <= ((LOW_ANTIGEN_MAX_VALUE+RESISTANT_ANTIGEN_MAX_VALUE)/2))
                          if (rand > (LOW_ANTIGEN_DEATH_PROBABILITY))
                             current_population_antigen(end+1) = initial_population_antigen(c);
                         end
                    elseif (initial_population_antigen(c) <= LOW_ANTIGEN_MAX_VALUE)
                          % if (rand > (LOW_ANTIGEN_DEATH_PROBABILITY))
                          if (rand > linearProbability(LOW_ANTIGEN_DEATH_PROBABILITY, HIGH_ANTIGEN_DEATH_PROBABILITY, ...
                                 ((LOW_ANTIGEN_MAX_VALUE+RESISTANT_ANTIGEN_MAX_VALUE)/2), LOW_ANTIGEN_MAX_VALUE, ...
                                 initial_population_antigen(c)))
                            current_population_antigen(end+1) = initial_population_antigen(c);
                         end
                    elseif (rand > (HIGH_ANTIGEN_DEATH_PROBABILITY))
                       current_population_antigen(end+1) = initial_population_antigen(c);
                     end
                end
                
                current_population_antigen=current_population_antigen.';
                current_population_size = size(current_population_antigen,1);
           
                initial_population_antigen=current_population_antigen;
              
            end
        
        population_reduction_real = 1 - (size(current_population_antigen,1)/BASE_POPULATION_SIZE);
        
        
        % Poliferation of individual cells (population growth)
        
            current_population_size = size(current_population_antigen,1);
            
            while (current_population_size < 0.99*BASE_POPULATION_SIZE)
            
                current_population_size = size(current_population_antigen,1);
                
                poliferation_probability = (BASE_POPULATION_SIZE / current_population_size) - 1;
                
                if (poliferation_probability > 1.0)
                    poliferation_probability = 1.0;
                end
                
                for c = 1:current_population_size
                    if ((1-rand) < poliferation_probability)
                           current_population_antigen(end+1) =  current_population_antigen(c);
                    end
                end
            
            end
 
        % Mutation
        
            for c = 1:size(current_population_antigen,1)
              if (current_population_antigen(c) <= RESISTANT_ANTIGEN_MAX_VALUE)
              elseif (current_population_antigen(c) <= LOW_ANTIGEN_MAX_VALUE)
                    if ((1.0-rand) <= LOW_ANTIGEN_MUTATION_PROBABILITY)
                        current_population_antigen(c) = rand*RESISTANT_ANTIGEN_MAX_VALUE;
                   end
              else
                 if ((1-rand) <= HIGH_ANTIGEN_MUTATION_PROBABILITY)
                        current_population_antigen(c) = LOW_ANTIGEN_MAX_VALUE*3/2-rand*LOW_ANTIGEN_MAX_VALUE*4/3;
                        if (current_population_antigen(c) <= RESISTANT_ANTIGEN_MAX_VALUE)
                            current_population_antigen(c) = LOW_ANTIGEN_MAX_VALUE + 1;
                        end 
                 end 
              end
            end
    
    end % End of code excluded for first iteration
    
    % Counting individuals from different groups
    
    resistant_cells_cardinality = 0;
    low_antigen_cell_cardinality = 0;
    high_antigen_cells_cardinality = 0;
    current_population_size=size(current_population_antigen,1);
    
    
    for c = 1:current_population_size
        if (current_population_antigen(c) <= RESISTANT_ANTIGEN_MAX_VALUE)
            resistant_cells_cardinality = resistant_cells_cardinality+1;
        elseif (current_population_antigen(c) <= LOW_ANTIGEN_MAX_VALUE)
            low_antigen_cell_cardinality = low_antigen_cell_cardinality+1;
        else
            high_antigen_cells_cardinality = high_antigen_cells_cardinality+1;
        end
    end
    
    % Saving to the file
    
        iteration_file_fullpath = data_folder_string+"/iteration_" + int2str(iteration) + "_current_population_antigen.mat";
        save(iteration_file_fullpath,'current_population_antigen');
        
        iteration_file_fullpath = data_folder_string+"/iteration_" + int2str(iteration) + "_scalars.mat";
        save(iteration_file_fullpath,'iteration', ...
            'RESISTANT_ANTIGEN_MAX_VALUE', 'LOW_ANTIGEN_MAX_VALUE', ...
            'LOW_ANTIGEN_DEATH_PROBABILITY', 'HIGH_ANTIGEN_DEATH_PROBABILITY', ...
            'LOW_ANTIGEN_MUTATION_PROBABILITY', 'HIGH_ANTIGEN_MUTATION_PROBABILITY', ...
            'resistant_cells_cardinality', 'low_antigen_cell_cardinality', 'high_antigen_cells_cardinality', ...
            'POPULATION_REDUCTION_GOAL', 'population_reduction_real');
    
    % Displaying of the histogram
    
        h = histogram(current_population_antigen, bins);
        set(gca, "XScale", "log");
        set(get(gca,'YLabel'),'Rotation',0);
        
        str = sprintf(['Resistant-A Max Value: %d\n Low-A Max Value: %d\n\n' ...
            'Low-A Death P %f\n High-A Death P %f\n\n' ...
            'Low-A Mutation P %f\n High-A Mutation P %f\n\n' ...
            'resistant-a cells cardinality %d\n low-a cells cardinality %d\n high-a cells cardinality %d\n\n' ...
            'Popul Reduct Goal %f\n popul reduct real %f\n\n'], ...
            RESISTANT_ANTIGEN_MAX_VALUE, LOW_ANTIGEN_MAX_VALUE, ...
            LOW_ANTIGEN_DEATH_PROBABILITY, HIGH_ANTIGEN_DEATH_PROBABILITY, ...
            LOW_ANTIGEN_MUTATION_PROBABILITY, HIGH_ANTIGEN_MUTATION_PROBABILITY, ...
            resistant_cells_cardinality, low_antigen_cell_cardinality, high_antigen_cells_cardinality, ...
            POPULATION_REDUCTION_GOAL, population_reduction_real);
        % title(str);
        ylabel(str);
        
        str = sprintf('iteration %d, population %d', iteration, current_population_size);
        xlabel(str);
        drawnow;
    
    iteration
    
    
    % Assignment initial population antigen from current population antigen
             
        initial_population_antigen=current_population_antigen;
                 
        if (PAUSE_ON_ITERATION)
            disp('PRESS ANY KEY to pass to the next iteration, press CTRL + C to finish the simulation');
            pause
        end
        
end % End of iteration