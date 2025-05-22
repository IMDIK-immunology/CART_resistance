% Parameters related to the presentation of histograms

    bins = 10.^((-100:500)/100);

% Other parameters

    BASE_POPULATION_SIZE = 10000000;
    
    RESISTANT_ANTIGEN_MAX_VALUE = 10;
    LOW_ANTIGEN_MAX_VALUE = 1000;

% Creation of a directory with date and time for the copied and created files

    data_folder_string_tmp = string(datetime);
    data_folder_string = "Base population "+strrep(data_folder_string_tmp,':','-');
    mkdir(data_folder_string);

% Copying the currently executed file to this directory

    executed_file_name = string(mfilename);
    copyfile(executed_file_name+".m",data_folder_string);

% Preparation of the base population
 
    z = randn(BASE_POPULATION_SIZE,1);
    
    % clear initial_population_antigen;
    initial_population_antigen=[];
    
    % Additon of only high and low-antigen cells
    for c = 1:size(z,1)
        new_element = z(c)*3000+10000;
        if (new_element > RESISTANT_ANTIGEN_MAX_VALUE) 
            initial_population_antigen(end+1) = new_element;
        end
    end
    
    % Optional addition of the single resistant cell
    initial_population_antigen(end+1) = RESISTANT_ANTIGEN_MAX_VALUE/2;
    
    initial_population_antigen=initial_population_antigen.';

% Displaying of the base population

    h = histogram(initial_population_antigen, bins);
    set(gca, "XScale", "log");
    xlabel('Base population');

% Saving the base population to file

    iteration_file_fullpath=data_folder_string+"/base_population_antigen.mat";
    save(iteration_file_fullpath,'initial_population_antigen');

