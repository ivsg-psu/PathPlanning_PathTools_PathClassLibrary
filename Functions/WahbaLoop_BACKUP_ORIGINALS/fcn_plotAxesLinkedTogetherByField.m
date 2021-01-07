function fcn_plotAxesLinkedTogetherByField
%fcn_plotAxesLinkedTogetherByField  Link all handles together
% Revision history:
% 2019_11_25 - first write of code by sbrennan@psu.edu
% 2019_11_26 - changed zoom type to xy, not just x. 

all_h = findobj;
num_axes = 0;
field_numbers = [];
for i_handle = 1:length(all_h)
    if strcmp(all_h(i_handle).Type,'figure')        
        num_axes = num_axes+1;
        field_numbers = [field_numbers;  mod(all_h(i_handle).Number,100)]; %#ok<AGROW>        
    end
end

field_numbers_no_repeats = unique(field_numbers);
for i_fieldNum = 1:length(field_numbers_no_repeats)
    current_field_number = field_numbers_no_repeats(i_fieldNum);
    num_axes = 0;
    linked_axes = [];
    for i_handle = 1:length(all_h)
        if strcmp(all_h(i_handle).Type,'figure')
            num_axes = num_axes+1;
            if current_field_number == mod(all_h(i_handle).Number,100)
                figure(all_h(i_handle).Number);
                linked_axes(num_axes) = gca;     
            end % Ends if statement for current_field_number
        end % Ends if statement for if current handle is a figure
    end % Ends for loop over all open handles
end % Ends for loop over each field number

linkaxes(linked_axes,'xy');
return

