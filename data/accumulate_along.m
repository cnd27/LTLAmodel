function output = accumulate_along(data,accums,rename_undefineds,accumulated_input)

if nargin==2
    rename_undefineds = false;
    accumulated_input = 1;
elseif nargin==3
    accumulated_input = 1;
end

%% remove undefines and nans
for i=1:length(accums) 
    if iscategorical(data.(accums{i}))
        if rename_undefineds
            data.(accums{i})(isundefined(data.(accums{i}))) = "undefined";
        else
            data(isundefined(data.(accums{i})),:) = [];
        end
    elseif isdatetime(data.(accums{i}))
        if length(accumulated_input)>1
            accumulated_input(isnat(data.(accums{i})),:) = [];        
            accumulated_input(data.(accums{i})<datetime(2020,1,1),:) = []; % remove dates that are too early
        end
        data(isnat(data.(accums{i})),:) = [];        
        data(data.(accums{i})<datetime(2020,1,1),:) = []; % remove dates that are too early
   else
        data(isnan(data.(accums{i})),:) = [];        
    end
end
if length(accumulated_input)>1
    data(isnan(accumulated_input),:) = [];
    accumulated_input(isnan(accumulated_input)) = [];
end

%% Find numbers for everything
[xsize,~] = size(data);
outputs = zeros(length(accums),xsize);
for i=1:length(accums)
    if ~isdatetime(data.(accums{i})) % Don't do this for dates because we want to include dates that don't have data
        [output_names{i},~,outputs(i,:)] = unique(data.(accums{i}));
        output.(accums{i}) = output_names{i};
    else
        outputs(i,:) = datenum(data.(accums{i}))-datenum(2020,1,1)+1;
        dates = datenum(2020,1,1):datenum(max(data.(accums{i})));
        output.(accums{i}) = datetime(dates,'ConvertFrom','datenum')';
    end
end


%% accumulate
accum = accumarray(outputs',accumulated_input);

output.dimensions = accums;
output.accum = accum;