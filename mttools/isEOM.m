function [eom_flag, eom_vec] = isEOM(serial_vec, isLastDayEoM)
%
%.........................................................................
% Create: January 20, 2020 (Masaki Tanaka)
%

%% 1. Check imput
nDay = size(serial_vec, 1);

%Check if the date vector is vartical. If not, transpose the vector.
if nDay == 1
    isHorizon  = 1;
    serial_vec = serial_vec';
    nDay       = size(serial_vec, 1);
else
    isHorizon = 0;
end

% Check if the date vector is given as "serial numbers" or "datetimes"
if isnumeric(serial_vec)
    isSerial = 1;
    date_vec = datetime(serial_vec, 'ConvertFrom','datenum');
elseif isdatetime(serial_vec)
    isSerial   = 0;
    date_vec   = serial_vec;
    serial_vec = datenum(date_vec);
end

if nargin == 1
    % Default: The last of the array is not end of month
    isLastDayEoM = 0;
end


%% 2. Check if each day is the end of months
month_vec     = month(date_vec);
eom_flag      = zeros(nDay, 1);   % Prepare memory for an array
eom_flag(end) = isLastDayEoM;     % Record if the last of the array is end of month

for i = 1 : nDay - 1
    if month_vec(i) ~= month_vec(i + 1)
        % If the next record of a day belongs to another month, the day is
        % recognized to be the end of month.
        eom_flag(i) = 1;
    end 
end

%% 3. Output

% Convert dummy variables to logical variables. 
eom_flag   = logical(eom_flag);
if isSerial
   % If input consists of serial numbers, returen serial numbers 
   eom_vec = serial_vec(eom_flag);
else
   % If input consists of datetimes, returen datetimes 
   eom_vec = date_vec(eom_flag);
end

if isHorizon
    eom_vec  = eom_vec' ;
    eom_flag = eom_flag';
end
end

