function cellName = mtExcelAddress(row, col)
% convert a numerical coordinate into an Excel cell address
% The first input corresponds the row number and the second the column
% number.
% @e.g.) (1,1) --> A1, (3,2) --> B3
%
%..........................................................................
% Create: July 2, 2019 (Masaki Tanaka, Bank of Japan)
%
%

%% Main code
assert(((row > 0) && (row <= 65536)),...
      'Row number cannot be negative or more than 65536.');
assert(((col > 0) && (col <= 256)),...
      'Column number cannot be negative or more than 65536.');

% Find the alphabet corresponding to the column number
cellName = strcat(cellAlphabet(col), num2str(row));

%% Support function
function alphabetOut = cellAlphabet(col)
% This support function finds the alphabet corresponding to the column
% number.

quotient  = floor((col-1)/26);
remainder = mod(col-1, 26);

if quotient <= 0
    str = char(remainder + 65);
elseif quotient < 26
    str = char(quotient - 1 + 65);
    str = strcat(str, char(remainder + 65));
else
    str = strcat(cellAlphabet(quotient),char(remainder + 65));
end

alphabetOut = str;
end

end