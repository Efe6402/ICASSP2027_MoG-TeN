function T = append_table(T,newrow)
%APPEND_TABLE Append a table row safely when T may be empty.
if isempty(T)
    T = newrow;
else
    T = [T;newrow];
end
end
