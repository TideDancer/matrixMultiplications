% list must be sorted, in ascending order !!!!
% can be implemented using binary search for high performance searching
% temporily use linear search for simplicity
% if search beyond last value, will return -1

function index = findInInterval(list, value);
if sort(list, 'ascend') == list
  if value <= list(1)
    index = 1;
    return;
  else
    for index = 2: length(list)
      if value <= list(index) 
        return;
      end
    end
    index = -1;
  end
end
