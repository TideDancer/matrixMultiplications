% list must be sorted, in ascending order !!!!
% implemented using binary search
% if not in range, will return -1

function index = findInInterval(list, value);
if sort(list, 'ascend') == list
  if value <= list(1)
    index = 1;
    return;
  elseif value >= list(end)
    index = length(list);
    return;
  else
    left = 1;
    right = length(list);
    mid = ceil((left+right)/2);
    while value <= list(mid-1) || value > list(mid)
      if value == list(mid-1)
        index = mid - 1;
        return;
      elseif value < list(mid-1)
        right = mid;
      else
        left = mid - 1;
      end
      mid = ceil((left+right)/2);
    end
    index = mid;
  end
else
  'list not in ascend order'
end
