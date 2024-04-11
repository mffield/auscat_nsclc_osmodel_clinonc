function [Cstat]=concordance_index(event_times,event_observed,predicted_event_times)


died_mask = logical(event_observed);
%# TODO: is event_times already sorted? That would be nice...
died_truth = event_times(died_mask);
[~,ix] = sort(died_truth); % argsort?
died_truth = died_truth(ix);
temp = predicted_event_times(died_mask);
died_pred = temp(ix);

censored_truth = event_times(~died_mask);
[~,ix] = sort(censored_truth);
censored_truth = censored_truth(ix);
temp = predicted_event_times(~died_mask);
censored_pred = temp(ix);

censored_ix = 1;
died_ix = 1;
times_to_compare = BTree(sort(unique(died_pred)));
num_pairs = 0;
num_correct = 0;
num_tied = 0;
%
%     # we iterate through cases sorted by exit time:
%     # - First, all cases that died at time t0. We add these to the sortedlist of died times.
%     # - Then, all cases that were censored at time t0. We DON'T add these since they are NOT
%     #   comparable to subsequent elements.
while 1
    has_more_censored = censored_ix < length(censored_truth);
    has_more_died = died_ix < length(died_truth);
    %# Should we look at some censored indices next, or died indices?
    if has_more_censored && (~has_more_died || (died_truth(died_ix) > censored_truth(censored_ix)))
        [pairs, correct, tied, next_ix] = handlePairs(censored_truth, censored_pred, censored_ix, times_to_compare);
        censored_ix = next_ix;
    elseif has_more_died && (~has_more_censored || (died_truth(died_ix) <= censored_truth(censored_ix)))
        [pairs, correct, tied, next_ix] = handlePairs(died_truth, died_pred, died_ix, times_to_compare);
        for j=1:length(died_pred(died_ix:next_ix))-1
            pred = died_pred(died_ix:next_ix);
            times_to_compare.insert(pred(j))
        end
        died_ix = next_ix;
    else
        assert(~(has_more_died || has_more_censored))
        break
    end
    
    num_pairs = num_pairs + pairs;
    num_correct = num_correct+ correct;
    num_tied = num_tied + tied;
end

% disp((num_correct + num_tied/2)/num_pairs)
% disp([num_correct, num_tied, num_pairs])
% disp(times_to_compare.counts)
% disp(times_to_compare.tree)

Cstat = (num_correct + num_tied/2)/num_pairs;
end


function [pairs, correct, tied, next_ix]=handlePairs(truth, pred, first_ix, times_to_compare)

next_ix = first_ix;

while (next_ix < length(truth)) && (truth(next_ix) == truth(first_ix))
    next_ix = next_ix + 1;
end
pairs = times_to_compare.len() * (next_ix - first_ix);
correct = 0; tied = 0;
for i=first_ix:(next_ix-1)
    [rank, count] = times_to_compare.rank(pred(i));
    correct = correct + rank;
    tied = tied + count;
end


end


