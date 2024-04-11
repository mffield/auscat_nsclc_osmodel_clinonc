function [CM, acc, predicted_outcome_dist, actual_outcome_dist] = model_scoring(p, y, predict_thresh,calibration_bins)

CM = zeros(2,2,length(predict_thresh));
for k=1:length(predict_thresh)
    py = ones(size(p)); % assign the class labels based on a chosen threshold (default to zero)
    py(p<predict_thresh(k)) = -1;
    py(p>=predict_thresh(k)) = 1;
    acc = sum(py(:)==y(:))/length(y); % compute the prediction accuracy
    [CM(:,:,k), ~] = confusionmat(y(:),py(:),'order',[1 -1]); % compute the confusion matrix
end
predicted_outcome_dist = histcounts(p,calibration_bins);
actual_outcome_dist = histcounts(p(y==1),calibration_bins);

end