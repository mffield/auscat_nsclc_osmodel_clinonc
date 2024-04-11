function pcorr = pcorr_central(model)

%retrieved stats now report/plot
model.client(1).pcorr
model.pcorr=model.client(1).pcorr;
modelfields = fieldnames(model.client(1).pcorr);
for i=2:length(model.client)
    for k=1:length(modelfields)
        model.pcorr.(modelfields{k}) = model.pcorr.(modelfields{k}) + model.client(i).pcorr.(modelfields{k});
    end
end

pcorr = (model.pcorr.n.*model.pcorr.xys - model.pcorr.xs.*model.pcorr.ys)./...
    (sqrt(model.pcorr.n.*model.pcorr.xs2 - model.pcorr.xs.^2).*sqrt(model.pcorr.n.*model.pcorr.ys2 - model.pcorr.ys.^2));
