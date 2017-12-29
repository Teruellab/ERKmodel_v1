% Compare MI (across EGF doses) for unpaired and paired-parameter cases

all_species = {'Raf','Phase1','MEK','Phase2','ERK','Phase3','Phase4'};
combos = nchoosek(1:length(all_species),2);


% Split out dose response
for idx = 1:21; % index of MEK-ERK pairing
    k_depth = 5;
    disp(all_species(combos(idx,:)));

    var1 = all_erk_paired;
    erk1 = cell(size(var1{idx},3),1);
    for i = 1:size(var1{idx},3)
        erk1{i} = var1{idx}(:,:,i)';
        erk1{i} = [ prctile(erk1{i},90); prctile(erk1{i},10); prctile(erk1{i},50)]; % Subsample
        %erk1{i} = sum(erk1{i}>5e6);
        %erk1{i} = erk1{i}([1000 3000 5000 7000],:);
    end
    cc_paired = getCC(erk1,k_depth,0);
    disp(['Paired parameter information: ', num2str(cc_paired),' bits'])

    var1 = all_erk_unpaired;
    erk1 = cell(size(var1{idx},3),1);
    for i = 1:size(var1{idx},3)
        erk1{i} = var1{idx}(:,:,i)';
        erk1{i} = [ prctile(erk1{i},90); prctile(erk1{i},10); prctile(erk1{i},50)]; % Subsample
        %erk1{i} = sum(erk1{i}>5e6);
        %erk1{i} = erk1{i}([1000 3000 5000 7000],:);
    end
    cc_unpaired = getCC(erk1,k_depth,0);
    disp(['Unpaired parameter information: ', num2str(cc_unpaired),' bits'])
    cc_diff = log(abs(2^cc_unpaired - 2^cc_paired))/log(2);
    if cc_unpaired>cc_paired
        cc_diff = [0 cc_diff];
    else
        cc_diff = [cc_diff 0];
    end
    disp(['DIFFERENCE: ', num2str(cc_diff)])
    disp('- - - - - - ')
end