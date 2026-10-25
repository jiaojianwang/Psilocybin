function results = permutation_lmm(cfg)
%PERMUTATION_LMM Test a linear session trend with within-subject permutations.
%
% Required cfg fields:
%   inputType   'nifti', 'matrix', or 'array'
%   outputDir   Directory for CSV, MAT, and optional NIfTI outputs
%
% For 'nifti' and 'matrix':
%   dataDir     Directory containing ses01, ses02, ses03, ...
%   subjectIds  Cell array or string array of subject IDs
%   sessions    Session numbers (default: 1:3)
%   atlasFile   Required for 'nifti'; nonzero values identify regions
%
% For 'array':
%   values      Observations-by-features numeric matrix
%   subjects    Subject ID for each observation
%   sessions    Numeric session code for each observation
%   featureNames  Optional feature names
%
% Optional: nPermutations (default 1000), seed (default 1),
%           filePattern (default '*sub-%s-ses%02d*.nii' or .txt;
%                        placeholders are subject ID and session number),
%           writeMaps (default true for NIfTI input).
%
% Each feature uses y ~ session + (1|subject), with numeric session codes
% (1 = baseline, 2 = 24 h, 3 = 1 week by default), maximum-likelihood
% fitting, and an upper-tail F permutation test. Only
% observed visits are shuffled within each subject. The reported p_perm is
% mean(F_permuted >= F_observed), matching the original script exactly.
% The manuscript threshold is uncorrected p_perm < 0.05. BH-FDR and
% Bonferroni values are supplied separately for transparency.
% Codes 1, 2, 3 represent equally spaced visit indices, not elapsed days.
%
% Example:
%   cfg.inputType = 'nifti';
%   cfg.dataDir = 'path/to/mk';
%   cfg.atlasFile = 'path/to/atlas.nii';
%   cfg.subjectIds = {'001', '002', '003'};
%   cfg.outputDir = 'path/to/mk_results';
%   results = permutation_lmm(cfg);
% For 92-by-92 FC text matrices, set inputType to 'matrix', dataDir to the
% matrix directory, and outputDir to a separate results directory.

    cfg = validate_config(cfg);
    switch lower(char(cfg.inputType))
        case 'nifti'
            [values, subjects, sessions, featureNames, labels, node1, node2, atlas] = load_nifti(cfg);
        case 'matrix'
            [values, subjects, sessions, featureNames, labels, node1, node2] = load_matrices(cfg);
            atlas = [];
        case 'array'
            values = double(cfg.values);
            subjects = string(cfg.subjects(:));
            sessions = double(cfg.sessions(:));
            if size(values, 1) ~= numel(subjects) || numel(subjects) ~= numel(sessions)
                error('values, subjects, and sessions must have the same number of observations.');
            end
            featureNames = "feature_" + string((1:size(values, 2))');
            if isfield(cfg, 'featureNames')
                featureNames = string(cfg.featureNames(:));
            end
            if numel(featureNames) ~= size(values, 2)
                error('featureNames must match the number of columns in values.');
            end
            labels = nan(size(values, 2), 1);
            node1 = labels;
            node2 = labels;
            atlas = [];
        otherwise
            error('inputType must be ''nifti'', ''matrix'', or ''array''.');
    end

    nFeatures = size(values, 2);
    fSession = nan(nFeatures, 1);
    betaSession = nan(nFeatures, 1);
    pModel = nan(nFeatures, 1);
    pPerm = nan(nFeatures, 1);
    nObserved = zeros(nFeatures, 1);
    nSubjects = zeros(nFeatures, 1);
    status = repmat("insufficient_data", nFeatures, 1);

    for k = 1:nFeatures
        use = isfinite(values(:, k)) & ~ismissing(subjects) & ~ismissing(sessions);
        y = values(use, k);
        subject = categorical(subjects(use));
        session = sessions(use);
        nObserved(k) = numel(y);
        nSubjects(k) = numel(categories(subject));
        if nObserved(k) < 4 || nSubjects(k) < 2 || numel(unique(session)) < 2
            continue;
        end

        data = table(y, session, subject);
        try
            [fSession(k), pModel(k), betaSession(k)] = session_test(data);
            if ~isfinite(fSession(k))
                error('The observed session F statistic is not finite.');
            end
        catch ME
            status(k) = "model_failed";
            warning('Feature %s: %s', char(featureNames(k)), ME.message);
            continue;
        end

        % Reproducible stream for each feature in this fixed feature order.
        rng(cfg.seed + k - 1, 'twister');
        baseSession = data.session;
        group = findgroups(data.subject);
        groupIndex = cell(max(group), 1);
        for g = 1:numel(groupIndex)
            groupIndex{g} = find(group == g);
        end
        fPermuted = nan(cfg.nPermutations, 1);
        for b = 1:cfg.nPermutations
            permutedSession = baseSession;
            for g = 1:numel(groupIndex)
                ix = groupIndex{g};
                permutedSession(ix) = baseSession(ix(randperm(numel(ix))));
            end
            data.session = permutedSession;
            try
                fPermuted(b) = session_test(data);
            catch ME
                error('Permutation %d failed for feature %s: %s', ...
                    b, char(featureNames(k)), ME.message);
            end
            if ~isfinite(fPermuted(b))
                error('Permutation %d gave a nonfinite F value for feature %s.', ...
                    b, char(featureNames(k)));
            end
        end
        pPerm(k) = mean(fPermuted >= fSession(k));
        status(k) = "ok";
        if mod(k, 10) == 0 || k == nFeatures
            fprintf('Completed %d/%d features.\n', k, nFeatures);
        end
    end

    valid = isfinite(pPerm);
    pBonferroni = nan(nFeatures, 1);
    qFDR = nan(nFeatures, 1);
    pBonferroni(valid) = min(1, nFeatures * pPerm(valid));
    qFDR(valid) = bh_fdr(pPerm(valid), nFeatures);
    results = table(featureNames(:), labels(:), node1(:), node2(:), ...
        nObserved, nSubjects, betaSession, fSession, pModel, pPerm, qFDR, pBonferroni, ...
        pPerm < 0.05, qFDR < 0.05, pBonferroni <= 0.05, status, ...
        'VariableNames', {'feature', 'label', 'node1', 'node2', ...
        'n_observations', 'n_subjects', 'beta_session', 'F_session', 'p_model', 'p_perm', ...
        'q_fdr', 'p_bonferroni', 'nominal_05', 'fdr_05', 'bonferroni_05', 'status'});

    if ~isfolder(cfg.outputDir)
        mkdir(cfg.outputDir);
    end
    writetable(results, fullfile(cfg.outputDir, 'permutation_lmm_results.csv'));
    settings = rmfield_if_present(cfg, {'values', 'subjects', 'featureNames'});
    save(fullfile(cfg.outputDir, 'permutation_lmm_results.mat'), 'results', 'settings');
    if ~isempty(atlas) && cfg.writeMaps
        write_stat_maps(cfg.outputDir, atlas, fSession, pPerm, qFDR, pBonferroni);
    end
end

function cfg = validate_config(cfg)
    if ~isstruct(cfg) || ~isfield(cfg, 'inputType') || ~isfield(cfg, 'outputDir')
        error('cfg must contain inputType and outputDir.');
    end
    if ~isfield(cfg, 'nPermutations'), cfg.nPermutations = 1000; end
    if ~isfield(cfg, 'seed'), cfg.seed = 1; end
    validateattributes(cfg.nPermutations, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(cfg.seed, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
    cfg.outputDir = char(cfg.outputDir);
    mode = lower(char(cfg.inputType));
    if ~ismember(mode, {'nifti', 'matrix', 'array'})
        error('Unknown inputType: %s', mode);
    end
    if strcmp(mode, 'array')
        for name = {'values', 'subjects', 'sessions'}
            if ~isfield(cfg, name{1}), error('Missing cfg.%s.', name{1}); end
        end
        validateattributes(cfg.values, {'numeric'}, {'2d'});
        validateattributes(cfg.sessions, {'numeric'}, {'vector', 'finite'});
    else
        for name = {'dataDir', 'subjectIds'}
            if ~isfield(cfg, name{1}), error('Missing cfg.%s.', name{1}); end
        end
        if ~isfield(cfg, 'sessions'), cfg.sessions = 1:3; end
        validateattributes(cfg.sessions, {'numeric'}, {'vector', 'integer', 'positive'});
        cfg.subjectIds = cellstr(string(cfg.subjectIds(:)));
        if isempty(cfg.subjectIds) || isempty(cfg.sessions)
            error('subjectIds and sessions cannot be empty.');
        end
        if numel(unique(cfg.subjectIds)) ~= numel(cfg.subjectIds) || ...
                numel(unique(cfg.sessions)) ~= numel(cfg.sessions)
            error('subjectIds and sessions must not contain duplicates.');
        end
        if ~isfield(cfg, 'filePattern')
            if strcmp(mode, 'nifti')
                cfg.filePattern = '*sub-%s-ses%02d*.nii';
            else
                cfg.filePattern = '*sub-%s-ses%02d*.txt';
            end
        end
        if strcmp(mode, 'nifti') && ~isfield(cfg, 'atlasFile')
            error('NIfTI input requires cfg.atlasFile.');
        end
    end
    if ~isfield(cfg, 'writeMaps'), cfg.writeMaps = strcmp(mode, 'nifti'); end
end

function [values, subjects, sessions, names, labels, node1, node2, atlas] = load_nifti(cfg)
    atlas = double(niftiread(cfg.atlasFile));
    info = niftiinfo(cfg.atlasFile);
    if numel(size(atlas)) ~= 3
        error('The atlas must be a 3-D NIfTI image.');
    end
    labels = unique(atlas(isfinite(atlas) & atlas ~= 0));
    if isempty(labels), error('The atlas has no nonzero labels.'); end
    regionIndex = arrayfun(@(x) find(atlas == x), labels, 'UniformOutput', false);
    [subjects, sessions] = observation_grid(cfg.subjectIds, cfg.sessions);
    values = nan(numel(subjects), numel(labels));
    nLoaded = 0;
    for row = 1:numel(subjects)
        file = input_file(cfg, subjects(row), sessions(row));
        if file == "", continue; end
        nLoaded = nLoaded + 1;
        volume = double(niftiread(char(file)));
        if ~isequal(size(volume), size(atlas))
            error('Image and atlas dimensions differ: %s', char(file));
        end
        for k = 1:numel(labels)
            voxels = volume(regionIndex{k});
            voxels = voxels(isfinite(voxels));
            if ~isempty(voxels), values(row, k) = mean(voxels); end
        end
    end
    if nLoaded == 0, error('No NIfTI files matched the configured subjects and sessions.'); end
    fprintf('Loaded %d/%d NIfTI scans.\n', nLoaded, numel(subjects));
    names = "roi_" + string(labels(:));
    node1 = nan(numel(labels), 1);
    node2 = node1;
    atlas = struct('image', atlas, 'info', info, 'regionIndex', {regionIndex});
end

function [values, subjects, sessions, names, labels, node1, node2] = load_matrices(cfg)
    [subjects, sessions] = observation_grid(cfg.subjectIds, cfg.sessions);
    values = [];
    node1 = [];
    node2 = [];
    nLoaded = 0;
    for row = 1:numel(subjects)
        file = input_file(cfg, subjects(row), sessions(row));
        if file == "", continue; end
        nLoaded = nLoaded + 1;
        matrix = readmatrix(char(file));
        if ~isnumeric(matrix) || size(matrix, 1) ~= size(matrix, 2)
            error('Expected a square numeric matrix: %s', char(file));
        end
        if isempty(values)
            n = size(matrix, 1);
            node1 = zeros(n * (n - 1) / 2, 1);
            node2 = node1;
            k = 0;
            for i = 2:n
                for j = 1:i-1
                    k = k + 1;
                    node1(k) = i;
                    node2(k) = j;
                end
            end
            index = sub2ind([n, n], node1, node2);
            values = nan(numel(subjects), numel(index));
        elseif size(matrix, 1) ~= n
            error('Matrix size differs across scans: %s', char(file));
        end
        values(row, :) = double(matrix(index)); %#ok<AGROW> preallocated after first scan
    end
    if isempty(values), error('No matrix files matched the configured subjects and sessions.'); end
    fprintf('Loaded %d/%d connectivity matrices.\n', nLoaded, numel(subjects));
    names = "edge_" + string(node1) + "_" + string(node2);
    labels = nan(numel(node1), 1);
end

function [subjects, sessions] = observation_grid(subjectIds, sessionNumbers)
    nSubjects = numel(subjectIds);
    nSessions = numel(sessionNumbers);
    subjects = strings(nSubjects * nSessions, 1);
    sessions = zeros(nSubjects * nSessions, 1);
    for s = 1:nSubjects
        rows = (s - 1) * nSessions + (1:nSessions);
        subjects(rows) = string(subjectIds{s});
        sessions(rows) = sessionNumbers(:);
    end
end

function file = input_file(cfg, subject, session)
    folder = fullfile(cfg.dataDir, sprintf('ses%02d', session));
    matches = dir(folder);
    pattern = ['^', regexptranslate('wildcard', ...
        sprintf(cfg.filePattern, char(subject), session)), '$'];
    matchedName = ~cellfun(@isempty, regexpi({matches.name}, pattern, 'once'));
    matches = matches(~[matches.isdir] & matchedName);
    if numel(matches) > 1
        error('Multiple files match subject %s in %s.', char(subject), folder);
    end
    file = "";
    if numel(matches) == 1
        file = string(fullfile(matches(1).folder, matches(1).name));
    end
end

function [f, p, beta] = session_test(data)
    model = fitlme(data, 'y ~ session + (1|subject)', 'FitMethod', 'ML');
    terms = anova(model, 'DFMethod', 'Residual');
    row = find(strcmp(terms.Term, 'session'), 1);
    if isempty(row), error('The session term is absent from the ANOVA table.'); end
    f = terms.FStat(row);
    p = terms.pValue(row);
    beta = model.Coefficients.Estimate(strcmp(model.Coefficients.Name, 'session'));
end

function q = bh_fdr(p, nTests)
    [sorted, order] = sort(p(:));
    m = numel(sorted);
    adjusted = sorted .* (nTests ./ (1:m)');
    adjusted = flipud(cummin(flipud(adjusted)));
    q = zeros(m, 1);
    q(order) = min(1, adjusted);
end

function write_stat_maps(folder, atlas, f, p, q, pBonferroni)
    info = atlas.info;
    info.Datatype = 'single';
    info.BitsPerPixel = 32;
    masks = {p < 0.05, q < 0.05, pBonferroni <= 0.05};
    names = {'session_F_nominal_05.nii', 'session_F_fdr_05.nii', ...
             'session_F_bonferroni_05.nii'};
    for m = 1:numel(masks)
        image = zeros(size(atlas.image), 'single');
        for k = find(masks{m})'
            image(atlas.regionIndex{k}) = single(f(k));
        end
        niftiwrite(image, fullfile(folder, names{m}), info);
    end
end

function value = rmfield_if_present(value, names)
    for k = 1:numel(names)
        if isfield(value, names{k}), value = rmfield(value, names{k}); end
    end
end
