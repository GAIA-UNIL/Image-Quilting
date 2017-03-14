function [y,i] = datasample(x,k,varargin)
%DATASAMPLE Randomly sample from data, with or without replacement.
%   Y = DATASAMPLE(DATA,K) returns K observations sampled uniformly at random,
%   with replacement, from the data in DATA.  If DATA is a vector, then Y is a
%   vector containing K elements selected from DATA. If DATA is a matrix, then
%   Y is a matrix containing K rows selected from DATA.  If DATA is an N-D
%   array, DATASAMPLE samples along its first non-singleton dimension.  DATA
%   may be a dataset array. Because the sample is taken with replacement, the
%   observations that DATASAMPLE selects from DATA may be repeated in Y.
%
%   Y = DATASAMPLE(DATA,K,DIM) returns a sample taken along dimension DIM of
%   DATA. For example, if DATA is a matrix and DIM is 2, Y contains a
%   selection of DATA's columns.  If DATA is a dataset array and DIM is 2, Y
%   contains a selection of DATA's variables.  Use DIM to ensure sampling
%   along a specific dimension regardless of whether DATA is a vector, matrix
%   or N-dimensional array.
%
%   Y = DATASAMPLE(DATA,K, 'PARAM1',val1, 'PARAM2',val2, ...) or Y =
%   DATASAMPLE(DATA,K,DIM, 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%   optional parameter name/value pairs to control how DATASAMPLE creates the
%   sample.  Parameters are:
%
%      'Replace' - select the sample with replacement if REPLACE is true (the
%                  default), or without replacement if REPLACE is false.  When
%                  sampling without replacement, the observations that
%                  DATASAMPLE selects from DATA are unique.
%
%      'Weights' - create a weighted sample using the positive weights in
%                  the vector W.
%
%   [Y,I] = DATASAMPLE(...) returns an index vector indicating which values
%   were sampled from DATA.  For example, Y = DATA(I) if DATA is a vector,
%   Y = DATA(I,:) if DATA is a matrix, etc.
%
%   DATASAMPLE uses RANDPERM and RANDI to generate random values and therefore
%   changes the state of MATLAB’s global random number generator.  Control
%   that generator using RNG.
% 
%   Y = DATASAMPLE(S,...) uses the random number stream S for random number
%   generation.
%
%   Examples:
%
%   Draw five unique values from the integers 1:10.
%      y = datasample(1:10,5,'Replace',false)
%
%   Generate a random sequence of the characters ACGT, with replacement,
%   according to specified probabilities.
%      seq = datasample('ACGT',48,'Weights',[0.15 0.35 0.35 0.15])
%
%   Select a random subset of columns from a data matrix.
%      X = randn(10,1000);
%      Y = datasample(X,5,2,'Replace',false)
%
%   Resample observations from a dataset array to create a bootstrap
%   replicate dataset.
%      load hospital
%      y = datasample(hospital,size(hospital,1))
%   
%   Use the second output to sample "in parallel" from two data vectors.
%      x1 = randn(100,1);
%      x2 = randn(100,1);
%      [y1,i] = datasample(x1,10)
%      y2 = x2(i)
%
%   See also RAND, RANDI, RANDPERM, RNG.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2011/04/16 06:41:04 $

nargs = nargin;

% Process the stream argument, if present.
defaultStream = isnumeric(x) || ~isa(x,'RandStream'); % simple test first for speed
if ~defaultStream
    % shift right to drop s from the argument list
    nargs = nargs - 1;
    s = x;
    x = k;
    if nargs > 1
        k = varargin{1};
        varargin(1) = [];
    end
end

if nargs < 2
    error(message('stats:datasample:TooFewInputs'));
elseif ~(isscalar(k) && isnumeric(k) && (k==round(k)) && (k >= 0))
    error(message('stats:datasample:InvalidK'));
end

% Pull out the DIM arg if present, or determine which dimension to use.
if nargs > 2 && ~ischar(varargin{1})
    dim = varargin{1};
    varargin(1) = [];
    nargs = nargs - 1;
else
    dim = find(size(x)~=1, 1); % first non-singleton dimension
    if isempty(dim), dim = 1; end
end
n = size(x,dim);

replace = true;
wgts = [];
if nargs > 2 % avoid parsing args if there are none
    pnames = {'replace' 'weights'};
    dflts =  { replace      wgts };
    [replace,wgts] = internal.stats.parseArgs(pnames,dflts,varargin{:});

    if ~isempty(wgts)
        if ~isvector(wgts) || length(wgts) ~= n
            error(message('stats:datasample:InputSizeMismatch'));
        else
            sumw = sum(wgts);
            if ~(sumw > 0) || ~all(wgts>=0) % catches NaNs
                error(message('stats:datasample:InvalidWeights'));
            end
        end
    end
    if ~isscalar(replace) || ~islogical(replace)
        error(message('stats:datasample:InvalidReplace'));
    end
end

% Sample with replacement
if replace
    if n == 0
        if k == 0
            i = zeros(0,1);
        else
            error(message('stats:datasample:EmptyData'));
        end
        
    elseif isempty(wgts) % unweighted sample
        if defaultStream
            i = randi(n,1,k);
        else
            i = randi(s,n,1,k);
        end
        
    else % weighted sample
        p = wgts(:)' / sumw;
        edges = min([0 cumsum(p)],1); % protect against accumulated round-off
        edges(end) = 1; % get the upper edge exact
        if defaultStream
            [~, i] = histc(rand(1,k),edges);
        else
            [~, i] = histc(rand(s,1,k),edges);
        end
    end
    
    % Sample without replacement
else
    if k > n
        error(message('stats:datasample:SampleTooLarge'));
        
    elseif isempty(wgts) % unweighted sample
        if defaultStream
            i = randperm(n,k);
        else
            i = randperm(s,n,k);
        end
        
    else % weighted sample
        if sum(wgts>0) < k
            error(message('stats:datasample:TooFewPosWeights'));
        end
        if defaultStream
            i = wswor(wgts,k);
        else
            i = wswor(s,wgts,k);
        end
    end
end

% Use the index vector to sample from the data.
if ismatrix(x) % including vectors abd including dataset
    if dim == 1
        y = x(i,:);
    elseif dim == 2
        y = x(:,i);
    else
        reps = [ones(1,dim-1) k];
        y = repmat(x,reps);
    end
else % N-D
    subs = repmat({':'},1,max(ndims(x),dim));
    subs{dim} = i;
    y = x(subs{:});
end
