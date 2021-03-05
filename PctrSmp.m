classdef PctrSmp < handle
properties
    IszRC
    PszRC
    overlapPixRC
    bCPoverlap
    maxNptch
    rndSd

    cpLookup

    IctrInd=cell(2,1);
    count=[0 0];
end
properties(Hidden=true)
    Lookup=cell(2,1)
    LookupTmp=cell(2,1)
end
properties(Access=private)
    bCPoverlapSET=-1 % bCPOverlap set by the end user.  If 0, always becomes bCPoverlap. If 1, bCPoverlap when cpLookup is valid
    kStart           % the actual anchor that gets sampled from first
    bReset=1         % flag for seed recet
    indLookup        % map of indeces for quickly transforming from logical indexing to integer indexing
    blankmap         % ones
    blankLookup      % iintial lookup
    bRnd=0           % flag for randomization
end
methods
    function obj=PctrSmp(IszRC,PszRC, overlapPixRC, bCPoverlap, maxNptch, rndSd)
    % Quickly sample as many patch centers as possible, with or without overlap between patches, from a stereo-image.
    %
    %Constructor:
    %   obj=PctrSmp(IszRC,PszRC, overlapPixRC, cpLookup, bCPoverlap, maxNptch, rndSd)
    %
    %Conventions:
    %   k or K                                  LorR.
    %                                                k=1 -> L
    %                                                k=2 -> R
    %   nK                                      Opposite of LorR
    %                                                if k==1, nK=2
    %                                                if k==2, nK=1
    %   cells                                   Group things into LandR
    %                                               eg = IctrIndl{1} -> IctrIndl for L
    %   fname
    %                                               eg = IctrIndl{2} -> IctrIndl for R
    %
    %IN:
    %   IszRC               int  [1x2]          Image size by rows and columns
    %   PszRC               int  [1x2]          Minimum distance to seperate sample centers by pixel rows and columns
    %   maxNptch            int  [1x1]          Max number of samples to be taken [1 x 1]
    %   overlapPixRC        int  [1x2]          Number of pixels patches are allowed to overlap vertically [1 x 2]
    %
    %   cpLookup      {2x1} num  [nPtch x 2]    Corresponding points, formatted s.t. indices of each cell corresponds to anchoring CP.
    %                                           Input an index
    %                                           For example, assuming a 10x10 image:
    %                                               indeces 1,1 in {1}(1,:) corresponds to L anchor at pixel index 1 (or subscript [1,1]
    %                                                              {1}(1,:) contains R CP (in subscripts) for that anchor point
    %                                               indeces 1,1 in {2}(1,:) corresponds to R anchor at pixel index 1 (or subscript [1,1]
    %                                                              {2}(1,:) contains L CP (in subscripts) for that anchor point
    %                                               indeces 1,1 in {1}(50,:) corresponds to L anchor at pixel index 50 (or subscript [10,5]
    %                                                              {1}(1,:) contains R CP (in subscripts) for that anchor point
    %
    %                                               {2}(1,:) corresponds to R anchor index 1 (or subscript [1,1]
    %   bCPoverlap          bool [1x1]          Whether sampled corresponding patch (eg R) of an anchor patch (L) can overlap with anchored patches when anchor is switched (R)
    %   rndSd               int  [1x1]          Random seed for sampling.
    %                                           0 -> non-random sampling that maximizes number of samples
    %OUT:
    %   IctrIndl      {2x1} int [nPtch x 1]     Indeces sampled for all images
    %   count
    %
    %NOTES:
    %   Constructor does not run sampling -> use obj.run (see: help PctrSmp.run
    %
        obj.IszRC=IszRC;
        obj.PszRC=PszRC;

        obj.indLookup=reshape(1:(obj.IszRC(1)*obj.IszRC(2)),obj.IszRC);
        obj.init_lookup();

        % OVERLAP PXI
        if ~exist('overlapPixRC','var') || isempty(overlapPixRC)
            obj.overlapPixRC=[0 0];
        elseif numel(obj.overlapPixRC) == 1
            obj.overlapPixRC=repmat(overlapPixRC,1,2);
        else
            obj.overlapPixRC=overlapPixRC;
        end

        if exist('bCPoverlap','var') && ~isempty(bCPoverlap)
            obj.bCPoverlapSET=repmat(obj.bCPoverlapSET,1,2);
        else
            obj.bCPoverlapSET=[1 1];
        end

        %maxNptch
        if ~exist('maxNptch','var') || isempty(maxNptch)
            obj.maxNptch=0;
        else
            obj.maxNptch=maxNptch;
        end

        % RND
        if ~exist('rndSd','var') || isempty(rndSd)
            obj.rndSd=0;
        else
            obj.rndSd=rndSd;
        end
        obj.bRnd=logical(obj.rndSd);

    end

    function obj=run(obj,PctrInd,PctrIndPre,K,cpLookup)
    %OPTIONAL IN:
    %     PctrInd    - Indeces of possible potential sample centers points {2 x 1}[N x 1]
    %                       i.e. map of regions that are ok to sample from
    %     PctrIndPre - Indeces of pre-ampled centers points {2 x 1}[N x 1]
    %                       Marked as invalid sampling regions according to overlap rules, but not included in final list of samples
    %     K          - which anchor to start with [1x1] defaults to 1 if no sampling has been done, otherwise defaults to 1 if 2 was ran last, and 2 if 1.
    %                       1 => L
    %                       2 => R
    %     cpLookup   - see constructor help
        if ~exist('PctrInd','var')  || (iscell(PctrInd) && all(cellfun(@isempty,PctrInd)))
            PctrInd=[];
        elseif ~iscell(PctrInd)
            PctrInd={PctrInd,[]};
        end
        if ~exist('PctrIndPre','var')  || (iscell(PctrIndPre) && all(cellfun(@isempty,PctrIndPre)))
            PctrIndPre=[];
        elseif ~iscell(PctrIndPre)
            PctrIndPre={PctrIndPre,[]};
        end

        if ~exist('K','var') || isempty(K)
            K=[];
        end

        if exist('cpLookup','var') && ~isempty(cpLookup)
            obj.set_cpLookup(cpLookup);
        end

        obj.main(PctrInd,PctrIndPre,K);
    end
    function obj=set_cpLookup(obj,cpLookup)
    %  set if new image
        if ~exist('cpLookup','var') || isempty(cpLookup)
            cpLookup={[],[]};
        elseif ~iscell(cpLookup)
            cpLookup={cpLookup,[]};
        end
        % CPLOOKUP
        obj.cpLookup=cpLookup;
        %obj.format_cpLookup(cpLookup);
    end
    function obj=reset_lookup(obj)
    % use if new image or new sampling conditions where overlap with previous samples is ok
        %obj.init_lookup();
        obj.Lookup{1}=obj.blankLookup;
        obj.Lookup{2}=obj.blankLookup;
    end
    function obj=reset_seed(obj)
    % set if new image to maintain reproducibility
        obj.bReset=1;
    end
    function obj=reset_samples(obj)
    % set if new image or new sampling condition in order to keep samples separate
        obj.IctrInd=cell(2,1);
        obj.count=[0 0];
    end
    function plot(obj)
        figure(1);
        imagesc(obj.Lookup{1});

        figure(2);
        imagesc(obj.Lookup{2});
    end
end
methods(Access=private)
%     N,M       - subscripts of valid sampling indeces. Used to apply 1s to surrounding region
%     Lookup    - specifies valid sampling points: valid == 0, invalid > 0
%     bRnd      - equals 1 if random seed (rndSd) is not empty
%                 if equals 1, will randomly select next sampled index
%                      different for each seed, but produces fewer points
%                 if equals 0, will sample first available point.
%                       produces more points, but will be the same each this code is ran
%     valPot    - Lookup, where pseudo logcial indeces==0 are converted to reguar indeces
%                 sampling occurs on from matrix
%     indLookup - indeces for each pixel for quick conversion from logical index
%                 [] -> auto-define within function - this slows down things considerably
%                       if looping over many images or sampling conditions. Define variable
%                       outside these loops if possible!
%% GET
    function obj=init_lookup(obj)
        % Create smpBorder that limits sampling, multiply smpRadius by 3 due to buffering in LRSI sampling

        right = obj.IszRC(2) - obj.PszRC(2)-1;
        left  = obj.PszRC(2) + 1;
        top   = obj.PszRC(1) + 1;
        bot   = obj.IszRC(1) - obj.PszRC(1)-1;

        obj.blankmap=zeros(obj.IszRC);
        obj.blankLookup=ones(obj.IszRC);

        %BORDER
        obj.blankLookup(  1:end,     1:left) = 0; %Left
        obj.blankLookup(  1:end, right:end ) = 0; %Right
        obj.blankLookup(  1:top,     1:end ) = 0; %Top
        obj.blankLookup(bot:end,     1:end ) = 0; %Bottom
        obj.reset_lookup();
    end
%% MAIN
    function obj=main(obj,PctrInd,PctrIndPre,K)

        obj.bCPoverlap=obj.bCPoverlap;
        if isempty(obj.cpLookup{1})
            obj.bCPoverlap(1)=0;
        else
            obj.bCPoverlap(1)=obj.bCPoverlapSET(1);
        end
        if isempty(obj.cpLookup{2})
            obj.bCPoverlap(2)=0;
        else
            obj.bCPoverlap(2)=obj.bCPoverlapSET(2);
        end

        if ~isempty(PctrInd)
            obj.apply_PctrInd(PctrInd);
        end

        if ~isempty(PctrIndPre)
            obj.pre_sample(PctrIndPre);
        end

        if isempty(K) && ~isempty(obj.kStart)
            obj.kStart=PcrSmp.get_nK(obj.kStart);
        elseif isempty(K)
            obj.kStart=1;
        else
            obj.kStart=K;
        end

        obj.sample();


    end
    function obj=apply_PctrInd(obj,PctrInd)

        for k = 1:length(PctrInd)
            if isempty(PctrInd{k}); continue; end

            Pctr=obj.blankmap; % all bad ones
            Pctr(PctrInd{k})=1; % add in good

            %figure(300)
            %subPlot([4 1],1,1)
            %imagesc(Pctr)
            %subPlot([4 1],2,1)
            %imagesc(obj.Lookup{k}>0)

            obj.LookupTmp{k}=(Pctr & obj.Lookup{k}>0);

            %subPlot([4 1],3,1)
            %imagesc(obj.LookupTmp{k})

            %subPlot([4 1],4,1)
            %imagesc(obj.Lookup{k})

            %drawnow
            %waitforbuttonpress

            % FOR DEBUGGING
            %pause(.1)
        end
    end

%% PRE
    function obj=pre_sample(obj,PctrIndPre)
        for k = 1:length(PctrIndPre)
            if isempty(PctrIndPre{k})
                continue
            end
            obj.pre_sample_fun(k,PctrIndPre{k});
        end

    end
    function obj =pre_sample_fun(obj,k,PctrIndPri)
        l=samplesPatchCentersMax3.get_nk(k);

        for i = 1:length(PctrIndPri)
            [N,M]=ind2sub(obj.IszRC,PctrIndPre(i));
            NM=obj.cpLookup{k}(PctrIndPri(i),:);
            [obj.Lookup{k},~,~]=update_lookup(obj.Lookup{k}, 0, N,M,1);
            [obj.Lookup{l},~,~]=update_lookup(obj.Lookup{l}, 0, NM(1),NM(2),1);
        end
    end
%% DAVE'S SPECIAL SAMPLING ALGORITHM
    function obj =sample(obj)

        if obj.bRnd && obj.bReset
            rng(obj.rndSd,'v4');
            obj.bReset=0;
        end

        bStop=[0 0];
        i=obj.kStart;
        while true
            cf=0;
            i=i+1;
            k=(mod(i,2)==0)+1;
            l=(mod(i,2)~=0)+1;
            if all(bStop(k))
                break
            end

            %VALID SAMPLING POINTS
            valPot=obj.indLookup(obj.LookupTmp{k}==1); % ALL VALID SAMPLING INDS % 49.9

            %RUN UNTIL NO MORE LEFT OR maxNptch HAS BEEN SATISFIED
            if isempty(valPot) || (obj.maxNptch>0 && obj.maxNptch<obj.count(k))
                bStop(k)=1;
                continue
            end

            % USE FIRST VALID POINT (L to R,T to B) UNLESS RANDOMIZING
            if obj.bRnd
                ind=randi(length(valPot),1);
            else
                ind=1;
            end
            val=valPot(ind); % index to  sample

            if ~obj.bCPoverlap
                [obj, cf]=obj.update_lookup_bi(val,k,l);
            else
                [N,M]=ind2sub(obj.IszRC,val);
                [obj.Lookup{k}, obj.count(k), cf]=obj.update_lookup(obj.Lookup{k},obj.count(k),N,M,0); %27.5
                obj.LookupTmp{k}=obj.LookupTmp{k} & obj.Lookup{k}; % 21
                %figure(301)
                %imagesc(obj.Lookup{k})
            end

            if cf==1
                obj.Lookup{k}(val)=1; % bad point, out of bounds
                continue
            end
            obj.IctrInd{k}(end+1,1)=val;
        end
    end
    function [vert,hori]=get_inds(obj,N,M)
        %REMOVE VALID POINTS FROM AROUND SAMPLED INDEX
        n(1)=obj.PszRC(1)-obj.overlapPixRC(1)*2;
        n(2)=obj.PszRC(1)-obj.overlapPixRC(1)*2;
        m(1)=obj.PszRC(2)-obj.overlapPixRC(1)*2;
        m(2)=obj.PszRC(2)-obj.overlapPixRC(1)*2;
        vert=cell2mat(arrayfun(@(y) y-n(1):y+n(2), N,'UniformOutput',false));
        hori=cell2mat(arrayfun(@(x) x-m(1):x+m(2), M,'UniformOutput',false));
        vert=round(reshape(vert,numel(vert),1));
        hori=round(reshape(hori,numel(hori),1));
    end
    function [obj,cf]=update_lookup_bi(obj,val,k,l)
        [N,M]=ind2sub(obj.IszRC,val);
        NM=obj.cpLookup{k}(val,:);
        if any(isnan(NM))
            cf=1;
            return
        end
        %cpLookup{1}(10,:) -> (where 10 is IND in A) returns nearest pixel CP in B in RC

        [obj.Lookup{l},obj.count(l),cf]=update_lookup(obj.Lookup{l},obj.count(l),NM(1),NM(2),0);
        if cf; return; end % inds out of bounds
        [obj.Lookup{k},obj.count(k),~ ]=update_lookup(obj.Lookup{k},obj.count(k),N,    M,    0);
        obj.LookupTmp{k}=obj.LookupTmp{k} & obj.Lookup{k};
    end
    function [Lookup,count,continueflag]=update_lookup(obj,Lookup,count,N,M,bNoExit)
        continueflag=0;
        [vert,hori]=obj.get_inds(N,M); % 6.3
        bBad=max(hori) > obj.IszRC(2) || max(vert) > obj.IszRC(1) || min(vert) < 1 || min(hori) < 1 || any(isnan(hori)) || any(isnan(vert));
        if ~bNoExit &&  bBad
            continueflag=1;
            return
        elseif bNoExit  && bBad
            ind=hori > obj.IszRC(2) | vert > obj.IszRC(1) | vert < 1 | hori < 1;
            hori(ind)=[];
            vert(ind)=[];
        end
        Lookup(vert,hori)=Lookup(vert,hori)-1; %93.2
        count=count+1;
    end
%% PLOT
end
methods(Static=true, Access=private)
    function nk=get_nk(k)
        if k == 1
            nk=2;
        elseif k==2
            nk=1;
        end
    end
end
end
