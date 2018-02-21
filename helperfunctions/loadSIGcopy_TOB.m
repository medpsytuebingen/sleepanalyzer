function varargout = loadSIGcopy_TOB(filename, varargin)
% loadSIGcopy statt loadSIG
% 
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% CRITICAL DEBUG: Berechnung des Startblocks korrigiert.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% ACHTUNG, dies ist eine Kopie! Das Original Liegt als private Funktion zur 
% Klasse datafile im entsprechenden Ordner
%
% WICHTIG: ein neues fclose!
% entscheidende Änderung: der Header wurde verändert. Unterschiede treten
% nur intern auf, es sollten keine Probleme Auftreten, das die Erweiterung
% des Ablagefile verändert wurde.
% fast Faktor 10 Zeitersparnis!
% leider noch weit entfernt von "direktem" einlesen. das liegt v.a. an der
% blockkweisen Konstruktion (-> Skip-Wert), somit müssen viel mehr Werte
% als benötigt von der Platte gezogen werden. Ich denke, das ist nur mit
% erhöhtem aufwand (oder erneutem Abspeichern der Daten in einzelen
% Kanälen)  möglich.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% lade daten aus .STS/.SIG Dateien. Als Unterfunktion für Klasse datafile
% gedacht.
%
% [dInfoS, header] = loadSIG(filename, [iChan])
%   liefert dInfoS, ein struct mit den entscheidenden Fileinformationen
%   (die in Klasse datafile benötigt werden) zurück, sowie den "vollen"
%   header der Datei.
%   Optional können die Namen / Namensanfänge der gesuchten Kanäle als cell
%   array iChan angegeben werden
%
% data = loadSIG(filename, chan, [dLen], [dStart])
%   liefert die in chan angegebenen Kanaldaten zurück. optional kann die
%   Anzahl der Datenpunkte pro Kanal begrenzt (dLen) sowie der Startpunkt
%   festgelegt (dStart) werden.
%
% See also
%   DATAFILE
%
% tk, feb. mmxi

% version 2011,10: tmod externalisiert.

% defaults
  headFlag = 0;
  hotChans = {'TL0', 'TL1', 'TR0', 'TR1'};
  refChans = {'Cb1', 'Cb2'};
  
% filename parsing
  if ~ischar(filename)
    error('datafile:rSIGin1a', 'filename has to be a string')
  end
  [iPath, iName, iExt] = fileparts( filename );
  if all(~strcmpi({'.STS', '.SIG'}, iExt))
    error('datafile:rSIGin1b', 'wrong file type')
  end
  if ~exist(fullfile(iPath, [iName '.STS']),'file') || ~exist(fullfile(iPath, [iName '.SIG']),'file')
    error('datafile:rSIGin1c', 'file %s.STS or %s.SIG not found', fullfile(iPath, iName), fullfile(iPath, iName))
  end

% distinguish header reading form data reading
  if nargin == 1 
    headFlag  = 1;
  elseif ischar(varargin{1}) || iscellstr(varargin{1})
    headFlag  = 1;
    hotChans  = varargin{1};
    if ischar(hotChans)
      hotChans = {hotChans};
    end
  end

% collect header information
  [~, dt, Info, header] = fuCreateSigHead( filename );
  df.name = filename;  % name of original file
  df.type = '.SIG';   % type of data
  df.chans= numel(Info); % number of channels in file

  refFlag   = cell( size(refChans) );
  hot       = cell( size(hotChans) );
  for c1 = 1:numel( refChans )
    refFlag{c1} = find(strncmp(refChans{c1}, header.Head.rnames, numel(refChans{c1})));
  end % get reference channels' numbers
  for c1 = 1:numel( hotChans )
    hot{c1}     = find(strncmp(hotChans{c1}, header.Head.rnames, numel(hotChans{c1})));
  end % get 'hot' channels' numbers
  refFlag   = horzcat(refFlag{:});
  hot       = horzcat(hot{:});
  
  if ~all(refFlag==3|refFlag==4)
    warning('datafile:rSIGref','Data was referenced to %s and %s, rereferencing to %s and %s', ...
              header.Head.rnames{refFlag}, refChans{:})
  else
    refFlag = [];
  end % check for correct referencing channels

  df.hot  = hot;     % indices of 'interesting' channels
  df.dt   = dt;      % scanning interval
  df.head = header.Head;   % file / header information
  df.tt   = header.tt;     % time covered by measurement [s]
%   df.ttf  = header.ttf;
%   df.ttc  = header.ttc;
  
  if headFlag
    varargout{1} = df;
    return
  end % 

%%%%%%%%%%%%%%%%%%%%%%%%
% read data for real
%%%%%%%%%%%%%%%%%%%%%%%%

% get interesting channels
  if isnumeric(varargin{1}) && ~isempty(varargin{1})
    oChan = floor(varargin{1}(:));
    oChan = tmod(oChan, df.chans);
  else
    oChan = 1;
  end
  
% get reading length  
  if nargin > 2 && isnumeric(varargin{2}) && ~isempty(varargin{2})
    dLen = max( floor(varargin{2}(1)), 1);
  else
    dLen = header.tt / df.dt;
  end

% get points to skip
  if nargin > 3 && isnumeric(varargin{3}) && ~isempty(varargin{3})
    dStart = max( floor(varargin{3}(1)), 0);
  else
    dStart = 0;
  end
  
% % scaling always inactive at the moment
% %   reasons: reduces performance + no scaling factor available
%   if scaleFlag
%     handleRead  = @(fid, N, rType , rSkip ) fread(fid, N, rType , rSkip ) * sFactor;
%   else
  handleRead  = @(fid, N, rType , rSkip ) fread(fid, N, rType , rSkip );
%   end
  
  fid = fopen( fullfile( iPath, [iName '.SIG']) ,'r');
  
  offset        = mod( dStart, header.bLen );
  N             = dLen + offset;
  
  ausg = NaN( dLen, numel(oChan) );
 
  rChan = [oChan,refFlag];
  
  startblock    = ceil( (dStart + 1)/ header.bLen ) -1;

  for c1 = 1:numel(rChan)
    
    % seek starting block
    fseek(fid, startblock* header.bSize + header.cBeg(rChan(c1)), 'bof' );

    % use matlab's build-in ability to skip datapoints while reading
    rType         = sprintf('%g*%s', header.bLen, header.dType);
    rSkip         = header.bSize - header.bLen * precision(header.dType);

    % read & remove unwanted data that was read from the block's
    % beginning
    dummy     = handleRead(fid, N, rType , rSkip );
    ausg(:,c1)= dummy( 1+offset : end );

  end
  
  fclose(fid);
  
  % rerefence Data if needed
  if ~isempty(refFlag)
    fprintf('data rereferenced to %s \n', refChans{:}) 
    dummy = sum( ausg(:,numel(oChan)+1:end) ,2 ) / numel(refFlag);
    ausg  = ausg( :, 1:numel(oChan) ) - repmat( dummy, [ numel(oChan), 1] );
  end
    
  
  varargout{1} = ausg;
  varargout{2} = df;
  varargout{3} = hotChans;

end  


%% Unterfunktionen
function arr = fuReadNArr(fid, varargin)
% reads arrays that come as (int16)=N, N*(intx). x = 8 per default
%
% arr = fuReadNArr(fid) reads the next int8 array, the first entry is the
%                       array's length
% arr = fuReadNArr(fid, format) as above but different format (int16, int32
%                       supported atm)
% arr = readNarr(fid, format, M), as above but M times (lengths may differ,
%                       thus return value is cell array)
%
% this structure is used in STS files,
% See also:
%   fuCreateSigHead, fuReadCStr, fuReadVStr

  form = 'int8';
  niter = 1;

  if nargin > 1
      if  ~isempty(varargin{1}) && isa(varargin{1},'char')
          switch lower(varargin{1})
              case {'int8', 'int16', 'int32', 'char'}
                  form = lower(varargin{1});
              otherwise
                  error('Stellate:ReadArrays', 'provided format is not supported');
          end
      end
  end

  if nargin > 2
      if  ~isempty(varargin{2}) && isnumeric(varargin{2}) && floor(varargin{2}(1)) > 0
          niter = varargin{2}(1);
      end
  end

  if niter == 1
    alen        = fread( fid, 1, 'int16');
    if alen < 1
      arr = [];
      return
    end
    arr = fread( fid, alen, form );

  else

    arr{niter} = {};
    for cou = 1: niter
      arr{cou} = [];
      alen        = fread( fid, 1, 'int16');
      if alen < 1
        arr(cou:end) = [];
        if isempty(arr)
          arr = [];
        end
        return
      end
      arr{cou} = fread( fid, alen, form );
    end

  end

end
%%
function fstr = fuReadVStr(fid, varargin)
% read strings from file, lenght is given by a int8 number heading the
% string
% 
% fstr = fuReadVStr(fid) reads a string fstr from the file defined by fid.
%           the string's lenght L is given by the next value to be read
%           from the file (thus L+1 bytes will be read)
% fstr = fuReadVStr(fid, N) as above but reads N strings consecutively. fstr
%           is a cell containing the strings
%
% this structure is used in STS files,
% See also:
%   fuCreateSigHead, fuReadCStr, fuReadNArr

  if ftell(fid) < 0
      error('Stellate:IOError1', 'file handle was passed improperly')
  end

  if      nargin > 1 && isnumeric( varargin{1}(1) ) && lower(varargin{1}(1)) > 0
      niter   = lower(varargin{1}(1));
  elseif  nargin > 1 && isnumeric( varargin{1}(1) ) && lower(varargin{1}(1)) == 0
      fstr    = [];
      return
  else
      niter   = 1;
  end
  odata{niter} = [];

  for cou = 1 : niter
      slen        = fread( fid, 1, 'uint8');
      % nothing read / end of file
      if isempty( slen )
          fstr  = [];
          return
      end
      if slen == 255 % more than max (int: -1)
          slen = fread(fid, 1, 'int16');
          if slen == 511
              fseek(fid, -2, 'cof');
              fseek(fid, -1, 'cof');
%             warning('no more  items')
              odata{cou}  = [];
              break
          else
              odata{cou}  = char(fread( fid, slen, 'char')');
          end
      elseif slen < 1
          odata{cou}  = [];
      else
          odata{cou}  = char(fread( fid, slen, 'char')');
      end
  end

  if niter == 1
      fstr = odata{1};
  else
      fstr = odata;
  end
end
%%
function fstr = fuReadCStr(fid, varargin)
% read special strings from file, lenght is given by a int16 number heading
% the string. Also, the string also begins with 'C', which is not returned.
% 
% fstr = fuReadCStr(fid) reads a string fstr from the file defined by fid.
%           the string's lenght L is given by the next value to be read
%           from the file (thus L+3 bytes will be read)
% fstr = fuReadCStr(fid, N) as above but reads N strings consecutively. fstr
%           is a cell containing the strings
%
% this structure is used in STS files,
% See also:
%   fuCreateSigHead, fuReadVStr, fuReadVStr


  if ftell(fid) < 0
      error('Stellate:IOError1', 'file handle was passed improperly')
  end

  if nargin > 1 && isnumeric( varargin{1}(1) ) && lower(varargin{1}(1)) > 0
      niter = lower(varargin{1}(1));
  else
      niter = 1;
  end

  odata{niter}  = [];
  for cou = 1 : niter
      slen        = fread(fid, 1, 'int16');
      if slen
          fseek(fid, 1, 'cof'); % skip the heading 'C'
          odata{cou}  = char(fread(fid, slen-1, 'char')');
      else
          odata{cou}  = [];
      end
  end

  if niter == 1
      fstr = odata{1};
  else
      fstr = odata;
  end
end
%% Header für Stellate-files (.SIG + .STS) anlegen
function  [tt, dt, Info, header] = fuCreateSigHead( file, varargin )
% soll automatisiert Hilfsdatei f�r SIG-Dateien erstellen.
%
% out ist 1 wenn file erstellt wurde,
%         0 sonst

  neuExt  = '.sigHeadNeu';
  
  [fPath, fName, ~] = fileparts(file);

  stellateHeader    = fullfile( fPath, [ fName, '.STS' ] );
  myStellateHeader  = fullfile( fPath, [ fName, neuExt ] );
  datFile           = fullfile( fPath, [ fName, '.SIG' ] );
  

  % checking if input file are valid
  if exist( myStellateHeader, 'file' )
    load( myStellateHeader, '-mat', 'tt', 'dt', 'Info', 'header');
    return
  elseif ~exist( stellateHeader, 'file' )
      error('klustaana:SigInfoWrongFileType', ...
            'file %s not found', file )
  elseif ~exist( file, 'file' )
      error( 'Klusta:SigFileNotFound','file %s not found', file )
  end

  dType   = 'int16';
  bLen    = 64;                   % kanaleinträge pro block %andere blocklänge als 64 möglich?
  dSize   = precision( dType );   % Byte
  hSize   = 8;                    % Byte, keine Ahnung, was drin ist.

  Head    = [];

  % read propriate header file in order to create our own
  fid       = fopen( stellateHeader, 'r' );
% gesamtlänge des Files bestimmen
  fseek(fid, 0, 'eof'); 
  lof       = ftell(fid); 
  fseek(fid, 0, 'bof');

% read distinguishable Header parts
  sfDatHeader
  sfDatPatInfo
  sfDatCodec
  sfDatCalibration
  sfDatPhysicalMontage
  sfDatRecordingMontage
  sfDatReformatedMontage
  sfDatStatusGroup
  sfDatStatusItem

  fclose(fid);
  
  % anzahl datenpunkte / blöcke ermitteln
  fid   = fopen( datFile, 'r' );
  fseek( fid, 0, 'eof');
  Nall  = ftell( fid );
  fclose(fid);

  bSize = hSize + channels * bLen * dSize;
  dpb   = channels * bLen;
  dt    = 1 / Abtastfrequenz;
  tt    = Nall / bSize / Abtastfrequenz * bLen; %[s]
% %   die Daten scheinen Zuverlässig in ganzen Böcken abgelegt zu werden,
% %   somi die die folgenden beiden Zeilen nicht notwendig
%   ttf   = floor(Nall/ bSize) / Abtastfrequenz * bLen; %[s] nur volle Blöcke
%   ttc   = ttf + (rem(Nall,bSize)~=0) * ( rem(Nall , bSize) - hSize ) / (bSize -hSize) / Abtastfrequenz * bLen;
  
  % define header and Info structs
  Info(numel(channels)).scale = NaN;
  
  for cC = (1 : channels)
    Info(cC).scale  = scaling(cC);          % hier! : naja, wo ist echte Skalierung verzeichnet?
    Info(cC).title  = channames{cC};
  end
%     header{cC}.Beg    = (cC -1) * bLen * dSize + hSize; % channels offset from data start (in data points!)
  header.cBeg   = (0 : channels-1) * bLen * dSize + hSize;
  header.tt     = tt;
  header.bSize  = bSize;
  header.bLen   = bLen;
  header.dType  = dType;
  header.dpb    = dpb;
  header.Head   = Head; % mehr infos als eigentlich notwendig
  
%   save( fullfile( fPath, [ fName, neuExt ] ), '-mat', 'tt', 'dt', 'Info', 'header');
  save( myStellateHeader, '-mat', 'tt', 'dt', 'Info', 'header');
% speichern ist noch langsamer als parsen!
  
  % subfunctions either reusable stuff or lots to read
  % most that is read is not used.
  function ssfOZCheck
    OZ = fread(fid, 2, 'int16');
    if any( OZ ~= [-1; 1] )
      error('Stellate: file does not specifications')% -1 1
    end
  end

  function sfDatHeader
%     ssfOZCheck
    if any(fread(fid, 3, 'int16') ~= [1, -1, 1]' ) || ...
        ~strcmp( fuReadCStr(fid), 'FileHeader' )
      fclose( fid );
      error('Stellate:wrongFileType', 'file %s does not meet specifications for .STS', file)
    end
    fread(fid, 1, 'int16');
    fuReadVStr(fid, 'char');        % Pat & institution
    fuReadVStr(fid, 'char');        % Pat & institution
    % schlauer machen!
%     fuReadNArr(fid, 'char');        % Pat & institution
    fread(fid, 5, 'int16');         % cryptic numbers (recording time included?)
    fuReadVStr(fid);                % devicename
  end
  
  function sfDatPatInfo
    ssfOZCheck
    fuReadCStr(fid);                % 'PatientInfo', could be verified
    fread(fid, 1, 'int16');         % another number (2)
    Head.neu = fuReadVStr(fid,4);              % various Pat.data including name
    fread(fid, 5, 'int16');         % cryptic numbers

    dreck = fread(fid, 8, 'int8');  % cryptic numbers
    if sum( dreck )
        warning('Stellate:FileIrregular1','irregular file 1');
        fseek(fid, -8, 'cof');
    end
    Head.Felder = fuReadVStr(fid, 400); % ChannelNames
  end

  function sfDatCodec
    ssfOZCheck
    fuReadCStr(fid);                % 'StdCoDec', could be verified
    fread(fid, 1, 'int16');         % another number (1)
    Head.NChans = fread(fid, 1, 'int16'); % number of channels
    channels    = Head.NChans;
    Head.bsize  = fread(fid, 1, 'int16'); % blocksize
    fread(fid, Head.NChans + 1, 'int16'); % bytes per channel? but nchan + 1
  end

  function sfDatCalibration
    ssfOZCheck
    fuReadCStr(fid);                % 'Calibration', could be verified
    fread(fid,  1, 'int16');        % another number (2)
    fread(fid, 10, 'int16');        % cryptic numbers. last modified?
    nc = fread(fid, 1, 'int16');    % number of channels
    fread(fid, nc * 6, 'int32');    % some calibration?
  end

  function  sfDatPhysicalMontage
    ssfOZCheck
    fuReadCStr(fid);                % 'PhysicalMontage', could be verified
    fread(fid, 1, 'int16');         % another number (4)
    fuReadVStr(fid);                % ableitung invasiv oder nicht?
    Head.pnam   = fread(fid, 2, 'int16')';    % number of possible channel names? (why twice?)
    Head.pnames = fuReadVStr(fid, Head.pnam(1)); % channel names

    fread(fid, 1, 'int16');         % another number (1)
    dum = fread(fid, 1, 'int16');   % 3 * nnam
    fread(fid, dum * 2, 'int32');   % what? == lots of data (6 *nnam *int32)
    fread(fid, 1, 'int16');         % another number (1)
  end

  function sfDatRecordingMontage
    ssfOZCheck
    fuReadCStr(fid);                % 'RecordingMontage', could be verified
    fread(fid, 1, 'int16');         % another number (1)
    Head.rM = fuReadVStr(fid);      % name Montage 1
    fread(fid, 1, 'int16');         % montChans?
    fuReadNArr(fid, 'int16', 5);    % 5 * montChans?
    scaling =fuReadNArr(fid, 'int32', 1);    % montChans?, scaling?

    Head.rx = fread(fid, 1, 'int16'); % another number (10)
    Head.rn = fread(fid, 1, 'int16');	% trace-bytes per block?
    Head.rf = fread(fid, 1, 'int16');	% 1000; scanning freq?
    Abtastfrequenz  = Head.rf;

    Head.rnam = fread(fid, 1, 'int16');   % n chan-names
    Head.rnames = fuReadVStr(fid, Head.rnam); % channel names
    channames   = Head.rnames;
    fuReadNArr(fid, 'int8', 1);     % montChans?. int8! content?
    fuReadNArr(fid, 'int16', 1);    % montChans?. content?
    fuReadNArr(fid, 'int8', 1);     % montChans?. int8! content?
    fread(fid, 1, 'int16');         % another number (1)
  end

  function sfDatReformatedMontage
    ssfOZCheck
    fuReadCStr(fid);                % 'ReformatedMontage', could be verified

    fread(fid, 1, 'int16');         % another number (4)
    fuReadVStr(fid);                % name Montage
    fread(fid, 1, 'int16');         % montChans
    fuReadNArr(fid, 'int16', 5);    % 5 * montChans?: scanning rate/ 2*small / 2*0
    fuReadNArr(fid, 'int32', 1);    % montChans?, scaling?

    fread(fid, 1, 'int16');         % another number (12)
    dum = fread(fid, 1, 'int16');   % 
    fuReadVStr(fid, dum);           % channel names
    dum = fread(fid, 1, 'int16');   % 
    fuReadVStr(fid, dum);           % channel names
    dum = fread(fid, 1, 'int16');   % another number (1)
    fuReadVStr(fid, dum);           % channel name short
    dum = fread(fid, 1, 'int16');   % another number (1)
    fuReadVStr(fid, dum);           % channel name long

    % things about the channel names
    % 1: channel names
    % 2: names again !?!
    % 3: describing the averaged channel
    % 4: describing the averaged channel

    xcl = fread(fid, 4, 'int8')';	% cryptic (74 -120 0 0)
    %%%%%% schluessel des Problems wahrscheinlich in den ersten beiden zaehlern der vorigen
    %%%%%% Reihe
    if xcl(3) == 4
      warning('Stellate:FileIrregular2','reading yet another montage')

      fuReadVStr(fid); % 2nd Montage name
  %     nc = 
      fread(fid, 1, 'int16');
      fuReadNArr(fid,'int16', 5);
      fuReadNArr(fid,'int32', 1);

      fread(fid, 1, 'int16');     % another number

      dum = fread(fid, 1, 'int16')';
      fuReadVStr(fid, dum);
      dum = fread(fid, 1, 'int16')';
      fuReadVStr(fid, dum);
      dum = fread(fid, 1, 'int16')';
      fuReadVStr(fid, dum);
      dum = fread(fid, 1, 'int16')';
      fuReadVStr(fid, dum);

      fread(fid, 4, 'int8');          % cryptic
    end
    fread(fid, 1, 'int16');         % another number
  end

  function sfDatStatusGroup
    ssfOZCheck
    fuReadCStr(fid);                  % 'StatusGroup', could be verified

    % read StatusGroups
    for couSG = 1: 9
      fread(fid, 1, 'int16');         % another number (3)
      fuReadVStr(fid,3);                %
      fread(fid, 1, 'int8');          % single number
      fread(fid, 19,'int8');          % numbers
      fread(fid, 1, 'int8');          % single number
      fread(fid, 3, 'int8');          % numbers
      fuReadVStr(fid);
      dum = fread(fid, 1, 'int16');   % possible Nstrings
      fuReadVStr(fid,dum);            %
      dum = fread(fid, 1, 'int16');   % possible Nstrings
      fuReadVStr(fid,dum);            %
      dum = fread(fid, 1, 'int16');   % possible Nstrings
      fuReadVStr(fid,dum);            %
      dum = fread(fid, 1, 'int16');   % possible Nstrings
      fuReadVStr(fid,dum);            %
      dum = fread(fid, 1, 'int16');   % possible Nstrings
      fuReadVStr(fid,dum);            %
      fread(fid, 2, 'int8');          % usually (15  -128)
    end
  end

  function sfDatStatusItem
    ssfOZCheck
    fuReadCStr(fid);                 % 'StatusItem', could be verified

    while lof > ftell(fid)
      fread(fid, 4, 'int16');         % 1 17/20 ? ?
      fread(fid, 10, 'int16');        % numbers various numbers, why 12? timestamp?
      fread(fid, 2, 'int16');         % 0 -1 !!!
      fuReadVStr(fid);
      fread(fid, 3, 'int8');          % 0 0 0 why?
      fread(fid, 1, 'int16');         % cryptic (15  -128) (sometimes (25 -128)?!?)
    end
  end

end
%% f�r Datentyp ben�tigte anzahl bytes zur�ckliefern
function  nBytes = precision( instr )
  if ~ischar( instr ) || isempty( instr )
    error('klustaana:PrecisionWrongInput', 'input to precision has to be a string' )
  end
  
  switch lower( instr )
    case 'int8'
      nBytes = 2;
    case 'int16'
      nBytes = 2;
    case 'int32'
      nBytes = 4;
    case 'int64'
      nBytes = 8;
    case 'uint8'
      nBytes = 2;
    case 'uint16'
      nBytes = 2;
    case 'uint32'
      nBytes = 4;
    case 'uint64'
      nBytes = 8;
    case 'float32'
      nBytes = 4;
    case 'char'
      nBytes = 1;
    otherwise
      error('klustaana:unknownPrecision', 'precision %s is  not defined', instr )
  end
end
function   out = tmod(B, d)
% modifizierte Version von mod(a,b). 
% Der Wertebereich ist nun [1:b] im Gegensatz zu [0:b-1] bei mod.
%
% TK, MMXI, Okt

  out = mod(B,d) + ~mod(B,d).*d;
end