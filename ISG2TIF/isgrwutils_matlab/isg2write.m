function fname = isg2write(geoid, phi, lam, nodata, meta, region, comments)
% Function to write an isg file (version 2.0)
%
% SYNTAX:
%   fname = isg2write(geoid, phi, lam, nodata, meta, region)
%   fname = isg2write(geoid, phi, lam, nodata, meta, region, comments)
%
% INPUT:
% - geoid      --> matrix (for gridded data) or vector (for sparse data) of
%                  geoid / quasi geoid values [m-ft]
% - phi        --> latitude [deg] / north [m-ft] vector of the knots of the
%                  grid or sparse points 
% - lam        --> longitude [deg] / east [m-ft] vector of the knots of the
%                  grid or sparse points 
% - nodata     --> value of nodata in the input data (it will be converted 
%                  to -9999 in the exported file)
% - meta       --> structure with the metadata values (see below for the
%                  possible field, the optional and not required can be not
%                  present in the input)
% - region     --> string with the name of the region covered by the model
% - comments   --> comment to be inserted befor the header (optional).
%                  Write "\n" in the string for a new line
%
% The possible field of meta, according to the format specifications 
% (www.isgeoid.polimi.it/Geoid/format_specs.html) are:
%   - model_name          [optional]     if not present the region name will be used
%   - model_year,         [optional]     if not present it will be left unknown
%   - model_type,         [mandatory]
%   - data_type',         [mandatory]    
%   - data_format',       [optional]     detected from the data shape, used only if the data are wrongly formatted
%   - data_units',        [mandatory]
%   - data_ordering',     [not required] cannot be controlled by the user
%   - ref_ellipsoid',     [optional]     if not present it will be left unknown
%   - ref_frame',         [optional]     if not present it will be left unknown
%   - height_datum',      [optional]     if not present it will be left unknown (optional for geometric or hybrid, not applicable for gravimetric)
%   - tide_system',       [optional]     if not present it will be left unknown
%   - coord_type',        [mandatory]    
%   - coord_units',       [mandatory]
%   - map_projection',    [optional]     only for projected coordinates (if not present it will be left unknown)
%   - EPSG_code',         [optional]     if not present it will be left unknown
%
% OUTPUT:
% - fname      --> name of the exported file
%
%
%   v1.0
%   Lorenzo Rossi
%   DICA - Politecnico di Milano
%   2020-08-04
%
%
%--------------------------------------------------------------------------
%   isg2write v1.0
%   Copyright (C) 2020 Lorenzo Rossi, Mirko Reguzzoni
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

% -------------------------------------------------------------------------
% metadata
%           keyword                   flag
metakey = {'model name',              1; % [x] if not present use region name
           'model year',              2; % [x] optional
           'model type',              1; % [x] mandatory
           'data type',               1; % [x] mandatory
           'data format',             2; % [x] optional (if not provided, detected from the data shape)
           'data units',              1; % [x] mandatory
           'data ordering',           0; % [x] to be ignored (set by default)
           'ref ellipsoid',           2; % [ ] optional
           'ref frame',               2; % [ ] optional
           'height datum',            4; % [ ] optional for geometric or hybrid, not applicable for gravimetric
           'tide system',             2; % [ ] optional
           'coord type',              1; % [x] mandatory
           'coord units',             1; % [x] mandatory (to be optimized)
           'map projection',          3; % [ ] only for projected coordinates
           'EPSG code',               2};% [ ] optional
% flags:
% - 0   --> ignored field (they are set by reading the data)
% - 1   --> mandatory field
% - 2   --> optional (it will be filled with --- if empty)
% - 3   --> mandatory only for projected coordinate system (otherwise ---)
% - 4   --> mandatory only for geometric and hybrid (otehrwiwise ---)
% - [x] --> indicates that the field has been used in a switch sentence
%           (without relying on the name given in the string)

% extract field names for the input meta structure
metaflag  = cell2mat(metakey(:,2));
metakey   = metakey(:,1);
metafield = cellfun(@(x) strrep(x, ' ', '_'), metakey, 'UniformOutput', false);

% verify if all the fields are given, otherwise add empty missing fields
for i = 1:length(metafield)
    if ~isfield(meta, metafield{i})
        meta.(metafield{i}) = '';
    end
end


% -------------------------------------------------------------------------
% verify the model name
if strcmp(meta.model_name, '')
    meta.model_name = region;
end

% convert the model year to a string (only if necessary)
if ~ischar(meta.model_year)
    meta.model_year = sprintf('%.0f', meta.model_year);
end

% check model type
switch lower(meta.model_type)
    case 'gravimetric'
        mts = 'grav';
        meta.model_type = 'gravimetric';
    case 'hybrid'
        mts = 'hybr';
        meta.model_type = 'hybrid';
    case 'geometric'
        mts = 'geom';
        meta.model_type = 'geometric';
    otherwise
		error('The field model_type must assume one of the following values: gravimetric, hybrid or geometric');
end

% check data type
switch lower(meta.data_type)
    case 'geoid'
        dts = 'G';
        meta.data_type = 'geoid';
    case 'quasi-geoid'
        dts = 'Q';
        meta.data_type = 'quasi-geoid';
    otherwise
		error('The field data_type must assume one of the following values: geoid or quasi-geoid');
end

% check the data format according to the data size (grid vs sparse)
if (size(geoid, 2) == 1 || size(geoid, 1) == 1) && length(geoid) == length(phi) && length(geoid) == length(lam) % data are in sparse point format
    meta.data_format = 'sparse';
    % add warning meta.data_format ignored
elseif size(geoid, 1) == length(phi) && size(geoid, 2) == length(lam)  % data are in grid point format
    meta.data_format = 'grid';
    % add warning meta.data_format ignored
else
    % rely on the user input only if the data shape is not clear (it should
    % not be necessary
    warning('The data_format has not been recognized from the shape of the input')
    switch lower(meta.data_format)
        case 'grid'
            meta.data_format = 'grid';     % to standardize the item
        case 'sparse'
            meta.data_format = 'sparse';   % to standardize the item
    end
end

% check the coordinate type and coordinate units chosen by the user
switch lower(meta.coord_type)
    case 'geodetic'
        meta.coord_type = 'geodetic';          % to standardize the item
		switch lower(meta.coord_units)
			case 'deg'
				meta.coord_units = 'deg';      % to standardize the item
			case 'dms'
				meta.coord_units = 'dms';      % to standardize the item
			otherwise
				error('The field coord_type must assume one of the following values: deg or dms')	
		end
		ndp = 6; % number of decimal places for checking / approximate the grid spacing
    case 'projected'
        meta.coord_type = 'projected';             % to standardize the item
		switch lower(meta.coord_units)
			case 'meters'
				meta.coord_units = 'meters';       % to standardize the item
			case 'feet'
				meta.coord_units = 'feet';         % to standardize the item
			otherwise
				error('The field coord_type must assume one of the following values: meters or feet')	
		end
		ndp = 3; % number of decimal places for checking / approximate the grid spacing
    otherwise
		error('The field coord_type must assume one of the following values: geodetic or projected')
        % meta.coord_type = 'geodetic';  % if no input provided the data are interpreted as geodetic / deg
end


% -------------------------------------------------------------------------
% create the coordinates vectors and verify the input data
switch meta.data_format
    case 'grid'
        % verify nodata and replace with -9999 (default nodata value)
        if isnan(nodata)
            geoid(isnan(geoid)) = -9999;
        else
            geoid(geoid == nodata) = -9999;
        end
        nodata = -9999;
        
        % latitude / north
        dphi = abs(unique(round(diff(phi)*10^(ndp+1)))/10^(ndp+1));
        if length(dphi) > 1
            if std(dphi) < 10^(-ndp-1)  % negligible differences (very small std of the possible dphi)
                dphi = mean(dphi);
            else                        % not negligible differences (ask to the user)
                fprintf('Possible values of dphi [deg] / dnorth [m-ft]: ')
                fprintf('%.*f ', ndp+1, dphi)
                fprintf('\n')
                clear dphi
                dphi = input('dphi [deg] / dnorth [m-ft]? ');
            end
        end
        
        % longitude / east
        dlam = abs(unique(round(diff(lam)*10^(ndp+1)))/10^(ndp+1));
        if length(dlam) > 1
            if std(dlam) < 10^(-ndp-1) % negligible differences (very small std of the possible dlam)
                dlam = mean(dlam);
            else                       % not negligible differences (ask to the user)
                fprintf('Possible values of dlam [deg] / deast [m-ft]: ')
                fprintf('%.*f ', ndp+1, dlam)
                fprintf('\n')
                clear dlam
                dlam = input('dlam [deg] / deast [m-ft]? ');
            end
        end
                
        % define grid limits
        phiMin = min(phi);
        phiMax = max(phi);
        lamMin = min(lam);
        lamMax = max(lam);
        
        % compute grid size
        nr = size(geoid,1);
        nc = size(geoid,2);
        
        % make the vector lam and phi column vector
        if size(lam,1) == 1
            lam = lam';
        end
        if size(phi,1) == 1
            phi = phi';
        end
        
        % check if the phi and lam vectores are in the default direction 
        % (N-to-S, W-to-E) according to the required acucracy
        if sum(round(sort(phi)*10^(ndp+1)) == round(phi(end:-1:1)*10^(ndp+1))) ~= length(phi(:))
            geoid = flipud(geoid);
        end
        if sum(round(sort(lam)*10^(ndp+1)) == round(lam*10^(ndp+1))) ~= length(lam(:))
            geoid = fliplr(geoid);
        end
        % set the correct data_ordering (ignoring the one given by the
        % user)
        meta.data_ordering = 'N-to-S, W-to-E';
        warning('The data_ordering has been set to %s by default. Input from the user was ignored!',  meta.data_ordering);
        
        % check dimensions
        if length(phi) ~= nr
            error('Wrong latitude dimension')
        end
        if length(lam) ~= nc
            error('Wrong longitude dimension')
        end
    case 'sparse'
        % verify nodata and replace with -9999 (default nodata value)
        if isnan(nodata)
            geoid(isnan(geoid)) = -9999;
        else
            geoid(geoid==nodata) = -9999;
        end
        nodata = -9999;
        
        % define the data ordering field
        switch meta.coord_type
            case 'geodetic'
                switch meta.data_type
                    case 'geoid'
                        meta.data_ordering = 'lat, lon, N';          % add a warning of ignored data_ordering field
                    case 'quasi-geoid'
                        meta.data_ordering = 'lat, lon, zeta';       % add a warning of ignored data_ordering field                        
                end
            case 'projected'
                switch meta.data_type
                    case 'geoid'
                        meta.data_ordering = 'east, north, N';       % add a warning of ignored data_ordering field
                    case 'quasi-geoid'
                        meta.data_ordering = 'east, north, zeta';    % add a warning of ignored data_ordering field                        
                end
        end
        warning('The data_ordering has been set to %s by default. Input from the user was ignored!',  meta.data_ordering);
        
        % boundary of the points
        phiMin = min(phi);
        phiMax = max(phi);
        lamMin = min(lam);
        lamMax = max(lam);
        
        % make the vectors lam, phi and geoid column vectors
        if size(geoid,1) == 1
            geoid = geoid';
        end
        if size(lam,1) == 1
            lam = lam';
        end
        if size(phi,1) == 1
            phi = phi';
        end
        % data size
        nr = size(geoid,1);
        nc = 3;
end

% -------------------------------------------------------------------------
% get current date (to append in the file name and to the properties of the
% file)
d = datevec(date);		

% build the filename
fname = sprintf('%s_%s_%s_%s%s_%04.0f%02.0f%02.0f.isg', region, meta.model_year, meta.model_name, mts, dts, d(1), d(2), d(3));

tic;
fprintf('Creating the file %s%s%s\n', pwd, filesep, fname);
fprintf('Writing the header...\n')

try
	fid = fopen(fname,'w+');     % write + read (overwrite if the file exist)
	nh = 0;                      % counter for the number of lines of comments + header
	
    % write the comment section of the file -----------------------------------
	if nargin > 6
		fprintf(fid, sprintf(comments));
		fprintf(fid, '\n\n');
		nh = nh + length(strfind(comments, '\n')) + 2;
	end
	fprintf(fid, 'Created by Matlab isg2write v1.0\n\n');
	nh = nh + 2;
	
	% write the header of the file --------------------------------------------
	fprintf(fid, 'begin_of_head ================================================\n'); nh = nh + 1;
	% textual values of the header
	for i = 1:length(metakey)
		% - 0   --> ignored field
		% - 1   --> mandatory field
		% - 2   --> optional
		% - 3   --> mandatory only for projected coordinate system
		% - 4   --> mandatory only for geometric and hybrid
		if (metaflag(i)==2 && isempty(meta.(metafield{i}))) ||...
		(metaflag(i)==3 && strcmp(meta.coord_type, 'geodetic')) ||...
		(metaflag(i)==4 && strcmp(meta.model_type, 'gravimetric')) ||...
		(metaflag(i)==4 && strcmp(meta.model_type, 'hybrid') && isempty(meta.(metafield{i}))) ||...
		(metaflag(i)==4 && strcmp(meta.model_type, 'geometric') && isempty(meta.(metafield{i})))
			fprintf(fid, '%-14s : %s\n', metakey{i}, '---');
		else
			fprintf(fid, '%-14s : %s\n', metakey{i}, meta.(metafield{i}));
		end
		nh = nh + 1;
	end
	
	% numerical values of the header
	switch meta.coord_type
		case 'geodetic'
			switch meta.coord_units
				case 'deg'
					fprintf(fid, 'lat min        = %11.6f\n',  phiMin); nh = nh + 1;
					fprintf(fid, 'lat max        = %11.6f\n',  phiMax); nh = nh + 1;
					fprintf(fid, 'lon min        = %11.6f\n',  lamMin); nh = nh + 1;
					fprintf(fid, 'lon max        = %11.6f\n',  lamMax); nh = nh + 1;
					switch meta.data_format
						case 'grid'
							fprintf(fid, 'delta lat      = %11.6f\n',  dphi); nh = nh + 1;
							fprintf(fid, 'delta lon      = %11.6f\n',  dlam); nh = nh + 1;
						case 'sparse'
							fprintf(fid, 'delta lat      = ---\n'); nh = nh + 1;
							fprintf(fid, 'delta lon      = ---\n'); nh = nh + 1;
					end
				case 'dms'
					fprintf(fid, 'lat min        = %11s\n',  deg2dms(phiMin)); nh = nh + 1;
					fprintf(fid, 'lat max        = %11s\n',  deg2dms(phiMax)); nh = nh + 1;
					fprintf(fid, 'lon min        = %11s\n',  deg2dms(lamMin)); nh = nh + 1;
					fprintf(fid, 'lon max        = %11s\n',  deg2dms(lamMax)); nh = nh + 1;
					switch meta.data_format
						case 'grid'
							fprintf(fid, 'delta lat      = %11s\n',  deg2dms(dphi)); nh = nh + 1;
							fprintf(fid, 'delta lon      = %11s\n',  deg2dms(dlam)); nh = nh + 1;
						case 'sparse'
							fprintf(fid, 'delta lat      = ---\n'); nh = nh + 1;
							fprintf(fid, 'delta lon      = ---\n'); nh = nh + 1;
					end
			end
		case 'projected'
			fprintf(fid, 'north min      = %11.3f\n',  phiMin); nh = nh + 1;
			fprintf(fid, 'north max      = %11.3f\n',  phiMax); nh = nh + 1;
			fprintf(fid, 'east min       = %11.3f\n',  lamMin); nh = nh + 1;
			fprintf(fid, 'east max       = %11.3f\n',  lamMax); nh = nh + 1;
			switch meta.data_format
				case 'grid'
					fprintf(fid, 'delta north    = %11.3f\n',  dphi); nh = nh + 1;
					fprintf(fid, 'delta east     = %11.3f\n',  dlam); nh = nh + 1;
				case 'sparse'
					fprintf(fid, 'delta north    = ---\n'); nh = nh + 1;
					fprintf(fid, 'delta east     = ---\n'); nh = nh + 1;
			end
	end
	fprintf(fid, 'nrows          = %11.0f\n',  nr); nh = nh + 1;
	fprintf(fid, 'ncols          = %11.0f\n',  nc); nh = nh + 1;
	fprintf(fid, 'nodata         = %11.4f\n',  nodata); nh = nh + 1;
	fprintf(fid, 'creation date  =  %02d/%02d/%04d\n', d(3), d(2), d(1)); nh = nh + 1;
	fprintf(fid, 'ISG format     = %11.1f\n', 2.0); nh = nh + 1;
	fprintf(fid,'end_of_head ==================================================\n'); nh = nh + 1;
	
	% return at the beginning of the file and write the header in the command
	% window
	frewind(fid);
	for i = 1: nh
		fprintf('\t%s\n', fgetl(fid));
	end
	
	% write the geoid / quasi-geoid ---------------------------------------
	switch meta.data_format
		case 'sparse' % the data are sparse
			switch meta.coord_type
				case 'geodetic'
				    % the order is lat, lon, N / zeta and the conversion to 
				    % dms is implemented (if required)
					switch meta.coord_units
						case 'deg'
							fprintf(fid, '%11.6f %11.6f %10.4f\n', [phi, lam, geoid]');
						case 'dms'
							% convert the coordinates
							phi = deg2dms(phi);
							lam = deg2dms(lam);
							geoid = reshape(sprintf('%11.4f', geoid), 11, length(geoid))';
							% write the output
							fprintf(fid, '%s', sprintf([phi, repmat(' ', length(phi), 1), ...
                                lam, repmat(' ', length(lam), 1), ...
								geoid, repmat('\n', length(geoid), 1)]'));
					end
				case 'projected'
				    % the order is east, north, N / zeta
					fprintf(fid, '%11.3f %11.3f %10.4f\n', [lam, phi, geoid]');
			end
		case 'grid'
		% the geoid is a grid stored from N-to-S and W-to-E by default
			fprintf(fid, [repmat('%10.4f ',1,nc-1), '%10.4f', '\n'], geoid');
	end
	fclose(fid);
catch e % if an error occur while writing the file, the file will be closed and the error rethrowed
	fclose(fid);
	rethrow(e)
end

fprintf('File written in %.1f sec\n', toc);
end

function dms = deg2dms(deg)
% To convert a vector of deg into the string dd°mm'ss". If the length of
% the vector is one the 0 in front are ignored, e.g. 0.033333 --> 02'00"
sg = sign(deg);                       % sign of the angle
deg = abs(deg);                       % absolute value
dd = fix(deg);                        % extract the degrees
mm = fix((deg - dd) .* 60);           % extract minutes
ss = (deg - dd - mm / 60) .* 3600;    % extract seconds

mm(round(ss) >= 60) = round(mm(round(ss) >= 60) + 1);
ss(round(ss) >= 60) = round(ss(round(ss) >= 60)) - 60;
mm(round(mm) >= 60) = round(dd(round(mm) >= 60) + 1);
mm(round(mm) >= 60) = round(mm(round(mm) >= 60)) - 60;

if length(deg) == 1
    if dd==0 && mm==0
        dms = sprintf('         %02.0f"', sg*ss);
    elseif dd==0
        dms = sprintf('    %02d''%02.0f"', sg*mm, ss);
    else
        dms = sprintf('%4d°%02d''%02.0f"', sg*dd, mm, ss);
    end
else
    dms = reshape(sprintf('%4d°%02d''%02.0f"',  [sg.*dd, mm, ss]'), 11,  length(deg))';
end
end