function [geoid, phi, lam, meta, comm] = isgread(fname, verbose)
% Function to load geoid/quasi-geoid models in ISG format up to version 2.0
%
% SYNTAX:
%      [geoid, phiGrid, lamGrid, meta, comm] = isgread(fname)
%
% INPUT:
% - fname   --> string containing the filename (full or relative path)
% - verbose --> boolean, true = verbose output. If no specified it is
%               assumed true
%
% OUTPUT:
% - geoid   --> matrix (gridded data) or vector (sparse) containing the 
%               geoid / quasi-geoid model [m-ft]
% - phi     --> vector containing latitude [deg] / north [m-ft] values of
%               the knots of the grid (gridded) or of the points (sparse)
% - lam     --> vector containing longitude [deg] / east [m-ft] values of
%               the knots of the grid (gridded) or of the points (sparse)
% - meta    --> structure containing the metadata information (only textual
%               values)
% - comm    --> string containing the comments before the header
%
%   v2.0
%   Lorenzo Rossi
%   DICA - Politecnico di Milano
%   2020-08-03
%
%--------------------------------------------------------------------------
%   isgread v2.0
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

% key of the isg v 2.0 format
metakey = {'model name';
           'model year';
           'model type';
           'data type';
           'data format';
           'data units';
           'data ordering';
           'ref ellipsoid';
           'ref frame';
           'height datum';
           'tide system';
           'coord type';
           'coord units';
           'map projection';
           'EPSG code'};
       
tic
if nargin == 1
    verbose = 1;
end

% verify if the file exist, othervise check if the string is provided
% without the extension
if exist(fname,'file')
    fID = fopen(fname,'r');     % open the file as read only
else
    % if file is not found verify if fname contains the etension, otherwise
    % add it
    if ~strcmp(fname(end-3:end),'.isg') && ~strcmp(fname(end-3:end),'.ISG')
        fname = [fname '.isg'];
        fID = fopen(fname,'r');     % open the file as read only
    else
        msg ='File does not exist!';
        error(msg);
    end
end

try % to close the file if errors occours
    hn = 1;                     % counter of number of header rows
    hr = fgetl(fID);            % string updated with content of each row of the header
    if verbose
        fprintf('Reading header of %s:\n', fname);
        fprintf('\t%s\n',hr);       % write the header line hn
    end
    
    % initialize the output metadata structure
    meta = struct();
    geom = struct();
    comm = '';
    
    while ischar(hr) && ~strncmp(hr, 'begin_of_head', length('begin_of_head'))
        comm = [comm hr '\n'];
        hr = fgetl(fID);          % read the new line of header
        if verbose
            fprintf('\t%s\n',hr); % write the line hn on the screen
        end
    end
    
    hr = fgetl(fID);          % read the new line of header (to skip_begin_of_head)
    
    while ischar(hr) && ~strncmp(hr,'end_of_head', length('end_of_head'))
        idc = strfind(hr,':');    % find ":" symbol in the read string
        ide = strfind(hr,'=');    % find "=" symbol in the read string
        if ~isempty(idc) % read header info whne there is ":"
            meta.(strrep(strtrim(hr(1:idc-2)), ' ', '_')) = strtrim(hr(idc+2:end));
        end
        if ~isempty(ide) % read header info when there is "="
            geom.(strrep(strtrim(hr(1:ide-2)), ' ', '_')) = strtrim(hr(ide+2:end));
        end
        hr = fgetl(fID);          % read the new line of header
        if verbose
            fprintf('\t%s\n',hr); % write the line hn on the screen
        end
        hn = hn + 1;              % increase number of header line conunter
    end
    
    % check the isg version of the loaded file. If the ISG format key is not
    % found in the ISG header, the version is assumed to be 1.0
    if  ~isfield(geom, 'ISG_format')
        warning('No version of the format, version 1.0 is assumed!\n')
        geom.ISG_format = '1.0';
    else
        if str2double(geom.ISG_format) > 2
            error('This function can be used to read ISG files up to version 2.0 of ISG format')
        end
    end
    
    % if last rows is not a char some errors in reading file happens
    if ~ischar(hr) 
        error('Error while reading the header of the file');
    end
    
    ndp             = 6;                                                                                           % number of decimal places (in case of geodetic coordinates) 
    geomfields      = {'lat_min', 'lon_min', 'lat_max', 'lon_max', 'delta_lat', 'delta_lon', 'nrows', 'ncols', 'nodata'};               % mandatory geometry fields (standard)
    geomfields_proj = {'north_min', 'north_min', 'east_max', 'east_max', 'delta_north', 'delta_east', 'nrows', 'ncols', 'nodata'};      % mandatory geometry fields (projected)
    
    % add the format version to the metadata
    meta.ISG_format = geom.ISG_format;
    
    % check the format version
    switch str2double(geom.ISG_format)% v 1.0 and v 1.01
        case {1, 1.01}
            % convert to number the information present in the header and
            % check if they were all provided
            for i = 1:length(geomfields)
                if isfield(geom, geomfields{i})
                    geom.(geomfields{i}) = str2double(geom.(geomfields{i}));
                else
                    error('The mandatory field %s is missing in the header of the file', geomfields{i});
                end
            end
            
            % compute grid spacing
            Dphi1 = (geom.lat_max - geom.lat_min) / geom.nrows;    % built latitude and longitude step
            Dlam1 = (geom.lon_max - geom.lon_min) / geom.ncols;
            if round(Dphi1*10^ndp) ~= round(geom.delta_lat*10^ndp) % verify latitude step up to ndp decimal place (computed vs value in the header)
                error('Error in latitude grid');
            end
            if round(Dlam1*10^ndp) ~= round(geom.delta_lon*10^ndp) % verify longitude step up to ndp decimal place (computed vs value in the header)
                error('Error in longitude grid');
            end
            % define the vectors of latitude and longitude of the grid
            phi   = (geom.lat_max-Dphi1/2 : -Dphi1 : geom.lat_min)';
            lam   = (geom.lon_min+Dlam1/2 :  Dlam1 : geom.lon_max)';
            
            % load geoid
            geoid = cell2mat(textscan(fID, repmat('%f', 1, geom.ncols),'Delimiter',' ','MultipleDelimsAsOne',1)); 
            geoid(geoid == geom.nodata) = nan; % replace nodata with nan
            
            % rename meta fields (according to v2.0 convention)
            meta.data_units = meta.units;
            meta = rmfield(meta, 'units');
            meta.ref_ellipsoid = meta.reference;
            meta = rmfield(meta, 'reference');
            
            % add missing fields (according to v2.0 convention)
            meta.data_ordering = 'N-to-S, W-to-E';
            meta.coord_type    = 'geodetic';
            meta.coord_units   = 'deg';
            meta.data_format   = 'grid';
            dd = fname(end-11:end-4);
            if ~isnan(str2double(dd))
                meta.creation_date = [dd(7:8) '/' dd(5:6) '/' dd(1:4)];
            end
            
        case 2 % v 2.0
            % store information of projected grid into the values
            % of the geodetic grid (i.e. deast will be stored in dlam)
            if strcmp(meta.coord_type, 'projected')
                for i = 1:length(geomfields)
                    if isfield(geom, geomfields_proj{i})
                        geom.(geomfields{i}) = geom.(geomfields_proj{i});
                        geom = rmfield(geom, geomfields_proj{i});
                    end
                end
                ndp = 3;
            end
            % convert to number the information present in the header
            % and check if they were all provided
            for i = 1:length(geomfields)
                if isfield(geom, geomfields{i})
                    if strcmp(meta.coord_units, 'dms') && i <= 6 % convert dms to deg (if required, for the first 6 fields)
                        geom.(geomfields{i}) = dms2deg(geom.(geomfields{i}));
                    else
                        geom.(geomfields{i}) = str2double(geom.(geomfields{i}));
                    end
                else
                    error('The mandatory field %s is missing in the header of the file', geomfields{i});
                end
            end
            switch meta.data_format
                case 'grid'               
                    % compute grid spacing
                    Dphi1 = (geom.lat_max - geom.lat_min) / (geom.nrows-1);  % build latitude / north step
                    Dlam1 = (geom.lon_max - geom.lon_min) / (geom.ncols-1);  % build longitude / east step
                    if round(Dphi1*10^ndp) ~= round(geom.delta_lat*10^ndp)   % verify latitude / north step up to ndp decimal place (computed vs value in the header)
                        error('Error in latitude grid');
                    end
                    if round(Dlam1*10^ndp) ~= round(geom.delta_lon*10^ndp)   % verify longitude / east step up to ndp decimal place (computed vs value in the header)
                        error('Error in longitude grid');
                    end
                    
                    % define the vectors of latitude and longitude of the grid
                    switch meta.data_ordering
                        case 'N-to-S, W-to-E'
                            phi   = (geom.lat_max : -Dphi1 : geom.lat_min)';
                            lam   = (geom.lon_min :  Dlam1 : geom.lon_max)';
                        case 'S-to-N, W-to-E'
                            phi   = (geom.lat_min :  Dphi1 : geom.lat_max)';
                            lam   = (geom.lon_min :  Dlam1 : geom.lon_max)';
                        case 'N-to-S, E-to-W'
                            phi   = (geom.lat_max : -Dphi1 : geom.lat_min)';
                            lam   = (geom.lon_max : -Dlam1 : geom.lon_min)';
                        case 'S-to-N, E-to-W'
                            phi   = (geom.lat_min :  Dphi1 : geom.lat_max)';
                            lam   = (geom.lon_max : -Dlam1 : geom.lon_min)';
                    end
                    
                    % load geoid model
                    geoid = cell2mat(textscan(fID, repmat('%f', 1, geom.ncols),'Delimiter',' ','MultipleDelimsAsOne',1)); 
                    geoid(geoid == geom.nodata) = nan; % replace nodata with nan
                    
                case 'sparse'
                    % load geoid model (as string to manage dms)
                    data = textscan(fID, repmat('%s', 1, geom.ncols), 'Delimiter',' ','MultipleDelimsAsOne',1);
                    
                    % reorder columns according to "data ordering" 
                    idx_comma = [0, strfind(meta.data_ordering, ','), length(meta.data_ordering)+1]; % find "," in the field (to separate the values)
                    if length(idx_comma) == 4
                        for i = 1:length(idx_comma)-1
                            switch strtrim(meta.data_ordering(idx_comma(i)+1:idx_comma(i+1)-1))
                                case {'lon', 'east'}
                                    if strcmp(meta.coord_units, 'dms')
                                        lam = cell2mat(cellfun(@(x) dms2deg(x), data{:,i}, 'UniformOutput', false));
                                    else
                                        lam = cell2mat(cellfun(@(x) str2double(x), data{:,i}, 'UniformOutput', false));
                                    end
                                case {'lat', 'north'}
                                    if strcmp(meta.coord_units, 'dms')
                                        phi = cell2mat(cellfun(@(x) dms2deg(x), data{:,i}, 'UniformOutput', false));
                                    else
                                        phi = cell2mat(cellfun(@(x) str2double(x), data{:,i}, 'UniformOutput', false));
                                    end
                                case {'N', 'zeta'}
                                    geoid = cell2mat(cellfun(@(x) str2double(x), data{:,i}, 'UniformOutput', false));
                            end
                        end
                    else
                        error('Wrong definition of the field "data ordering". The number of column must be 3');
                    end
            end
            % add information about the format and the date to the metadata
            meta.ISG_format = geom.ISG_format;
            meta.creation_date = geom.creation_date;
        otherwise
            error('The format %s is not supported by this function', geom.ISG_format);
    end
    fclose(fID); % close the file
catch e
    fclose(fID);
    rethrow(e);
end

% add missing fields (according to v 2.0 format conventions)
for i = 1:length(metakey)
    if ~isfield(meta, strrep(metakey{i}, ' ', '_'))
        meta.(strrep(metakey{i}, ' ', '_')) = '---';
    end
end
% if verbose mode is on, write on the screen some information about the
% loaded file
if verbose
    fprintf('\nFile %s was correctly read in %.1f sec\n\n', fname, toc);
    fprintf('\tmodel name:\t\t\t\t%s\n', meta.model_name);
    fprintf('\tmodel type:\t\t\t\t%s\n', meta.model_type);
    fprintf('\tdata units:\t\t\t\t%s\n', meta.data_units);
    fprintf('\tref ellipsoid:\t\t\t%s\n', meta.ref_ellipsoid);
    fprintf('\tgeoid size (lat x lon):\t%d x %d\n', size(geoid,1), size(geoid,2));
    fprintf('\tmin geoid latitude:\t\t%.10f\n', min(phi(:)));
    fprintf('\tmax geoid latitude:\t\t%.10f\n', max(phi(:)));
    fprintf('\tmin geoid longitude:\t%.10f\n', min(lam(:)));
    fprintf('\tmax geoid longitude:\t%.10f\n', max(lam(:)));
end

end

function deg = dms2deg(dms)
% To convert a vector of string dd°mm'ss" into deg. It accepts as input
% only one string to be converted

idx_dd = strfind(dms, '°');
idx_mm = strfind(dms, '''');
idx_ss = strfind(dms, '"');

% deg
if isempty(idx_dd)
    dd = 0;
    idx_dd = 0;
    sg = 1;
else
    dd = str2double(dms(1:idx_dd-1));
    sg = sign(dd+1e-16);
    dd = abs(dd);
end

% arc-min
if isempty(idx_mm)
    mm = 0;
    idx_mm = 0;
else
    mm = str2double(dms(idx_dd+1:idx_mm-1));
    sg = sg * sign(mm+1e-16);
    mm = abs(mm);
end

% arc-sec
ss = str2double(dms(idx_mm+1:idx_ss-1));
sg = sg * sign(ss+1e-16);
ss = abs(ss);
deg = (dd + mm / 60 + ss / 3600) * sg; 
end