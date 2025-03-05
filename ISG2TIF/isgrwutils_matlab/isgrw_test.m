%--------------------------------------------------------------------------
%  Copyright (C) 2020 Lorenzo Rossi, Mirko Reguzzoni
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

%% ------------------------------------------------------------------------
% see Example 1 on the specifications of version 2.0
% (http://www.isgeoid.polimi.it/Geoid/format_specs.html)

% string of comments
comm = ['These data are freely available under a Creative Commons Attribution 4.0\n', ...
        'International Licence (CC BY 4.0)\n', ...
        '\n', ...
        'When using the data, please cite:\n', ...
        'A. Name1, B. Name2 (year). Title. Version 1.0. GFZ Data Services.\n', ...
        'http://doi.org/10.5880/isg.2020.001\n', ...
        '\n', ...
        'The original data were provided by C. Name3 (email of dd/mm/yyyy to ISG).\n', ...
        'The present file is distributed by ISG.\n', ...
        '\n', ...
        'This is an example.\n', ...
        'Here some information about model computation can be provided.\n', ...
        '\n', ...
        'Bibliographic reference:\n', ...
        'D. Name4, E. Name5 (year). Title. Journal, Volume(Number), pp. xxx-yyy.'];

% cell array with the header information
metakey = {'model name',              'EXAMPLE';
           'model year',              2020;
           'model type',              'gravimetric';
           'data type',               'geoid';
           'data units',              'meters';
           'data format',             'grid';             % to be ignored by the function
           'data ordering',           'N-to-S, W-to-E';   % to be ignored by the function
           'ref ellipsoid',           'GRS80';
           'ref frame',               'ITRF2014';
           'height datum',            '';
           'tide system',             'mean-tide';
           'coord type',              'geodetic';
           'coord units',             'dms';
           'map projection',          '';
           'EPSG code',               '7912'};

% convert the cell to structure
meta = cell2struct(metakey(:,2), cellfun(@(x) strrep(x, ' ', '_'), metakey(:,1), 'UniformOutput', false), 1);

% geoid data
geoid = [30.1234 31.2222 32.3456 33.4444 34.5678 36.6666;
         41.1111 42.2345 43.3333 44.4567 45.5555 46.6789;
         51.4321 52.9753 53.6543 54.8642     nan     nan;
         61.9999 62.8888 63.7777 64.6666     nan     nan];

% coordinates definition
phiGrid = (41  : -20/60 : 40)';
lamGrid = (120 :  20/60 : 121+40/60)';

% write the file
fname = isg2write(geoid, phiGrid, lamGrid, nan, meta, 'TEST', comm);

% reload the file
[geoid_red, phiGrid_read, lamGrid_read, meta_read, comm_read] = isgread(fname);


%% ------------------------------------------------------------------------
% see Example 2 on the specifications of version 2.0
% (http://www.isgeoid.polimi.it/Geoid/format_specs.html)

% string of comments
comm = ['These data are freely available under a Creative Commons Attribution 4.0\n', ...
        'International Licence (CC BY 4.0)\n', ...
        '\n', ...
        'When using the data, please cite:\n', ...
        'A. Name1, B. Name2 (year). Title. Version 1.0. GFZ Data Services.\n', ...
        'http://doi.org/10.5880/isg.2020.001\n', ...
        '\n', ...
        'The original data were provided by C. Name3 (email of dd/mm/yyyy to ISG).\n', ...
        'The present file is distributed by ISG.\n', ...
        '\n', ...
        'This is an example.\n', ...
        'Here some information about model computation can be provided.\n', ...
        '\n', ...
        'Bibliographic reference:\n', ...
        'D. Name4, E. Name5 (year). Title. Journal, Volume(Number), pp. xxx-yyy.'];

% cell array with the header information
metakey = {'model name',              'EXAMPLE';
           'model year',              '2020';
           'model type',              'gravimetric';
           'data type',               'geoid';
           'data units',              'meters';
           'data format',             'grid';             % to be ignored by the function
           'data ordering',           'N-to-S, W-to-E';   % to be ignored by the function
           'ref ellipsoid',           'GRS80';
           'ref frame',               'ITRF2014';
           'height datum',            '';
           'tide system',             'mean-tide';
           'coord type',              'geodetic';
           'coord units',             'deg';
           'map projection',          '';
           'EPSG code',               '7912'};

% convert the cell to structure
meta = cell2struct(metakey(:,2), cellfun(@(x) strrep(x, ' ', '_'), metakey(:,1), 'UniformOutput', false), 1);

% geoid data
geoid = [30.1234 31.2222 32.3456 33.4444 34.5678 36.6666;
         41.1111 42.2345 43.3333 44.4567 45.5555 46.6789;
         51.4321 52.9753 53.6543 54.8642  444444  444444;
         61.9999 62.8888 63.7777 64.6666  444444  444444];

% coordinates definition
phiGrid = (41  : -20/60 : 40)';
lamGrid = (120 :  20/60 : 121+40/60)';

% write the file
fname = isg2write(geoid, phiGrid, lamGrid, 444444, meta, 'TEST', comm);

% reload the file
[geoid_red, phiGrid_read, lamGrid_read, meta_read, comm_read] = isgread(fname);

%% ------------------------------------------------------------------------
% see Example 2 on the specifications of version 2.0
% (http://www.isgeoid.polimi.it/Geoid/format_specs.html)

% string of comments
comm = ['These data are freely available under a Creative Commons Attribution 4.0\n', ...
        'International Licence (CC BY 4.0)\n', ...
        '\n', ...
        'When using the data, please cite:\n', ...
        'A. Name1, B. Name2 (year). Title. Version 1.0. GFZ Data Services.\n', ...
        'http://doi.org/10.5880/isg.2020.001\n', ...
        '\n', ...
        'The original data were provided by C. Name3 (email of dd/mm/yyyy to ISG).\n', ...
        'The present file is distributed by ISG.\n', ...
        '\n', ...
        'This is an example.\n', ...
        'Here some information about model computation can be provided.\n', ...
        '\n', ...
        'Bibliographic reference:\n', ...
        'D. Name4, E. Name5 (year). Title. Journal, Volume(Number), pp. xxx-yyy.'];

% cell array with the header information
metakey = {'model name',              'EXAMPLE';
           'model year',              '2020';
           'model type',              'gravimetric';
           'data type',               'geoid';
           'data units',              'meters';
           'data format',             'grid';             % to be ignored by the function
           'data ordering',           'N-to-S, W-to-E';   % to be ignored by the function
           'ref ellipsoid',           'GRS80';
           'ref frame',               'ITRF2014';
           'height datum',            '';
           'tide system',             'mean-tide';
           'coord type',              'geodetic';
           'coord units',             'deg';
           'map projection',          '';
           'EPSG code',               '7912'};

% convert the cell to structure
meta = cell2struct(metakey(:,2), cellfun(@(x) strrep(x, ' ', '_'), metakey(:,1), 'UniformOutput', false), 1);

% geoid data
geoid = [30.1234 31.2222 32.3456 33.4444 34.5678 36.6666;
         41.1111 42.2345 43.3333 44.4567 45.5555 46.6789;
         51.4321 52.9753 53.6543 54.8642     nan   -9999;
         61.9999 62.8888 63.7777 64.6666     nan     nan];

% coordinates definition
phiGrid = (41  : -20/60 : 40)';
lamGrid = (120 :  20/60 : 121+40/60)';
[lam, phi] = meshgrid(lamGrid, phiGrid);

lam = lam(~isnan(geoid(:)));
phi = phi(~isnan(geoid(:)));
geoid = geoid(~isnan(geoid(:)));

% write the file
fname = isg2write(geoid, phi, lam, -9999, meta, 'TEST', comm);

% reload the file
[geoid_red, phi_read, lam_read, meta_read, comm_read] = isgread(fname);