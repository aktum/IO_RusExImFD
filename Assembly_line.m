%% Data import
% 43 countries+rest of the world, 56 sectors, 44*56=2464
wiot = importdata('WIOD_data/WIOT2006_Nov16_ROW.csv');
%% Intermediate supply-use table
Z = wiot.data(1:2464, 1:2464);
%% Value added
VA = wiot.data(2470,1:2464);
%% Output
X = wiot.data(2472,1:2464);
%% Final demand
FD = wiot.data(1:2464,2465:2684);
%% taxes less subsidies on intermediate products
TLS_intermediate = wiot.data(2467,1:2464);
%% taxes less subsidies on final products
TLS_final = wiot.data(2467,2465:2684);
%% Structural Indicators - based on the values presented in ICIO

% -----------------------------------
% A - Transaction Matrix
% W - Value Added row vector
% X - Output row vector
% -----------------------------------

%% Value added
W = TLS_intermediate + VA; 
%% Value added as a share of Gross Output, by industry, percentage
V = W./X;
V(isnan(V)) = 0 ;
%% Leontief inverse
A = Z./X;
X_2464 = repmat(X,2464,1);
A2 = Z./X_2464;
A(isnan(A)) = 0;
t = isequal(A,A2);
B = inv(eye(2464)-A);

%% Leontief with off-diagonal block elements equal to zero
DiagonalZeroBlocks = kron(diag(diag(ones(44))), ones(56,56));
B_offdiag_c_zero = DiagonalZeroBlocks.*B;
clear DiagonalZeroBlocks

%% Leontief with diagonal block elements equal to zero
DiagonalZeroBlocks = kron(ones(44)-diag(diag(ones(44))), ones(56,56));
B_diag_c_zero = DiagonalZeroBlocks.*B;
clear DiagonalZeroBlocks

%% GRTR_INT, EXGR_INT
tic
% diag(diag(A)) returns a diagonal matrix composed of the diagonal
% elements of the original matrix A.
% kron returns the Kronecker tensor product of matrices A and B.
% If A is an m-by-n matrix and B is a p-by-q matrix, then kron(A,B)
%is an m*p-by-n*q matrix formed by taking all possible products between
% the elements of A and the matrix B.
DiagonalZeroBlocks = kron(ones(44)-diag(diag(ones(44))), ones(56,56));
% element-wise multiplication
GRTR_INT = DiagonalZeroBlocks.*Z;
% Summing every 36 column for each country
EXGR_INT = reshape(sum (reshape(GRTR_INT',56, [])), [], ...
    size(GRTR_INT, 1));
EXGR_INT = (EXGR_INT)';
clear DiagonalZeroBlocks
toc

%% GRTR_FNL, EXGR_FNL
tic
% diag(diag(A)) returns a diagonal matrix composed of the diagonal
% elements of the original matrix A.
% kron returns the Kronecker tensor product of matrices A and B.
% If A is an m-by-n matrix and B is a p-by-q matrix, then kron(A,B)
%is an m*p-by-n*q matrix formed by taking all possible products between
% the elements of A and the matrix B.
DiagonalZeroBlocks = kron(ones(44)-diag(diag(ones(44))), ones(56,5));
DiagonalZeroBlocks(:,221:end) = [];
% element-wise multiplication
GRTR_FNL = DiagonalZeroBlocks.*FD;
% Summing every 36 column for each country
EXGR_FNL = reshape(sum (reshape(GRTR_FNL',5, [])), [], ...
    size(GRTR_FNL', 2));
EXGR_FNL = (EXGR_FNL)';
clear DiagonalZeroBlocks
toc

% -----------------------------------
%% Country-industry exports to world
% -----------------------------------
EXGR = sum(EXGR_INT,2)+sum(EXGR_FNL,2);

% -----------------------------------
%% Intermediate Imports
% -----------------------------------

% Imports of intermediates by country (column) from country-industry
IMGR_INT = EXGR_INT;
% Imports of final goods by country (column) from country-industry
IMGR_FNL = EXGR_FNL;
% Total imports by country (column) from country-industry
IMGR = IMGR_INT + IMGR_FNL;

%% Gross trade balance, by partner country
% Exports by country to other countries
% Sum every 56 rows of IMGR matrix
EXGR_INT_cp = reshape(sum (reshape(IMGR_INT,56, [])), [], ...
    size(IMGR_INT, 2));
EXGR_FNL_cp = reshape(sum (reshape(IMGR_FNL,56, [])), [], ...
    size(IMGR_FNL, 2));
EXGR_cp = EXGR_INT_cp + EXGR_FNL_cp;
IMGR_cp = EXGR_cp';
BALGR = EXGR_cp - IMGR_cp;

%% EXGRpSH: Gross exports, partner shares, by industry, percentage

EXGRpSH = IMGR./EXGR.*100;

%% EXGR_DVA: Domestic value added content of gross exports, USD million

EXGR_DVA_ciworld = (V*B_offdiag_c_zero*diag(EXGR))';

%% EXGR_DVASH: Domestic value added share of gross exports, percentage

EXGR_DVASH_ci = EXGR_DVA_ciworld./EXGR.*100;

%% EXGR_DDC: Direct domestic industry value added content of gross
% exports, USD million

% Local Leontief inverse
% Leontief with off-diagonal block elements equal to zero
DiagonalZeroBlocks = kron(diag(diag(ones(44))), ones(56,56));
B_local = inv((eye(2464)-A).*DiagonalZeroBlocks);
clear DiagonalZeroBlocks

EXGR_DDC_c = diag(V)*diag(diag(B_local))*EXGR;
