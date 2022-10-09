classdef MatrixSOPOT
%   MatrixSOPOT.m
%       Class for SOPOT arithmetic implementation. The matrices are
%       approximated column by column using the MPGBP SOPOT
%       approximation method. This class overloads some usual matlab 
%       operators and facilitate the SOPOT implementation.
%
%
%   Syntax: obj = MatrixSOPOT(x, spt_num, wordlength, max_pot)
%
%
%   Input:
%       . x: 2-D array to be approximated.
%       . spt_num: Total amount of active signed power of two terms.
%       . wordlength: Number of possible powers of two.
%       . max_pot: Highest possible power of two.
%
%
%   Output:
%       . Object with the following properties:
%           - Value: SOPOT approximated value.
%           - Wordlength: Number of possible powers of two.
%           - MaxNumSPT: Total amount of active signed power of two terms.
%           - MaxPot: Highest possible power of two.
%
%
%   Author: Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
%

    properties
        % ----------------------------------------------------------------
        %                   ##### PUBLIC PROPERTIES #####
        % ----------------------------------------------------------------
        Value
        RealValue
        Wordlength
        MaxNumSPT
        MaxPot
        SOPOTNum
    end
    properties (Access = private)
        % ----------------------------------------------------------------
        %                  ##### PRIVATE PROPERTIES #####
        % ----------------------------------------------------------------
        
        NumRows
        NumCols

    end
    methods
        % ----------------------------------------------------------------
        %                   ##### PUBLIC METHODS #####
        % ----------------------------------------------------------------
        function obj = MatrixSOPOT(varargin)
            % ------------------------------------------------------------
            %                       CONSTRUCTOR
            % ------------------------------------------------------------
            if nargin == 2
            % ------------------------------------------------------------
            % * PRIVATE CONSTRUCTOR:
                sopot = varargin{1};
                ipt = varargin{2};
                [row_num, col_num, word_len] = size(sopot);
                % --------------------------------------------------------
                % * PROPERTIES ALLOCATION:
                obj.RealValue = ipt.real_val;
                obj.Value = MatrixSOPOT.recover(sopot, ipt.max_pot);
                obj.Wordlength = word_len;
                obj.MaxNumSPT = ipt.max_num_spt;
                obj.MaxPot = ipt.max_pot;
                obj.NumRows = row_num;
                obj.NumCols = col_num;
                obj.SOPOTNum = sopot;
            elseif nargin == 4
            % ------------------------------------------------------------
            % * PUBLIC CONSTRUCTOR:
                val = varargin{1};
                max_num_spt = varargin{2};
                word_len = varargin{3};
                max_pot = varargin{4};
                sopot = MatrixSOPOT.mpgbp(val, max_num_spt, word_len, max_pot);
                [row_num, col_num] = size(val);
                % --------------------------------------------------------
                % * PROPERTIES ALLOCATION
                obj.RealValue = val;
                obj.Value = MatrixSOPOT.recover(sopot, max_pot);
                obj.Wordlength = word_len;
                obj.MaxNumSPT = max_num_spt;
                obj.MaxPot = max_pot;
                obj.NumRows = row_num;
                obj.NumCols = col_num;
                obj.SOPOTNum = sopot;
            end
        end
        % ----------------------------------------------------------------
        %                    *** OPERATOR OVERLOAD ***
        % ----------------------------------------------------------------
        function result = plus(o1, o2)
            % operator: +
            pre_result = [o1.SOPOTNum] + [o2.SOPOTNum];
            ipt.real_val = [o1.RealValue] + [o2.RealValue];
            ipt.max_pot = o1.MaxPot;
            ipt.max_num_spt = o1.MaxNumSPT;
            result_aux = zeros(o1.NumRows, o1.NumCols, o1.Wordlength);
            for i = 1:o1.NumRows
                for j = 1:o1.NumCols
                    result_aux(i, j, :) = MatrixSOPOT.reapprox(pre_result(i, j, :),...
                        o1.MaxNumSPT);
                end
            end
            result = MatrixSOPOT(result_aux, ipt);
        end
        function result = minus(o1, o2)
            % operator: -
            pre_result = [o1.SOPOTNum] - [o2.SOPOTNum];
            ipt.real_val = [o1.RealValue] - [o2.RealValue];
            ipt.max_pot = o1.MaxPot;
            ipt.max_num_spt = o1.MaxNumSPT;
            result_aux = zeros(o1.NumRows, o1.NumCols, o1.Wordlength);
            for i = 1:o1.NumRows
                for j = 1:o1.NumCols
                    result_aux(i, j, :) = MatrixSOPOT.reapprox(pre_result(i, j, :),...
                        o1.MaxNumSPT);
                end
            end
            result = MatrixSOPOT(result_aux, ipt);
        end
        function result = times(o1, o2)
            % operator: .*
            ipt.real_val = [o1.RealValue] .* [o2.RealValue];
            ipt.max_pot = o1.MaxPot;
            ipt.max_num_spt = o1.MaxNumSPT;
            if o1.NumRows == o2.NumRows && o1.NumCols == o2.NumCols
                % this case object 1 and 2 have same dimension
                result_aux = zeros(o2.NumRows, o2.NumCols, o2.Wordlength);
                for i = 1:o2.NumRows
                    for j = 1:o2.NumCols
                        k = MatrixSOPOT.base_product(o1.SOPOTNum(i, j, :), o2.SOPOTNum(i, j, :));
                        l = MatrixSOPOT.reapprox(k, o1.MaxPot);
                        result_aux(i, j, :) = l(o1.MaxPot+1:o1.Wordlength+o1.MaxPot);
                    end
                end
            elseif o2.NumRows == 1 && o2.NumCols == 1
                % this case one of them is a constant.
                result_aux = zeros(o1.NumRows, o1.NumCols, o1.Wordlength);
                for i = 1:o1.NumRows
                    for j = 1:o1.NumRows
                        k = MatrixSOPOT.base_product(o2.SOPOTNum, o1.SOPOTNum(i, j, :));
                        l = MatrixSOPOT.reapprox(k, o1.MaxPot);
                        result_aux(i, j, :) = l(o1.MaxPot+1:o1.Wordlength+o1.MaxPot);
                    end
                end
            else
                result_aux = zeros(o2.NumRows, o2.NumCols, o2.Wordlength);
                for i = 1:o2.NumRows
                    for j = 1:o2.NumRows
                        k = MatrixSOPOT.base_product(o1.SOPOTNum, o2.SOPOTNum(i, j, :));
                        l = MatrixSOPOT.reapprox(k, o1.MaxPot);
                        result_aux(i, j, :) = l(o1.MaxPot+1:o1.Wordlength+o1.MaxPot);
                    end
                end
            end
            result = MatrixSOPOT(result_aux, ipt);
        end
        function result = mtimes(o1, o2)
            % operator: *
            if o1.NumCols == o2.NumRows
                pre_result = zeros(o1.NumRows, o2.NumCols, o1.Wordlength);
                for k = 1:o1.NumRows
                    for i = 1:o2.NumCols
                        aux = zeros(o1.NumCols, 2*o1.Wordlength - 1);
                        for j = 1:o1.NumCols
                            aux(j, :) = MatrixSOPOT.base_product(o1.SOPOTNum(k, j, :), o2.SOPOTNum(j, i, :));
                            aux2 = sum(aux, 1);
                            aux3 = MatrixSOPOT.reapprox(aux2, o1.MaxPot);
                            pre_result(k, i, :) = aux3(o1.MaxPot+1:o1.Wordlength+o1.MaxPot);
                        end
                    end
                end
            end
            ipt.real_val = [o1.RealValue] * [o2.RealValue];
            ipt.max_pot = o1.MaxPot;
            ipt.max_num_spt = o1.MaxNumSPT;
            result = MatrixSOPOT(pre_result, ipt);
        end
        function [] = rdivide()
            % operator: ./
        end
        function [] = mrdivide()
            % operator: /
        end
        function [] = transpose()
            % operator: '
        end
    end
    methods (Static = true, Access = private)
        % ----------------------------------------------------------------
        %                  PRIVATE AND STATIC METHODS
        % ----------------------------------------------------------------
        function sopot = mpgbp(val, max_num_spt, word_len, max_pot)
            % ------------------------------------------------------------
            %               MPGBP FOR SOPOT APPROXIMATION
            % ------------------------------------------------------------
            [row_num, col_num] = size(val);
            sopot = zeros(row_num, col_num, word_len);
            if nnz(val) ~= 0
                P = floor(sqrt(row_num*col_num));
                residue = val(:);
                col_aux = zeros(row_num*col_num, word_len);
                while norm(residue) > 1e-12
                    v_rm = zeros(row_num*col_num, 1);
                    temp_val = residue;
                    [~, idx_temp] = sort(abs(temp_val), 'descend');
                    m_idx = idx_temp(1:P);
                    for i = 1:P
                        v_rm(m_idx(i)) = 1*sign(temp_val(m_idx(i)));
                    end
                    norm_v_rm = norm(v_rm);
                    v_rn = v_rm./norm_v_rm;
                    ip = (residue'*v_rn)./norm_v_rm;
                    km = ceil(log2(3/(4*ip)));
                    kfactor = 2^-km;
                    residue = residue - kfactor*v_rm;
                    if (km+max_pot+1) <= word_len && km >= -max_pot
                        for i = 1:P
                            col_aux(m_idx(i), km+max_pot+1) = col_aux(m_idx(i), km+max_pot+1) + v_rm(m_idx(i));
                        end
                    end
                end
                for row = 1:row_num
                    col_aux(row, :) = MatrixSOPOT.reapprox(col_aux(row, :), max_num_spt);
                end
                sopot = reshape(col_aux, row_num, col_num, word_len);
            end
        end
        function col_aux = reapprox(x, max_num_spt)
            % ------------------------------------------------------------
            %             SPT TERM REDUCTION FOR COLUMNS
            % ------------------------------------------------------------
            word_len = length(x);
            if ndims(x) == 3
                col_aux(:, :) = x(1, 1, :);
            else
                col_aux = x;
            end
            count = false;
            while norm(col_aux, inf) > 1
                for i = 2:word_len
                    if col_aux(i-1) > 0 && col_aux(i) < 0
                        col_aux(i-1) = col_aux(i-1) - 1;
                        col_aux(i) = col_aux(i) + 2;
                    elseif col_aux(i-1) < 0 && col_aux(i) > 0
                        col_aux(i-1) = col_aux(i-1) + 1;
                        col_aux(i) = col_aux(i) - 2;
                    end
                end
                for i = word_len:-1:2
                    col_aux(i-1) = col_aux(i-1) +...
                        sign(col_aux(i))*floor(abs(col_aux(i))/2);
                    col_aux(i) = col_aux(i) - 2*...
                        sign(col_aux(i))*floor(abs(col_aux(i))/2);
                end
                for i = 2:word_len
                    if col_aux(i-1) > 0 && col_aux(i) < 0
                        col_aux(i-1) = col_aux(i-1) - 1;
                        col_aux(i) = col_aux(i) + 2;
                    elseif col_aux(i-1) < 0 && col_aux(i) > 0
                        col_aux(i-1) = col_aux(i-1) + 1;
                        col_aux(i) = col_aux(i) - 2;
                    end
                end
                % --------------------------------------------------------
                % * OVERFLOW:
                if (abs(col_aux(1)) >= 2 && count == true) || (abs(col_aux(1)) > 2)
                    col_aux = sign(col_aux(1))*ones(1, word_len);
                elseif abs(col_aux(1)) == 2 && count == false
                    col_aux2 = sign(col_aux(1))*ones(1, word_len);
                    col_aux(1) = col_aux(1) - sign(col_aux(1));
                    col_aux = col_aux + col_aux2;
                    count = true;
                end
            end
            % ------------------------------------------------------------
            % * ROUNDOFF:
            nz = nnz(col_aux);
            while nz > max_num_spt
                for i = word_len:-1:1
                    if col_aux(i) ~= 0
                        col_aux(i) = 0;
                        nz = nnz(col_aux);
                        break
                    end
                end
            end
        end
        function result = base_product(x, y)
            % ------------------------------------------------------------
            %                BASIC SOPOT PRODUCT
            % ------------------------------------------------------------
            word_len = length(x);
            result = zeros(word_len, 2*word_len-1);
            for i = 1:word_len
                for j = 1:word_len
                    result(i, i+j-1) = x(i)*y(j);
                end
            end
            result = sum(result, 1);
        end
        function result = recover(x, max_pot)
            % ------------------------------------------------------------
            %           RECOVER FROM SOPOT REPRESENTATION
            % ------------------------------------------------------------
            [num_row, num_col, word_len] = size(x);
            m = zeros(num_row*num_col*word_len, 1);
            bit_plane = flip(ceil(find(flip(x, 3)~=0)./(num_row*num_col))+(max_pot-word_len));
            n = x(x~=0);
            m(x~=0) = n(:).*2.^bit_plane;
            aux = reshape(m, num_row, num_col, word_len);
            result = sum(aux, 3);
        end
    end
end
