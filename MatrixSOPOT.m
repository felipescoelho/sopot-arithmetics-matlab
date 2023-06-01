classdef MatrixSOPOT
    % MatrixSOPOT       SOPOT numeric object
    %
    %   Syntax:
    %       xSOPOT = MatrixSOPOT(x, maxNumSPT, wordLength, maxPower)
    %
    %   Description:
    %       MatrixSOPOT(x, maxNumSPT, wordLength, maxPower) returns a SOPOT
    %       object with approximated value x, total amount of active
    %       signed-power-of-two terms equals to maxNumSPT, the word length
    %       in the fixed-point-like implementation is defined by
    %       wordLength, and maxPower defines the highest possible power of
    %       two in the SOPOT.
    %       This code enables a fixed-point-like implementation of the
    %       SOPOT arithmetics in MATLAB.
    %
    %   A MatrixSOPOT object has the following properties:
    %
    %       Value       - SOPOT value, converted to MATLAB's double
    %       WordLength  - Word length of the object, stored as an integer
    %       MaxNumSPT   - Total amount of active signed-power-of-two terms,
    %                     stored as an integer
    %       MaxPower    - Highest possible power in the SOPOT
    %
    %   Author:
    %       Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
    %
    %   P.S.:
    %       If you use this code, please cite our works on the subject.
    %   
    
    properties
        Value
        WordLength
        MaxNumSPT
        MaxPower
    end
    properties (Access = private)
        RealPart
        ImagPart
        NumRows
        NumCols
    end
    methods
        function obj = MatrixSOPOT(varargin)
            % MatrixSOPOT       Object constructor
            %
            %   Syntax:
            %       obj = MatrixSOPOT(x, maxNumSPT, wordLength, maxPower)
            %       obj = MatrixSOPOT(realSOPOTArray, imaginarySOPOTArray,
            %           dataInfo)
            %
            %   Description:
            %       MatrixSOPOT(x, maxNumSPT, wordLength, maxPower) returns
            %       an object with approximated value x, total amount of
            %       active signed-power-of-two terms eqaults to maxNumSPT,
            %       the word length in the fixed-point-like implementation
            %       is defined by wordLength, and maxPower defines the
            %       highest possible power of two in the SOPOT.
            %
            %       MatrixSOPOT(realSOPOTArray, imaginarySOPOTArray,
            %       dataInfo) works as a private constructor, defining a
            %       new object with real part equals to realSOPOTArray,
            %       imaginary part equals to imaginarySOPOTArray, and
            %       informations on the implementation in dataInfo.
            
            if nargin == 3  % Private constructor
                realPart = varargin{1};
                imagPart = varargin{2};
                dataInfo = varargin{3};
                maxPower = dataInfo.MaxPower;
                wordLength = dataInfo.WordLength;
                maxNumSPT = dataInfo.MaxNumSPT;
                [numRows, numColumns, ~] = size(realPart);
                if nnz(imagPart) ~= 0
                    obj.Value = MatrixSOPOT.recover(realPart, maxPower) ...
                        + 1j*MatrixSOPOT.recover(imagPart, maxPower);
                else
                    obj.Value = MatrixSOPOT.recover(realPart, maxPower);
                end
                obj.RealPart = realPart;
                obj.ImagPart = imagPart;
                obj.WordLength = wordLength;
                obj.MaxNumSPT = maxNumSPT;
                obj.MaxPower = maxPower;
                obj.NumRows = numRows;
                obj.NumCols = numColumns;
            elseif nargin == 4  % Public Constructor
                inputArray = varargin{1};
                maxNumSPT = varargin{2};
                wordLength = varargin{3};
                maxPower = varargin{4};
                [numRows, numColumns] = size(inputArray);
                realPart = MatrixSOPOT.mpgbp(real(inputArray), ...
                    maxNumSPT, wordLength, maxPower);
                if ~isreal(inputArray)
                    imagPart = MatrixSOPOT.mpgbp( ...
                        imag(inputArray), maxNumSPT, wordLength, maxPower);
                else
                    imagPart = zeros(numRows, numColumns, wordLength);
                end
                if nnz(imagPart) ~= 0
                    obj.Value = MatrixSOPOT.recover(realPart, ...
                        maxPower) + 1j*MatrixSOPOT.recover( ...
                        imagPart, maxPower);
                else
                    obj.Value = MatrixSOPOT.recover(realPart, maxPower);
                end
                obj.RealPart = realPart;
                obj.ImagPart = imagPart;
                obj.WordLength = wordLength;
                obj.MaxNumSPT = maxNumSPT;
                obj.MaxPower = maxPower;
                obj.NumRows = numRows;
                obj.NumCols = numColumns;
            end
        end
        
        function varargout = size(obj)
            % SIZE  Array size
            %   
            %   sz = size(A)
            %
            %   Description:
            %       size(obj) returns a row vector whose elements are the
            %       lengths of the corresponding dimensions of the SOPOT
            %       array (obj).
            
            arraySize = [obj.NumRows obj.NumCols];
            if nargout == 0 || nargout == 1
                varargout{1} = arraySize;
            elseif nargout == 2
                varargout{1} = arraySize(1);
                varargout{2} = arraySize(2);
            end
        end
        
        function arrayLength = length(obj)
            % length    Length of largest array dimension.
            %
            %   Syntax:
            %       arrayLength = length(obj)
            %
            %   Description:
            %       length(obj) returns the length of the largest dimension
            %       of the SOPOT array (obj).
            
            arrayLength = max(obj.NumRows, obj.NumCols);
        end
        
        function result = isreal(objA)
            % isreal        Determine whether an array his complex
            %
            %   Syntax:
            %       result = isreal(objA)
            %
            %   Description:
            %       result = isreal(objA) returns logical 1 (true) when
            %       SOPOT array does not have an imaginary part, and
            %       logical 0 (false) otherwise. ISREAL returns logical 1
            %       (true) for any SOPOT that have zero imaginary part.
            
            if objA.ImagPart == 0
                result = true;
            else
                result = false;
            end
        end
        
        function resultObj = uplus(objA)
            % uplus, +      Unary plus
            %
            %   Syntax:
            %       resultObj = +objA
            %       resultObj = uplus(objA)
            %
            %   Description:
            %       resultObj = +objA returns array objA and stores it in
            %       resultObj.
            
            resultObj = objA;
        end
        
        function resultObj = uminus(objA)
            % uminus, -     Unary minus
            %
            %   Syntax:
            %       resultObj = -objA
            %       resultObj = uminus(objA)
            %
            %   Description:
            %       resultObj = -objA negates the values of the elements in
            %       objA and stores the result in resultObj.
            
            realPart = - objA.RealPart;
            imagPart = - objA.ImagPart;
            dataInfo.WordLength = objA.WordLength;
            dataInfo.MaxNumSPT = objA.MaxNumSPT;
            dataInfo.MaxPower = objA.MaxPower;
            resultObj = MatrixSOPOT(realPart, imagPart, dataInfo);
        end

        function resultObj = plus(objA, objB)
            % plus, +       Add SOPOT numbers
            %
            %   Syntax:
            %       resultObj = objA + objB
            %       resultObj = plus(objA, objB)
            %
            %   Description:
            %       resultObj = objA + objB returns a new object from the
            %       addition of the values in objA and objB by adding
            %       corresponding elements.
            %
            %       The sizes of objA and objB must be the same or be
            %       compatible. If the sizes of objA and objB are
            %       compatible, then the two arrays implicitly expand to
            %       match each other.
            
            complexTest = ~isreal(objA) || ~isreal(objB);
            if isequal(size(objA), size(objB))  % Same size
                [numRows, numCols] = size(objA);
                resultReal = zeros(numRows, numCols, objA.WordLength);
                resultImag = zeros(numRows, numCols, objA.WordLength);
                if complexTest
                    for row = 1:numRows
                        for col = 1:numCols
                            resultReal(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                objA.RealPart(row, col, :), ...
                                objB.RealPart(row, col, :), ...
                                objA.MaxPower, objA.WordLength);
                            resultImag(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                objA.ImagPart(row, col, :), ...
                                objB.ImagPart(row, col, :), ...
                                objA.MaxPower, objA.WordLength);
                        end
                    end
                else
                    for row = 1:numRows
                        for col = 1:numCols
                            resultReal(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                objA.RealPart(row, col, :), ...
                                objB.RealPart(row, col, :), ...
                                objA.MaxPower, objA.WordLength);
                        end
                    end
                end
            elseif isscalar(objA)  % objA is scalar
                [numRows, numCols] = size(objB);
                resultReal = zeros(numRows, numCols, objA.WordLength);
                resultImag = zeros(numRows, numCols, objA.WordLength);
                if complexTest
                    for row = 1:numRows
                        for col = 1:numCols
                            resultReal(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                objA.RealPart(1, 1, :), ...
                                objB.RealPart(row, col, :), ...
                                objA.MaxPower, objA.WordLength);
                            resultImag(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                objA.ImagPart(1, 1, :), ...
                                objB.ImagPart(row, col, :), ...
                                objA.MaxPower, objA.WordLength);
                        end
                    end
                else
                    for row = 1:numRows
                        for col = 1:numCols
                            resultReal(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                objA.RealPart(1, 1, :), ...
                                objB.RealPart(row, col, :), ...
                                objA.MaxPower, objA.WordLength);
                        end
                    end
                end
            elseif isscalar(objB)  % objB is scalar
                [numRows, numCols] = size(objA);
                resultReal = zeros(numRows, numCols, objA.WordLength);
                resultImag = zeros(numRows, numCols, objA.WordLength);
                if complexTest
                    for row = 1:numRows
                        for col = 1:numCols
                            resultReal(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                objA.RealPart(row, col, :), ...
                                objB.RealPart(1, 1, :), ...
                                objA.MaxPower, objA.WordLength);
                            resultImag(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                objA.ImagPart(row, col, :), ...
                                objB.ImagPart(1, 1, :), ...
                                objA.MaxPower, objA.WordLength);
                        end
                    end
                else
                    for row = 1:numRows
                        for col = 1:numCols
                            resultReal(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                objA.RealPart(row, col, :), ...
                                objB.RealPart(1, 1, :), ...
                                objA.MaxPower, objA.WordLength);
                        end
                    end
                end
            elseif ismatrix(objA) && iscolumn(objB)
                % ObjA is matrix and objB is column
                [numRowsA, numColsA] = size(objA);
                [numRowsB, ~] = size(objB);
                if numRowsA == numRowsB
                    resultReal = zeros(numRowsA, numColsA, objA.WordLength);
                    resultImag = zeros(numRowsA, numColsA, objA.WordLength);
                    if complexTest
                        for row = 1:numRowsA
                            for col = 1:numColsA
                                resultReal(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                    objA.RealPart(row, col, :), ...
                                    objB.RealPart(row, 1, :), ...
                                    objA.MaxPower, objA.WordLength);
                                resultImag(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                    objA.ImagPart(row, col, :), ...
                                    objB.ImagPart(row, 1, :), ...
                                    objA.MaxPower, objB.WordLength);
                            end
                        end
                    else
                        for row = 1:numRowsA
                            for col = 1:numColsA
                                resultReal(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                    objA.RealPart(row, col, :), ...
                                    objB.RealPart(row, 1, :), ...
                                    objA.MaxPower, objA.WordLength);
                            end
                        end
                    end
                else
                    message = 'Arrays have incompatible sizes for this' ...
                        + ' operation';
                    error(message)
                end
            elseif iscolumn(objA) && ismatrix(objB)
                % ObjA is column and objB is matrix
                [numRowsA, ~] = size(objA);
                [numRowsB, numColsB] = size(objB);
                if numRowsA == numRowsB
                    resultReal = zeros(numRowsB, numColsB, objA.WordLength);
                    resultImag = zeros(numRowsB, numColsB, objA.WordLength);
                    if complexTest
                        for row = 1:numRowsB
                            for col = 1:numColsB
                                resultReal(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                    objB.RealPart(row, col, :), ...
                                    objA.RealPart(row, 1, :), ...
                                    objA.MaxPower, objA.WordLength);
                                resultImag(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                    objB.ImagPart(row, col, :), ...
                                    objA.ImagPart(row, 1, :), ...
                                    objA.MaxPower, objB.WordLength);
                            end
                        end
                    else
                        for row = 1:numRowsB
                            for col = 1:numColsB
                                resultReal(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                    objB.RealPart(row, col, :), ...
                                    objA.RealPart(row, 1, :), ...
                                    objA.MaxPower, objA.WordLength);
                            end
                        end
                    end
                else
                    message = 'Arrays have incompatible sizes for this' ...
                        + ' operation';
                    error(message)
                end
            elseif iscolumn(objA) && isrow(objB)
                % objA is column and objB is row
                numRows = length(objB);
                numCols = length(objA);
                resultReal = zeros(numRows, numCols, objA.WordLength);
                resultImag = zeros(numRows, numCols, objA.WordLength);
                if complexTest
                    for row = 1:numRows
                        for col = 1:numCols
                            resultReal(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                objB.RealPart(1, col, :), ...
                                objA.RealPart(row, 1, :), ...
                                objA.MaxPower, objA.WordLength);
                            resultImag(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                objB.ImagPart(1, col, :), ...
                                objA.ImagPart(row, 1, :), ...
                                objA.MaxPower, objB.WordLength);
                        end
                    end
                else
                    for row = 1:numRows
                        for col = 1:numCols
                            resultReal(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                objB.RealPart(1, col, :), ...
                                objA.RealPart(row, 1, :), ...
                                objA.MaxPower, objA.WordLength);
                        end
                    end
                end
            elseif isrow(objA) && iscol(objB)
                % objA is row and objB is column
                numRows = length(objA);
                numCols = length(objB);
                resultReal = zeros(numRows, numCols, objA.WordLength);
                resultImag = zeros(numRows, numCols, objA.WordLength);
                if complexTest
                    for row = 1:numRows
                        for col = 1:numCols
                            resultReal(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                objA.RealPart(1, col, :), ...
                                objB.RealPart(row, 1, :), ...
                                objA.MaxPower, objA.WordLength);
                            resultImag(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                objA.ImagPart(1, col, :), ...
                                objB.ImagPart(row, 1, :), ...
                                objA.MaxPower, objB.WordLength);
                        end
                    end
                else
                    for row = 1:numRows
                        for col = 1:numCols
                            resultReal(row, col, :) = MatrixSOPOT.sum_sopots( ...
                                objA.RealPart(1, col, :), ...
                                objB.RealPart(row, 1, :), ...
                                objA.MaxPower, objA.WordLength);
                        end
                    end
                end
            else
                error('Arrays have incompatible sizes for this operation.')
            end
            dataInfo.WordLength = objA.WordLength;
            dataInfo.MaxPower = objA.MaxPower;
            dataInfo.MaxNumSPT = objA.MaxNumSPT;
            resultObj = MatrixSOPOT(resultReal, resultImag, dataInfo);
        end
        
        function result = minus(o1, o2)
            % operator: -
            
            dummy = -o2;
            
            result = o1+dummy;
        end
        
        function result = times(o1, o2)  % TODO
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
        
        function result = mtimes(o1, o2)  % TODO
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
        
        function [] = rdivide()  % TODO
            % operator: ./
        end
        
        function [] = mrdivide()  % TODO
            % operator: /
        end
        
        function [] = transpose()  % TODO
            % operator: '
        end
    end
    
    methods (Static = true, Access = private)
        
        function resultSOPOT = sum_sopots(SOPOTA, SOPOTB, maxNumSPT, ...
                wordLength)
            % Sum two SOPOTs
            %
            % Returns a SOPOT element containing an already reduced form of
            % SOPOT.

            resultSOPOT = zeros(size(SOPOTA));
            auxVar = 0;
            for indexer = wordLength:-1:1
                auxVar = auxVar + SOPOTA(:, :, indexer) ...
                    + SOPOTB(:, :, indexer);
                if abs(auxVar) ~= 2
                    resultSOPOT(:, :, indexer) = auxVar;
                    auxVar = 0;
                elseif indexer > 1
                    auxVar = sign(auxVar);
                else
                    resultSOPOT(1, 1, :) = reshape([ ...
                        sign(auxVar)*ones(1, maxNumSPT) ...
                        zeros(1, wordLength-maxNumSPT)], size(SOPOTA));
                end
            end
            indexCount = wordLength;
            while nnz(resultSOPOT) > maxNumSPT
                resultSOPOT(:, :, indexCount) = 0;
                indexCount = indexCount - 1;
            end
        end
        
        function sopot = mpgbp(inputArray, maxNumberSPT, wordLength, maxPower)
            % mpgbp
            [numberRows, numberColumns] = size(inputArray);
            sopot = zeros(numberRows, numberColumns, wordLength);
            if nnz(inputArray) ~= 0
                P = floor(sqrt(numberRows*numberColumns));
                residue = inputArray(:);
                col_aux = zeros(numberRows*numberColumns, wordLength);
                while norm(residue) > 1e-12
                    v_rm = zeros(numberRows*numberColumns, 1);
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
                    if (km+maxPower+1) <= wordLength && km >= -maxPower
                        for i = 1:P
                            col_aux(m_idx(i), km+maxPower+1) = col_aux(m_idx(i), km+maxPower+1) + v_rm(m_idx(i));
                        end
                    end
                end
                for row = 1:numberRows
                    col_aux(row, :) = MatrixSOPOT.reapprox(col_aux(row, :), maxNumberSPT);
                end
                sopot = reshape(col_aux, numberRows, numberColumns, wordLength);
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
        function result = recover(x, maxPower)
            % ------------------------------------------------------------
            %           RECOVER FROM SOPOT REPRESENTATION
            % ------------------------------------------------------------
            [num_row, num_col, word_len] = size(x);
            m = zeros(num_row*num_col*word_len, 1);
            bit_plane = flip(ceil(find(flip(x, 3)~=0)./(num_row*num_col))+(maxPower-word_len));
            n = x(x~=0);
            m(x~=0) = n(:).*2.^bit_plane;
            aux = reshape(m, num_row, num_col, word_len);
            result = sum(aux, 3);
        end
    end
end


% EoF
