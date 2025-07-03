function barcodes = generate_barcodes(numGenes, numBits, numOnes, writecsv)
% Generates barcodesMap to represent a specified number of genes    
%   numGenes - the number of genes you want to create barcode
%   numBits  - number of bits
%   numOnes  - number of 1s in the barcode
%   writecsv - flag to indicate 


    % create an empty array that stores the barcode
    barcodesMap = containers.Map('KeyType','char','ValueType','uint32');
    barcodes = table('Size', [numGenes,2], 'VariableTypes', ["string", "string"], 'VariableNames', ["genes", "barcodes"]);

    %for each genes, do 
    for i = 1:numGenes
        barcodeChar = generate_single_barcode(numBits, numOnes);
        while(isKey(barcodesMap, barcodeChar) == 1)
             barcodeChar = generate_single_barcode(numBits, numOnes);
        end
        barcodesMap(barcodeChar) = i; % store the barcode to map.container
        barcodes(i,:) = {i, barcodeChar};
    end
    
    if(writecsv)
        %write the barcode to csv file
        key = keys(barcodesMap);
        gene = values(barcodesMap);
        file = fopen('codebook.csv', 'wt');
        for i = 1:length(barcodesMap)
            fprintf(file, '%i, %s\n', gene{i}, key{i});
        end
        fclose(file);
    end
end


% generate the position for the 1 bits
function barcodeChar = generate_single_barcode(numBits, numOnes)
        r = randperm(numBits, numOnes);
        % generate the  barcode 
        barcode = zeros(1,numBits);
        for i = 1:numOnes
            barcode(r(i)) = 1;
        end
        barcodeChar = num2str(barcode);
%         barcodeChar = barcodeChar(find(~isspace(barcodeChar)));
end


    
    


