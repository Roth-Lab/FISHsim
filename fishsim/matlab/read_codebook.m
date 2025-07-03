
function barcodes = read_codebook(filename)
    codebookI = readtable(filename, 'Format', '%s %s %s'); 
    codebook = codebookI{:, {'name', 'barcode'}};
    name2 = convertCharsToStrings(codebook(:,1));
    barcode = convertCharsToStrings(codebook(:,2)); 
    barcodes = table('Size', [size(barcode,1),2], 'VariableTypes', ["string", "string"], 'VariableNames', ["genes", "barcodes"]);
    barcodes.genes = name2;
    barcodes.barcodes = barcode;
end