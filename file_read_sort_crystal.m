clear all, close all

lammpsdata = importdata('final_data_crystallization_0_.data');

natoms = 10000;
start = 38;
rowstart = start;
rowend = natoms+start-1;

indexmatrix = zeros(natoms,4);
for i = 1:natoms
        indexmatrix((i),1) = str2double(cell2mat(lammpsdata.textdata([i+rowstart-1],1)));
        indexmatrix((i),2) = str2double(cell2mat(lammpsdata.textdata([i+rowstart-1],5)));
        indexmatrix((i),3) = str2double(cell2mat(lammpsdata.textdata([i+rowstart-1],6)));
        indexmatrix((i),4) = str2double(cell2mat(lammpsdata.textdata([i+rowstart-1],7)));
end

sortindexmatrix = sortrows(indexmatrix);
positionmatrix = [];
positionmatrix(:,:) = sortindexmatrix(:,[2:4]);




