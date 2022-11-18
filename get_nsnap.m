function nsnap = get_nsnap(RunName,RunNbr,nbp,path_to_data)
% Before using this function make sure to create the txt file containing the total number of line in each fra file by using the command: wc -l FILENAME >> OUTPUT.TXT 

	nsnap = load( [ path_to_data RunName '_' RunNbr '.txt' ] );
	nsnap = nsnap(:,1)/(2*nbp);

end

