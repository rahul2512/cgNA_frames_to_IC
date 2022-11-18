function uvwp = mergeuvwp(runName, nfile)

tmp = load([ runName '_uvwp_nj_5scale_1.mat' ]) ;
names = fieldnames(tmp) ;
nfields = length(names) ;

uvwp = struct(tmp) ;

for i = 2:nfile

tmp = load([ runName '_uvwp_nj_5scale_' num2str(i) '.mat' ]) ;

for j = 1:nfields
  
  uvwp.(names{j}) = [uvwp.(names{j}) ; tmp.(names{j})] ;
  
end

end

save([ runName '_uvwp_nj_5scale.mat' ], '-struct', 'uvwp' ) ;

end

