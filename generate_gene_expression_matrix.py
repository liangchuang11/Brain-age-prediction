import abagen

data_dir = 'file_path_Allen_data'
atlas_file = 'augmented Schaefer-1016.nii'
expression = abagen.get_expression_data(atlas_file, missing='centroids',data_dir=data_dir)
expression.to_csv('gene_expression_matrix.csv')
