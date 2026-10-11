=================================
CTX-HPF SmartSeq sequencing data
=================================

RNA sequencing data of single cells isolated from mouse cortex and hippocampal formation.
The data set was described in Yao et al., Cell 2021

Gene expression data matrices
	expression_matrix_SSv4.hdf5
		Contains the (row, column) matrix of read counts obtained for each (gene, nucleus) based on 

Sample information (metadata_SSv4.csv)
	sample_name
		Unique sample identifier
	donor_sex_label
		Sex of donor
	region_label
		Region from which sample was dissected
	platform_label
		Single cell sequencing methods
	cluster_label
		Cell type cluster name (Level4)
	subclass_label
		Cell type subclass name (Level3)
	neighborhood_label
		Cell type neighborhood name (Level2)
	class_label
		Broad cell type class (e.g. GABAergic, Glutamatergic or Non-neuronal) (Level1)
	external_donor_name_label
		Unique donor identifier
	full_genotype_label
		Genotype of donor
	facs_population_plan_label
		Gating criteria used to select cells for sorting
	injection_roi_label
		If injected, region of brain in which injection took place
	injection_materials_label
		If injected, material that was injected
	injection_method_label
		If injected, method used for injection
	cell_type_accession_label
		Unique cell type identified according to cell type nomenclature (https://portal.brain-map.org/explore/classes/nomenclature)
	cell_type_alias_label
		Cluster label like above
	cell_type_alt_alias_label
		Alternative cell type label if type is known under different alias.
	cell_type_designation_label
	cortical_layer_label
		Cortical layer, if dissected from specific cortical layer
	expected_i1_fastq
		expected fastq file name

Sample group mapping
    sample_group_mapping.csv
        Contains the mapping from each GEO Sample record to multiple sample_names.



