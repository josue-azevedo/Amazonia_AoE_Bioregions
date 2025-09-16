This repository contains R scripts used in the analyses for the paper:

Manifold, megadiverse and multiscaled: Endemism and Regionalization patterns of Squamate Reptiles in Amazonia

The scripts are organized to perform complementary steps of the analyses, including calculation of areas of endemism (Biotic Elements), bioregionalization, model selection, and clade age estimation.

⸻

Repository structure:

	•	BE_End_All_github.r
Script for the calculation of biotic elements (part of Areas of Endemism), based on species occurrence data in the form of polygons.

	•	Rcluster_Am_regionalization_All_species_V_github.R
Script for bioregionalization analyses, using Rcluster. It generates dendrograms of spatial relationships among regions and plots the results on maps.

	•	model_selecion_varImp_github.R
Script for environmental variable processing and model selection (Multinomial and GLS models).

	•	AoE_clade_age_github.R
Script for the calculation of Clade Age metrics, linking phylogenetic information with spatial patterns.
