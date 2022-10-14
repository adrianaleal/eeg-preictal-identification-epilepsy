# Unsupervised EEG Preictal Interval Identification in Patients with Drug-resistant Epilepsy

__Code to reproduce results reported in our paper published as:__

Leal, A., Curty, J., Lopes, F., Pinto, M.F., Oliveira, A., Sales, F., Ruano, M.G., Dourado, A., Henriques, J., and Teixeira, C.A. "Unsupervised EEG Preictal Interval Identification in Patients with Drug-resistant Epilepsy." (2022)

__How to run the code__

(1) For feature extraction run main_feature_extraction.m script and read the corresponding description.

(2) For feature reduction run main_feature_reduction.m script and read the corresponding description.

(3) For clustering run main_clustering.m script and read the corresponding description.


__Folder named "FiguresReducedData" contains figures depicting reduced data and clustering solution before each seizure's onset (similar to Figure 4 in the paper)__

__Excel file named "results_eeg_preictal_unsupervised_learning.xlsx" contains the results for:__ 

- Feature reduction (including final parameters obtained for t-SNE and UMAP).
- Clustering (including the final clustering methods and the Dunn's Index value). 
- Preictal interval inspection (including the starting time and duration of the preictal interval and the correlation with the sleep-wake cycle)
