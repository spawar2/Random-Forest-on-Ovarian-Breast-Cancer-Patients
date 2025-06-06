# Random-Forest-on-Ovarian-Breast-Cancer-Patients. Date created/updated: December, 9, 2024.
[R: GEOquery, Biobase, preprocessCore, multiClust].
<img width="290" alt="AMMS" src="https://github.com/spawar2/Random-Forest-on-Ovarian-Breast-Cancer-Patients/assets/25118302/9e3af9bd-0d0d-440f-91de-8bf982217fa8">

[Keynote Speaker, 4th International Conference on Applied Mathematics and Simulation, 2 February, 2021. Machine learning application in genomics, by Shrikant Pawar, Keynote speaker and general chair,](http://www.4th-amms.org/com.html) Github, [*2021: 13],[**3].
Video presentation: https://campuspress.yale.edu/shrikantpawar/files/2024/05/6th-world-conference.pptx
https://www.youtube.com/watch?v=Y6skvhHVR2w&ab_channel=ShrikantPawar
Claflin University, Orangeburg, South Carolina, United States of America (USA). 
https://www.claflin.edu/
https://www.claflin.edu/academics-research/schools-departments/school-of-natural-sciences-and-mathematics/department-of-mathematics-computer-science/computer-science

Common cancer biomarkers of breast and ovarian types identified through artificial intelligence, Publication: Wileys: Chemical Biology and Drug Design, Publication date: May 15, 2020, IF=2.9, Shrikant Pawar, Tuck Onn Liew, Stanam, A., Dr. Lahiri, collaboration with Dr. Chandrajit. Lahiri, Sunway University, Malaysia .(Cited by 10)^^^ DOI: 10.1111/cbdd.13672, Volume: 96, Issue: 3, Pages: 995-1004.††
<img width="686" alt="Screenshot 2025-05-11 at 2 05 54 PM" src="https://github.com/user-attachments/assets/9e537e91-c690-448f-b738-6c5f6f4dc715" />

https://campuspress.yale.edu/shrikantpawar/files/2025/05/Common_cancer_biomarkers_of_breast_and_o.pdf

README.md: Breast, Ovarian, Colon, Lung cancer Microarray data read, quantile  Normalization, data Test-Train Split, Neural, cluster_analysis function for KMEANS ANALYSIS, Hierarchial, Random Forest, confusion matrix, accuracy, sensitivity, specificity, precision, recall, confusion matrix, log-loss, and area under curve and receiver operating characteristic, AUC-ROC evaluation.
selected function(getGEO, normalize.quantiles, merge, cluster_analysis, hclust, Kmeans, mas5, rowMeans, randomForest, survfit, chisq.test, pData, rep, colnames, factor, eBayes, decideTests, topTable, read.tree, plot, str, write.tree, library, setwd, ReadAffy, exprs, read.csv, read.delim, write.table, roundPhylogram, unroot, str, write.tree, RMA, read.table).

Testing: table(testing$V2,pred_test) Prediction_test alive dead alive 214 5 dead 31 11 ((214+11)/(nrow(testing)))*100 [1] 86.2069.
