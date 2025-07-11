# Depth Functions for Object Data

This project contains code to reproduce the results of the paper named 'Metric Oja Depth, New Statistical Tool for
Estimating the Most Central Objects' by Zamanifarizhandi, V. and Virta, J.

The goal of the project was to develope a novel measure of statistical depth, the metric Oja depth applicable to any object data. 

All code here is written in R and the main used packages are CovTools, ICtest, parallel, fda, ggplot2 and GA.

## The reason why this project is useful for you
If you're working with complex datasets, such as matrices or images, the following main functions can help you to identify the most central objects within your data. Simply calculate the distances between objects and apply them as an input in our metric functions.

### The main metric functions are:

MLD_cpp, MHDVJ_cpp, MSD_cpp, MOD2_cpp, MOD3_cpp,

coded in Rcpp. To use them you need to install from source the R-package metricDepth in the subfolder "All Metric depth function in C++".

In addition the project contains several helper functions for data generation, computing p-value for real dataset examples and the simulations performed in the paper.

### Great! Then How to Use this Project!

After calling all necessary libraries, you need to fetch 'All Metric depth function' file by following line in your R script:

`source("YOUR/PATH/All Metric depth functions.R")`

Next step is to import your dataset and compute the distance between objects in it. For example if you have functional data use metric.lp(YOUR FD DATA) to acheive this. Finally based on our suggestions in paper choose which depth functions meets your needs the best and insert metric.lp(YOUR FD DATA) as an input in chosen depth function. Now you have found the Most central Object in your dataset!

You can see the more detailed examples in [**Real Data example**](https://github.com/vidazamani/Depth-functions-for-Object-Data/tree/15c63bb935e6261f982a280901ab1a48eae834b9/Real%20data%20examples) Folder in the current repository.


### Authors

Zamanifarizhandi, V. and Virta, J.


### References
1. Tukey, J.W., 1975. Mathematics and the picturing of data, in: James,
R. (Ed.), Proceedings of the International Congress of Mathematicians,
Canadian Mathematical Congress. pp. 523–531.

2. Dai, X., Lopez-Pintado, S., 2023. Tukey’s depth for object data. Journal
of the American Statistical Association 118, 1760–1772.

3. Virta, J., 2023b. Spatial Depth for Data in Metric Spaces. arXiv preprint
arXiv:2306.09740.

4. Cholaquidis, A., Fraiman, R., Gamboa, F., Moreno, L., 2023. Weighted
Lens Depth: Some Applications to Supervised Classification. Canadian
Journal of Statistics 51, 652–673.

5. Geenens, G., Nieto-Reyes, A., Francisci, G., 2023. Statistical Depth in
Abstract Metric Spaces. Statistics and Computing 33, 46.
