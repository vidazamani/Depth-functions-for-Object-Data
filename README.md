# Depth Functions for Object Data

This project contains code to reproduce the results of the paper named 'Metric Oja Depth, New Statistical Tool for
Estimating the Most Central Objects' by Zamanifarizhandi, V. and Virta, J.

The goal of the project was to develop a novel measure of statistical depth, the metric Oja depth applicable to any object data. 

All code here is written in R and the main used packages are CovTools, ICtest, parallel, fda, ggplot2 and GA.

## The reason why this project is useful for you
If you're working with complex datasets, such as matrices or images, the following main functions can help you to identify the most central objects within your data. Simply calculate the distances between objects and apply them as an input in our metric depth functions.

### The main metric depth functions are:

MLD_cpp, MHDVJ_cpp, MSD_cpp, MOD2_cpp, MOD3_cpp,

coded in Rcpp. To use them you need to install from source the R-package MetricDepth_1.1 in the subfolder [**All Metric depth function in C++**](https://github.com/vidazamani/Depth-functions-for-Object-Data/tree/main/All%20Metric%20depth%20function%20in%20C%2B%2B).

In addition the project contains several helper functions for data generation, computing p-value for real dataset examples and the simulations performed in the paper.

### Replicating the examples in the paper

The following list shows which files can be used to replicate which simulation study in the manuscript:


### Great! Then How to Use this Project for my own data?

After calling all necessary libraries, and installing MetricDepth_1.1, the next step is to import your dataset and compute the distance between objects in it. For example if you have functional data use metric.lp(YOUR FUNCTIONAL DATA) to achieve this. Finally based on our suggestions in the paper, choose which depth function meets your needs the best and use your computed distance matrix as an input in the chosen depth function, e.g., MOD3_cpp. Now you have found the depths and, among them, the Most central Object in your dataset!

You can see the more detailed examples in [**Real Data example**](https://github.com/vidazamani/Depth-functions-for-Object-Data/tree/15c63bb935e6261f982a280901ab1a48eae834b9/Real%20data%20examples) folder in the current repository.




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
