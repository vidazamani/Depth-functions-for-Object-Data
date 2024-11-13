# Depth Functions for Object Data

This project contains code to reproduce the results of the paper named 'Metric Oja Depth, New Statistical Tool for
Estimating the Most Central Objects' by Zamanifarizhandi, V. and Virta, J.

The goal of the project was to develope a novel measure of statistical depth, the metric Oja depth applicable to any object data. 

All code here is written in R and the main used packages are CovTools, ICtest, parallel, fda, ggplot2 and GA.

## Why is this project useful for you?!
If you're working with complex datasets, such as matrices or images, the following main functions can help you identify the most central objects within your data. Simply calculate the distances between objects and apply them as an input in our metric functions.

### The main metric functions are:

MLD, MHD, MSD, MOD2, MOD3

In addition the project contains several helper functions for data generation, computing p-value for real dataset examples and the simulations performed in the paper.

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
