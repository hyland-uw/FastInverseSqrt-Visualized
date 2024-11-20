# Testing and instrumenting versions of the fast inverse square root

![An artistic deconstruction](/plots/header.png)

This project is a work in progress. Plot is made with deconstructed.c data pathway.

The code and data here are components of a larger investigation into the history and re-use of the Fast Inverse Square Root, including [the most famous implementation found in Quake III Arena](https://en.wikipedia.org/wiki/Fast_inverse_square_root). The larger project website is here, at [0x5f37642f.com](https://0x5f37642f.com/).

## Mode of operation

Run `make` in the base directory to generate csvs in the data directory. Sampling parameters are set in the sampling-harness.h file.

The project works only on single precision `float` variables, not doubles (As the integer representation is different)

### Compilation note

This project is on an M1 mac, so the makefile has specific settings to use OpenMP. Adjust those for your system.

## Organization

The project is organized by file name. Data is generated by c files, stored as csvs, and operated on by R files. An example is approximated.c -> approximated.csv -> approximated.R. The graphing library ggplot2 is used for most plots, with a legacy file (FISR-plotting.R) that contains base R plots.

All c files use `sampling-harness.h` which contains parameters for generation and common functions.

### Approximated

Computes errors of historical FRSR style approximations over a range of floats, including Quake III's FISR. Data appears broken down by specific approximation:

| ISR_function | input | reference | initial | final |
| --- | --- | --- | --- | --- |
| BlinnISR | 0.3209153 | 1.765244 | 1.858170 | 1.763802 |
| QuakeISR | 0.3209153 | 1.765244 | 1.790600 | 1.764695 |
| withoutDivISR | 0.3209153 | 1.765244 | 1.809832 | 1.763541 |
| optimalFISR | 0.3209153 | 1.765244 | 1.608169 | 1.765233 |

### Deconstructed

Replaces the usual iteration limit of 1-2 Newton-Raphson iterations with iteration to a tolerance, which supports plotting the space for random inputs and magic constants.

| input | reference | initial | final | iters | magic |
| --- | --- | --- | --- | --- | --- |
| 0.574197 | 1.319683 | 0.341227 | 1.319682 | 7 | 1580741785 |
| 0.186377 | 2.316348 | 1.789107 | 2.316348 | 4 | 1594125894 |
| 0.181007 | 2.350457 | 0.591144 | 2.350457 | 7 | 1580466735 |
| 0.629445 | 1.260437 | 0.655565 | 1.260436 | 5 | 1589142732 |

This format allows us to explore "bad" magic constants which produce poor approximations. The number of iterations to reach the right answer is a gross measure of the poor fit of the approximation.

#### Enumerated
Does the unusual job of enumerating a "best" magic constant for a given float. Imagine the world's least efficient lookup table.

| input | reference | initial | final | magic |
| --- | --- | --- | --- | --- |
| 0.5134103 | 1.3956220 | 1.3957580 | 1.3956220 | 1597267869 |
| 0.7665635 | 1.1421570 | 1.1421160 | 1.1421570 | 1597263768 |
| 0.5133603 | 1.3956900 | 1.3957830 | 1.3956900 | 1597267666 |
| 1.2038200 | 0.9114215 | 0.9114995 | 0.9114215 | 1597399913 |

For a given float (the dataset is not ordered by input) we compute the result of many different magic constants. For an input float there is a "best" magic constant which we determine.

### Optimized
performs a grid search of the Newton Raphson constants (~1.5 and ~0.5) over a range of floats and a specific magic constant. Used to show the role of the NR approximation step.

| input | initial | halfthree | halfone | error |
| --- | --- | --- | --- | --- |
| 0.1954819 | 2.301005 | 1.495 | 0.495 | 0.0002761815 |
| 0.1954819 | 2.301005 | 1.495 | 0.496 | 0.0013290450 |
| 0.1954819 | 2.301005 | 1.495 | 0.497 | 0.0023820130 |
| 0.1954819 | 2.301005 | 1.495 | 0.498 | 0.0034350870 |

A basic grid search is used, but the code can be modified for a more sophisticated search. 1.5 and 0.5 are optimal over most of the range, which makes older deviations from that interesting.

### Sliced
Maps the performance of a range of magic constants across sets of floats to visualize slices of the output.

| input | error | magic |
| --- | --- | --- |
| 1.712697 | 1.515474e-03 | 1597466973 |
| 1.712697 | 5.446439e-03 | 1597826541 |
| 1.712697 | 9.883186e-05 | 1597165521 |
| 1.712697 | 3.928469e-03 | 1597712212 |

In contrast to the above, these slices compute a range of outputs for a given float. So for an input `1.712697` many thousands of magic constants are used and the approximation generated. This can be used for plotting curves.

## Future directions

Ideally the end point for this is an R package containing data as well as functions to call the underlying C code at will. What needs to happen for that is:

1. [Convert the existing C code to C++](https://legalizeadulthood.wordpress.com/2007/05/18/refactoring-convert-c-to-c/) to allow the use of [Rcpp](http://dirk.eddelbuettel.com/code/rcpp.html) which has a much better interface than between  R and C.
2. Refactor the existing R and C code from the current mode of batch process to CSV to to an R function like generate_timelines() which would call the underlying C++ code with passed parameters.
3. Generate exemplary datasets and store them in .Rdata format.
4. Convert the whole project into a package with named exported functions and data.

Future needs for the project apart from that:
* Actual documentation
* Once packaging is finished, potentialy moving the project to CRAN

## License
I have not yet chosen a blanket license for these but each of the individual versions are licensed under a variety of terms. Quake III's source code is licensed under the GPL, while fdlibm (which is where the Kahan-Ng softsqrt was published in fixed form) is under a license which may [loosely be described](https://lists.fedoraproject.org/archives/list/legal@lists.fedoraproject.org/thread/2T6RANNIF652RMGG725LNRKT63ALAPN4/) as "MIT". Before borrowing please check the individual example licenses.
