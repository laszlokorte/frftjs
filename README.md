# Fast Fractional Fourier Transform in JavaScript

expects signal as array of shape `[[re, im], [re, im], [re, im], ...]` with the length being an integer power of 2.

Base on https://nalag.cs.kuleuven.be/research/software/FRFT/ which in turn is based on:

* R. Tao, G. Liang, X. Zhao.
An efficient FPGA-based implementation of fractional Fourier transform algorithm.
J. Signal Processing Systems, 60(1):47--58, 2010. DOI 10.1007/s11265-009-0401-0

* A. Bultheel.
A two-phase implementation of the fractional Fourier transform.
Report TW 588, Dept. Computer Science, K.U.Leuven, March 2011. 