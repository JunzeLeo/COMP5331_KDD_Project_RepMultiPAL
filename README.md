# COMP5331_KDD_Project_RepMultiPAL

## Setup 
This algorithm is based on matrix profile provided by the python interface of [scamp](https://github.com/zpzim/SCAMP).
So you need a python and pyscamp package before start:
```
pip install pyscamp
```

## Quick start
### 1.code

'demoRun.m' is used to implement the multi-way join algorithms and visualize the final results.

'RepMultiPAL.m' is the improved version of the original MultiPAL algorithm.

'correlation_functions' includes three proposed evaluation metrics, i.e. DCO, kDDTW and IPD.

### 2. benchmark_data

The datasets used in our experiments, i.e. Sinusoidal, Birds, Trace and PAMAP.

The '3.1Benchmark Dataset' section in the report describes more detailed introduction and statistics for each dataset.

### 3. case_study_data

'11stocks_volume.mat' and '11stocks_close.mat' include 11 stock time series sequences related to the the number of shares traded and the closing price respectively.

'10region_covid.mat' includes 10 time series sequences related to the death number of Covid-19 in different counties in US.

The '3.2 Stock Dataset' and '3.3 COVID-19 Dataset' sections in the report describes more detailed introduction and statistics for each dataset.

