# Cloud Advection Model
Implementation of the Cloud Advection Model as described by Ranalli et al. Simulates the smoothing of irradiance timeseries by a spatially distributed plant based on cloud advection. For more details see these references:

* J. Ranalli, E.E.M. Peerlings and T. Schmidt, “Cloud Advection and Spatial Variability of Solar Irradiance,” _2020 47th IEEE Photovoltaic Specialists Conference (PVSC)_, 2020. [https://doi.org/10.1109/PVSC45281.2020.9300700](https://doi.org/10.1109/PVSC45281.2020.9300700)
* J. Ranalli and E.E.M. Peerlings, "Cloud Advection Model of Solar Irradiance Smoothing by Spatial Aggregation". _(forthcoming, under review by Journal of Renewable and Sustainable Energy)_.

### camsmoothing.py
The code to run the model.

### demo.py
A demonstration of the model in action based upon a uniformly distributed 1-d plant.

### demo_2d.py
A demonstration of the model based upon a random 5km by 5km distribution of points that makeup a plant. As the model relies on a 1-d representation of the plant, these points are projected onto a line representing the cloud motion over the reference point. 

### livermore.csv
Sample data is provided in ```livermore.csv```. This dataset is a CSV file reproduction of the demo data included with the Wavelet Variability Model. Please refer to the source repository for the WVM in PVLIB MATLAB for info and license [https://github.com/sandialabs/wvm](https://github.com/sandialabs/wvm). 

## Dependencies
Running the model requires ```numpy```. The demo additionally requires ```matplotlib``` for the plots.