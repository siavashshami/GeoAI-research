# GRACE TWS Trend Assessment in Google Earth Engine

This README describes only the provided Google Earth Engine script. The script loads GRACE land mass-grid data for a user-defined date range and region of interest (ROI), displays the monthly terrestrial water storage anomaly time series, computes the pixel-wise linear trend in cm/year, maps the period mean and trend, and exports the mean image to an Earth Engine asset.

---

## Overview

Configuration at the top of the script:

    var START_DATE = ee.Date('2009-01-01');
    var END_DATE   = ee.Date('2010-01-01');        // exclusive
    var BAND_NAME  = 'lwe_thickness_csr';
    var EXPORT_CRS = 'EPSG:4326';                  // change to e.g. EPSG:32639 (Iran)
    var EXPORT_SCALE = 1000;                       // metres — realistic for GRACE
    var ASSET_ID   = 'users/YOUR_USERNAME/GRACE_mean_TWS_2009'; // <-- edit this
    var MAP_PALETTE = ['red', 'yellow', 'green', 'cyan', 'purple', 'blue', 'black'];

The script requires an existing `table` object representing the ROI. It centers the map on the ROI at zoom 7 and adds the ROI boundary in red as `ROI`.

---

## Region of Interest

The ROI is defined as:

    var roi = table.geometry();

The script adds the ROI boundary to the map using:

    table.style({color: 'red', fillColor: '00000000', width: 2})

---

## GRACE Collection

The script loads:

    ee.ImageCollection('NASA/GRACE/MASS_GRIDS/LAND')

It filters by `START_DATE` and `END_DATE` and selects `BAND_NAME`, which is set to:

    'lwe_thickness_csr'

The script prints:

- The filtered GRACE collection
- The number of monthly images
- Metadata of the first image

---

## Monthly Time Series

The script creates a monthly time-series chart using:

    ui.Chart.image.series({
      imageCollection: grace,
      region          : roi,
      reducer         : ee.Reducer.mean(),
      scale           : 1000,
      xProperty       : 'system:time_start'
    })

Chart title:

    GRACE Terrestrial Water Storage Anomaly (lwe_thickness_csr)

Axis labels:

- x-axis: `Date`
- y-axis: `Equivalent water thickness (cm)`

The chart uses a line width of 2, point size of 3, and interpolates null values.

---

## Pixel-Wise Linear Trend

The script computes a per-pixel linear trend in cm/year.

For each image, time is converted to years since `START_DATE`:

    var t = ee.Date(img.get('system:time_start'))
              .difference(START_DATE, 'year');

A constant band `t` is added to each image, and the original image properties are copied.

The trend is computed using:

    withTime.select(['t', BAND_NAME])
            .reduce(ee.Reducer.linearFit());

The `scale` band from `linearFit` is renamed:

    tws_trend_cm_per_yr

The trend image is printed to the Console.

---

## Map Layers

The script adds two map layers.

### Mean TWS

The period mean is computed as:

    grace.mean().clip(roi)

Display parameters:

    min: -8.554
    max: -1.054
    palette: MAP_PALETTE

Layer name:

    Mean TWS (cm)

### TWS Linear Trend

The trend image is clipped to the ROI.

Display parameters:

    min: -5
    max: 5
    palette: ['red', 'white', 'blue']

Layer name:

    TWS trend (cm/yr)

---

## Export

The script exports the mean TWS image to an Earth Engine asset:

    Export.image.toAsset({
      image       : meanTWS,
      description : 'GRACE_mean_TWS_2009',
      assetId     : ASSET_ID,
      region      : roi,
      scale       : EXPORT_SCALE,
      crs         : EXPORT_CRS,
      maxPixels   : 1e13
    });

The export uses:

- `description`: `GRACE_mean_TWS_2009`
- `assetId`: value of `ASSET_ID`
- `region`: ROI
- `scale`: `EXPORT_SCALE`
- `crs`: `EXPORT_CRS`
- `maxPixels`: `1e13`

The script starts the export task but does not run it automatically. The user must start the task from the Earth Engine Tasks panel.

---

## How to Run

1. Open the script in the Google Earth Engine Code Editor.
2. Define the `table` object for the ROI.
3. Edit `ASSET_ID` to a valid Earth Engine asset path for your user account.
4. Adjust `START_DATE`, `END_DATE`, `BAND_NAME`, `EXPORT_CRS`, and `EXPORT_SCALE` as needed.
5. Click Run.
6. View the printed collection metadata, monthly time-series chart, trend image, and map layers in the Console and Map panels.
7. Start the export task from the Tasks panel if the mean TWS asset is needed.

The script displays results in the Console and Map and contains one `Export.image.toAsset` call. It does not compute RUSLE factors or soil-loss statistics.