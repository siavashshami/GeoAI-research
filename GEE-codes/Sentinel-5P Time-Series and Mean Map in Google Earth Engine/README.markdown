# Sentinel-5P Time-Series and Mean Map in Google Earth Engine

This README describes two Google Earth Engine scripts: `Pollution_Chart.js` and `Pollution_Map.js`. Together they load Sentinel-5P TROPOMI NO₂ data for a user-defined region of interest (`table`), build daily mean images, chart the NO₂ column number density time series, display the period-mean NO₂ map clipped to the ROI, and export the mean image to an Earth Engine asset.

---

## Requirements

Both scripts require an existing `table` object in the Code Editor representing the region of interest (ROI). The scripts use:

    ee.ImageCollection("COPERNICUS/S5P/OFFL/L3_NO2")

and select the band:

    NO2_column_number_density

---

## Script 1: `Pollution_Chart.js`

### Overview

This script filters Sentinel-5P NO₂ images for a one-month period, converts the collection into one mean image per day, and charts the daily time series over the ROI.

### Date Range

The script defines:

    var start = ee.Date('2019-03-11');
    var end   = ee.Date('2019-04-09');
    var dateRange = ee.DateRange(start, start.advance(1, 'month'));

The NO₂ collection is filtered with `filterDate(dateRange)` and the `NO2_column_number_density` band is selected.

### Property Copying

Each image is multiplied by 1 and retains its `system:time_start` property. A formatted `date` property is also added:

    image.date().format('yyyy-MM-dd * HH-mm-ss')

### Daily Mean Images

The function `dailyImages(images, startdate, enddate)` builds one mean image per day:

- It computes the number of days using `end.difference(startdate, 'day').round()`.
- It loops over `ee.List.sequence(1, count, 1)`.
- For each day, it filters the collection between `date_start` and `date_end` and computes `.mean()`.
- Each output image is tagged with `system:time_start` equal to the day's date.

The result is stored as:

    var imagePerDay = ee.ImageCollection(dailyImages(modC, start, end));

and printed to the Console.

### Time-Series Chart

The script charts the daily series with:

    ui.Chart.image.series({
      imageCollection: imagePerDay,
      region: table,
      reducer: ee.Reducer.median(),
      scale: 1000,
      xProperty: 'system:time_start'
    })

Chart options:

- Title: `NO2_column_number_density`
- x-axis: `Date`
- y-axis: `mol/m^2`
- Line width: 1
- Point size: 3
- Red trend line option included in the options object

The chart is printed to the Console.

---

## Script 2: `Pollution_Map.js`

### Overview

This script displays the Sentinel-5P NO₂ collection over the ROI, charts the mean time series, maps the period-mean NO₂ clipped to the ROI, and exports the mean image to an Earth Engine asset.

### Collection

The script loads:

    ee.ImageCollection("COPERNICUS/S5P/OFFL/L3_NO2")
      .select('NO2_column_number_density')
      .filterDate('2019-03-11', '2019-04-09');

The ROI geometry is stored as:

    var london_geo = table.geometry();

### Time-Series Chart

The script prints a mean time-series chart using:

    ui.Chart.image.series(collection, london_geo, ee.Reducer.mean(), 30)

### Mean NO₂ Map

The period mean is computed and clipped to the ROI:

    var london_mean = collection.mean().clip(table);

It is added to the map with:

    Map.addLayer(london_mean, {
      region: table,
      min: 0, max: 0.002,
      opacity: 1.0,
      palette: ["black", "blue", "purple", "cyan", "green", "yellow", "red"]},
      'Mean temperature');

Note: the layer name is `Mean temperature` as written in the script, although the data represent NO₂ column number density.

### Export

The clipped mean image is exported to an Earth Engine asset:

    Export.image.toAsset({
      image: london_mean,
      description: 'test',
      assetId: 'test',
      region: table,
      scale: 100,
      crs: 'EPSG:32639',
      maxPixels: 1e13
    });

The export uses:

- `description`: `test`
- `assetId`: `test`
- `region`: `table`
- `scale`: `100`
- `crs`: `EPSG:32639`
- `maxPixels`: `1e13`

The user must start the export task from the Tasks panel and edit `assetId` to a valid path for their account.

---

## How to Run

1. Open either script in the Google Earth Engine Code Editor.
2. Define the `table` object for the ROI.
3. For `Pollution_Chart.js`, adjust `start` and `end` if a different period is needed, then click Run and view the daily series chart.
4. For `Pollution_Map.js`, adjust the date range, map visualization parameters, and export settings (`assetId`, `crs`, `scale`) as needed, then click Run.
5. View the printed chart and map layers in the Console and Map panels.
6. Start the export task from the Tasks panel if the mean NO₂ asset is needed.

`Pollution_Chart.js` produces a daily mean NO₂ time-series chart. `Pollution_Map.js` produces a mean NO₂ time-series chart, a clipped mean NO₂ map layer, and one `Export.image.toAsset` call. Neither script computes RUSLE factors or soil-erosion statistics.