# RUSLE-Based Soil Erosion Analysis in Google Earth Engine

This README describes only the provided Google Earth Engine script. The script computes RUSLE soil loss for a user-defined AOI and year range, displays RUSLE factor layers, prints land-cover histograms, performs sanity checks, prints factor and RUSLE statistics, computes multi-year land-cover erosion statistics, and displays charts.

---

## Overview

The script implements the RUSLE equation:

    A = R × K × LS × C × P

Configuration at the top of the script:

    var zone = ee.FeatureCollection('projects/ee-siavashshami/assets/Cs');
    var startYear = 2016;
    var endYear   = 2016;
    var years = ee.List.sequence(startYear, endYear);
    var climateZoneName = 'Cs';
    var analysisScale = 2000;

The script centers the map on the AOI at zoom 5 and adds the AOI boundary in red as `AOI: Cs`.

---

## Land Cover

MODIS `MCD12Q1` `LC_Type1` is filtered by year and AOI, clipped to the AOI, and reclassified from 17 original classes into 12 merged classes:

    EvergreenForest, DeciduousForest, MixedForest, Shrublands,
    Savanna, Grasslands, Wetlands, Agriculture, Urban,
    SnowIce, Barren, Water

The reclassification remaps:

    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
    to
    [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 8, 9, 8, 10, 11, 12]

Classes 0, 254, and 255 are masked out. The merged image is renamed `merged_LC`.

For the first year, the script prints:
- Original land-cover image
- Merged land-cover image
- Original land-cover histogram (17 classes)
- Merged land-cover histogram (12 classes)

The script adds:
- `LC Original (17 classes)` layer
- `LC Merged (12 classes)` layer

---

## RUSLE Factors

### R-Factor: CHIRPS Daily

Method: CHIRPS daily precipitation → monthly sums → annual rainfall → Modified Fournier Index (MFI) → R-factor.

The script:
- Filters `UCSB-CHG/CHIRPS/DAILY` by AOI and year.
- Computes monthly sums using `temporalMonthlySum`.
- Computes annual rainfall `P` and `P_safe = P.max(1)`.
- Computes monthly squares, sums them, and divides by `P_safe` to get `MFI`.
- Computes:

    R = 1.735 * pow(10, 1.5 * log10(MFI) - 0.8188)

- Sets R to 0 where annual rainfall equals 0.
- Clips R to the AOI.

Reference: Renard et al. (1997).

### K-Factor: OpenLandMap USDA Texture + NRCS K Values

The script uses `OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02`, selects `b0`, clips to the AOI buffered by 10000 m, resamples bilinear, and reprojects to `EPSG:4326` at `analysisScale`.

It remaps USDA texture classes 1–12 to K values:

    [0.41, 0.45, 0.31, 0.34, 0.36, 0.29, 0.28, 0.26, 0.20, 0.25, 0.17, 0.15]

K is masked to `K > 0` and `K <= 0.69`. The function returns:
- `K`
- `texture_class`
- `K_valid_mask`

Reference: Renard et al. (1997), USDA NRCS.

### C-Factor: NDVI-Based, Hemisphere-Aware

The script uses MODIS `MOD13Q1` NDVI.

It computes:
- Annual NDVI median for the year, scaled by 0.0001 and clamped 0–1.
- Hemisphere from AOI bounds centroid latitude.
- Growing-season NDVI:
  - Southern hemisphere: November–December of `startYear` merged with January–May of `endYear`.
  - Northern hemisphere: April–October of `startYear`.
- Growing-season NDVI scaled by 0.0001 and clamped 0–1.

C-factor components:

    C_annual  = exp(-2 * ndvi / (1 - ndvi))
    C_growing = exp(-2 * ndvi / (1 - ndvi))

Final C:

    C_scientific = C_annual * 0.3 + C_growing * 0.7

clamped to `0.001–0.95` and renamed `C_scientific`.

The function also returns:
- `growing_NDVI_raw`
- `annual_NDVI_raw`

Reference: Van der Knijff et al. (2000).

### LS-Factor: Moore & Burch (1986) with HydroSHEDS

The script uses:
- `USGS/SRTMGL1_003` DEM
- `WWF/HydroSHEDS/15ACC` flow accumulation

DEM is clipped to AOI, resampled bilinear, and reprojected to `EPSG:4326` at `analysisScale`.

Slope is computed with `ee.Terrain.slope`, then converted to radians.

Flow accumulation is converted to number of upstream cells plus current cell:

    N = flowAcc_raw / 100 + 1

Specific catchment area:

    As = N * cellSize_m

where `cellSize_m = analysisScale`.

L-factor:

    L = As > 1 ? pow(As / 22.13, 0.4) : 1

S-factor:

    S = s > 0.0896 ? (16.8 * sin(s) - 0.50) : (10.8 * sin(s) + 0.03)

where `s` is slope in radians.

Final LS:

    LS = L * S

renamed `LS_moore_burch`.

The function returns:
- `LS_moore_burch`
- `L_moore_burch`
- `S_rusle`
- `As_m`
- `N_cells`

Reference: Moore & Burch (1986); McCool et al. (1987).

### P-Factor: Slope-Based Uniform Expression

The script derives:

    slopePercent = S_rusle * 100

and renames it `slope_percent`.

P-factor is computed with the expression:

    (s < 2)  ? 0.10 :
    (s < 5)  ? 0.12 :
    (s < 8)  ? 0.14 :
    (s < 12) ? 0.25 :
    (s < 16) ? 0.30 :
    (s < 20) ? 0.50 : 1.00

where `s = slopePercent.divide(100)`.

The result is renamed `P_scientific` and clipped to the merged land-cover geometry. The function also returns `P_slope_only`.

Reference: Wischmeier & Smith (1978).

---

## Final RUSLE Calculation

For the first year, the script computes:

    R_first
    K_first
    cFactor_first
    ls_first
    P_first

Final RUSLE:

    rusle_first = R_first
      * K_first.K
      * ls_first.LS_moore_burch
      * cFactor_first.C_scientific
      * P_first.P_scientific

renamed `RUSLE_ton_ha`.

---

## Map Layers Added

The script adds:

- `AOI: Cs`
- `LC Original (17 classes)`
- `LC Merged (12 classes)`
- `R-factor 2016`
- `K-factor SCIENTIFIC (NRCS) 2016`
- `USDA Texture Classes`
- `C-factor 2016 (Van der Knijff)`
- `LS-factor Moore&Burch CORRECTED 2016`
- `P-factor UNIFORM (slope-based) 2016`
- `COMPLETE RUSLE Soil Loss 2016 (t/ha/yr)`

---

## Sanity Checks

The script performs two sanity checks and prints results and charts.

### 1. Slope Class vs. Mean Soil Loss

Slope breakpoints:

    [0, 2, 5, 8, 12, 16, 20, 100]

Slope class names:

    ['0-2%', '2-5%', '5-8%', '8-12%', '12-16%', '16-20%', '>20%']

The script creates a slope class image from `slope_percent`, masks class 0, groups mean RUSLE by slope class, prints a table, and displays a column chart titled `Mean Soil Loss by Slope Class`.

Expected behavior: soil loss should increase with slope.

### 2. NDVI Class vs. Mean Soil Loss

NDVI breakpoints:

    [-1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0]

NDVI class names:

    ['<0.1', '0.1-0.2', '0.2-0.3', '0.3-0.4', '0.4-0.5', '0.5-0.6', '0.6-0.7', '>0.7']

The script uses annual NDVI scaled by 0.0001, creates NDVI classes, masks class 0, groups mean RUSLE by NDVI class, prints a table, and displays a column chart titled `Mean Soil Loss by NDVI Class`.

Expected behavior: soil loss should decrease with increasing NDVI.

---

## Statistics Printed

The script prints min, max, mean, and standard deviation for:

- R-factor
- K-factor
- C-factor
- LS-factor
- P-factor

For total RUSLE, it prints:

- Min
- Max
- Mean
- Median
- Standard deviation

It also prints AOI centroid latitude and hemisphere.

---

## Multi-Year Land Cover Erosion Statistics

The script computes RUSLE for each year in `years = sequence(startYear, endYear)` using `computeRUSLEForYear`.

For each year and each merged land-cover class, `getLandCoverErosionStats` computes:

- `year`
- `climate_zone`
- `landcoverClass`
- `pixel_count`
- `total_class_area_km2`
- `total_class_area_ha`
- `mean_erosion_ton_ha_yr`
- `total_erosion_ton_yr`

Mean erosion is set to 0 for:

    Water, SnowIce, Wetlands, Urban

The script builds `allYearStats` and prints `MULTI-YEAR LAND COVER EROSION RESULTS:` using `allYearStats.limit(10)`.

---

## Charts

The script displays two grouped charts:

### Area Chart

Title: `Land Cover Area by Class (km²) - Cs (2016–2016)`

Grouped column chart using:
- x: `year`
- y: `total_class_area_km2`
- series: `landcoverClass`

### Erosion Rate Chart

Title: `Mean Soil Erosion Rate by Land Cover (t/ha/yr) - Cs (2016–2016)`

Grouped line chart using:
- x: `year`
- y: `mean_erosion_ton_ha_yr`
- series: `landcoverClass`

---

## How to Run

1. Open the script in the Google Earth Engine Code Editor.
2. Adjust the AOI asset path, `startYear`, `endYear`, `climateZoneName`, and `analysisScale`.
3. Click Run.
4. View map layers, printed histograms, sanity-check tables and charts, factor statistics, multi-year land-cover erosion results, and grouped charts in the Console and Map panels.

The script displays results in the Console and Map. It does not contain an export call.