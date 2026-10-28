Map.centerObject(table)
Map.addLayer (table)
var collection = ee.ImageCollection("COPERNICUS/S5P/OFFL/L3_NO2")
  .select('NO2_column_number_density')
  .filterDate('2019-03-11', '2019-04-09');

var london_geo = table.geometry()

print(ui.Chart.image.series(collection, london_geo, ee.Reducer.mean(), 30));
///////////////////////////////////////////////////////////////////////////


Map.centerObject(table)

Map.addLayer (table)
var london_geo = table.geometry()
//Clip to roi
var london_mean = collection.mean().clip(table);

// Add clipped image layer to the map.
Map.addLayer(london_mean, {
  region: table,
  min: 0, max: 0.002,
  opacity: 1.0,
  palette: ["black", "blue", "purple", "cyan", "green", "yellow", "red"]},
  'Mean temperature');
 
  Export.image.toAsset({
  image: london_mean, 
  description: 'test' , 
  assetId: 'test' , 
  region: table, 
  scale: 100, 
  crs:  'EPSG:32639' , 
  maxPixels: 1e13
  
});