// Import image collection
var NO2 = ee.ImageCollection("COPERNICUS/S5P/OFFL/L3_NO2");

var start = ee.Date('2019-03-11');
var end = ee.Date('2019-04-09')
var dateRange = ee.DateRange(start, start.advance(1, 'month'));

var modL = NO2.filterDate(dateRange).select('NO2_column_number_density');

var modC = modL.map(function(image) {
  return image.multiply(1)
    .copyProperties(image, ['system:time_start'])
    .set({date: image.date().format('yyyy-MM-dd * HH-mm-ss')});
});

// create one image per day
function dailyImages (images,startdate,enddate){
   startdate = ee.Date(startdate);
   enddate = ee.Date(enddate);
   var originDate = startdate;
   var count = end.difference(startdate, 'day').round();
   var listOfYear = ee.List.sequence(1, count, 1);
   return ee.ImageCollection(listOfYear.map(function(i){
     i = ee.Number(i);
     var date_start = originDate.advance(i.add(-1), 'day');
     var date_end = originDate.advance(i, 'day');
     var meanImage = images.filterDate(date_start,date_end).mean();
     return meanImage
     .set('system:time_start', date_start.format('yyyy-MM-dd'))
   }));
 }
 
 var imagePerDay = ee.ImageCollection(dailyImages (modC,start,end));
 print(imagePerDay)


// Chart the time-series
var temp_trend = ui.Chart.image.series({
  imageCollection: imagePerDay,
  region: table,
  reducer: ee.Reducer.median(),
  scale:1000,
   xProperty: 'system:time_start'})
  .setOptions({
    title: 'NO2_column_number_density',
    hAxis : {title : 'Date'},
    lineWidth : 1,
    pointSize : 3 ,
    trend_no2 : {1:{color:'red'
   }},
     vAxis : {title : 'mol/m^2'}});
print(temp_trend);
