/*=============================================================================

==== LIFE IP 4 NATURA - STEP 7/7 =====
Natalia Verde, AUTH/DUTH, 2020

BRIEF DESCRIPTION:
This GEE script classifies an image stack using a Random Forest algorithm. Exports the classified image to Google Drive.

HOW TO USE:
1. set the paths from the images for classification from assets (all tiles).
2. set path from assets, for the AOI shapefile, to the "AOI" variable (also used in pervious step)
3. set other properties you want in the "PARAMETER SETTINGS" section.
4. in the tasks tab, click RUN in order to export the classified image to Google Drive.

*/

//====================================================================================

//===============PARAMETER SETTINGS===================================================

// ---------- SET IMAGE TO CLASSIFY (All tiles) ----------

// (Monthly)
// only one tile is provided, for a working example of the script: 
var I1 = ee.Image('users/n-verde/shared/LIFE-IP_4_NATURA/LIFE_OB_S_1_IC4classif');
// if you have generated your own images in the previous steps, change this section accordingly,
// and also change variable "tempCol"

// // (Seasonal)
// var I1 = ee.Image('users/n-verde/GR_OB_S_1_IC4classif');
// var I2 = ee.Image('users/n-verde/GR_OB_S_2_IC4classif');
// var I3 = ee.Image('users/n-verde/GR_OB_S_3_IC4classif');
// var I5 = ee.Image('users/n-verde/GR_OB_S_5_IC4classif');
// var I6 = ee.Image('users/n-verde/GR_OB_S_6_IC4classif');
// var I7 = ee.Image('users/n-verde/GR_OB_S_7_IC4classif');
// var I8 = ee.Image('users/n-verde/GR_OB_S_8_IC4classif');
// var I9 = ee.Image('users/n-verde/GR_OB_S_9_IC4classif');
// var I10 = ee.Image('users/n-verde/GR_OB_S_10_IC4classif');
// var I11 = ee.Image('users/n-verde/GR_OB_S_11_IC4classif');
// var I12 = ee.Image('users/n-verde/GR_OB_S_12_IC4classif');
// var I15 = ee.Image('users/n-verde/GR_OB_S_15_IC4classif');
// var I16 = ee.Image('users/n-verde/GR_OB_S_16_IC4classif');

// // (Monthly)
// var I1 = ee.Image('users/n-verde/GR_OB_1_IC4classif');
// var I2 = ee.Image('users/n-verde/GR_OB_2_IC4classif');
// var I3 = ee.Image('users/n-verde/GR_OB_3_IC4classif');
// var I5 = ee.Image('users/n-verde/GR_OB_5_IC4classif');
// var I6 = ee.Image('users/n-verde/GR_OB_6_IC4classif');
// var I7 = ee.Image('users/n-verde/GR_OB_7_IC4classif');
// var I8 = ee.Image('users/n-verde/GR_OB_8_IC4classif');
// var I9 = ee.Image('users/n-verde/GR_OB_9_IC4classif');
// var I10 = ee.Image('users/n-verde/GR_OB_10_IC4classif');
// var I11 = ee.Image('users/n-verde/GR_OB_11_IC4classif');
// var I12 = ee.Image('users/n-verde/GR_OB_12_IC4classif');
// var I15 = ee.Image('users/n-verde/GR_OB_15_IC4classif');
// var I16 = ee.Image('users/n-verde/GR_OB_16_IC4classif');

// ---------- SET AOI ----------
var AOI = ee.FeatureCollection('users/n-verde/shared/LIFE-IP_4_NATURA/GR_boundingBox_buff1km'); 

// ---------- SET TRAINED CLASSIFIER ----------
var RF = require('users/n-verde/auth-shared:LIFE-IP_4_NATURA/6-RFtraining-VI_RUN-EXPORT')

// ---------- IMAGE EXPORT DETAILS ----------
var nameSuffix = 'LIFE';

// ---------- OTHER SETTINGS ----------
var scale = 10;
var prop = 'class_L2';  // the property in the feature collection from which to train model

//====================================================================================

//---------------FUNCTIONS------------------------------

////////////////////////////////////////////////////////////
// This function creates a mask for the area of a FC
function maskPolygons(image) {
  var mask = ee.Image.constant(0).int32();
  mask = mask.paint(featCol, 1);
  var newimage = image.updateMask(mask.not());
  return newimage;
}

////////////////////////////////////////////////////////////
// This function creates a mask outside the area of a FC
function maskOut(image) {
  var mask = ee.Image.constant(0).int32();
  mask = mask.paint(featCol, 1);
  var newimage = image.updateMask(mask);
  return newimage;
}

/////////////////////////////////////////////////////////////
// This function removes useless stdDev per object features from the image for classification.
function removeUselessStdDevs(image) {
  
  // remove features from image
  // texture indices (stdDevs of objects)
  var initialTextureList = ee.List(['CorB2','PANTEX','CorBCI','CorVV','CorVH']);
  var season = ee.String('S2');
  var S2textureList = initialTextureList.map(function(x) { return season.cat(x)});
  var season = ee.String('S3');
  var S3textureList = initialTextureList.map(function(x) { return season.cat(x)});
  var season = ee.String('S4');
  var S4textureList = initialTextureList.map(function(x) { return season.cat(x)});
  var allSeasontextureList_final = S2textureList.cat(S3textureList).cat(S4textureList); // stdDevs
  // print(allSeasontextureList_final)
  
  // DEM features (stdDevs of objects)
  var demList_final = ee.List(['elev','','slope','TRASP']);
  
  var allStdDevFeatures = allSeasontextureList_final.cat(demList_final);
  
  var imageWithoutStdDevs = image.select(image.bandNames().removeAll(allStdDevFeatures));
  
  return imageWithoutStdDevs;
}

//====================================================================================

//---------------MAIN PROGRAM--------------------------------

print('start');

// Mosaic tiles to create Image
// only one tile is provided, for a working example of the script: 
var tempCol = ee.ImageCollection.fromImages([I1]);
// var tempCol = ee.ImageCollection.fromImages([I1,I2,I3,I5,I6,I7,I8,I9,I10,I11,I12,I15,I16]);

var IC4classif = tempCol.mosaic();

print('IC4classif bands', IC4classif.bandNames())

// ---------- RUN CLASSIFIER TO CLASSIFY IMAGE STACK ----------
var classified = IC4classif.classify(RF.classifier).uint8();
// print('classified image:', classified)

//-------------- DISPLAYS --------------------------------------

Map.addLayer(IC4classif, {}, 'IC', false)

// CLASSIFIED WITH 21 CLASSES
// Define an SLD style of discrete intervals to apply to the classified product
var sld_intervals =
'<RasterSymbolizer>' +
  ' <ColorMap  type="intervals" extended="false" >' +
    '<ColorMapEntry color="#ff0000" quantity="1" label="Dense to medium dense Urban Fabric"/>' +
    '<ColorMapEntry color="#ff7d7d" quantity="2" label="Low density Urban Fabric"/>' +
    '<ColorMapEntry color="#ffffa8" quantity="3" label="Arable land"/>' +
    '<ColorMapEntry color="#cccc00" quantity="4" label="Permanent crops"/>' +
    '<ColorMapEntry color="#99ff99" quantity="5" label="Temperate deciduous forests"/>' +
    '<ColorMapEntry color="#00ff00" quantity="6" label="Mediterranean deciduous forests"/>' +
    '<ColorMapEntry color="#00cc66" quantity="7" label="Floodplain forests"/>' +
    '<ColorMapEntry color="#006600" quantity="8" label="Temperate mountainous coniferous forests"/>' +
    '<ColorMapEntry color="#009900" quantity="9" label="Mediterranean coniferous forests"/>' +
    '<ColorMapEntry color="#4f6228" quantity="10" label="Mediterranean sclerophyllous forests"/>' +
    '<ColorMapEntry color="#d1f616" quantity="11" label="Mixed Forest"/>' +
    '<ColorMapEntry color="#e2f5c1" quantity="12" label="Grasslands"/>' +
    '<ColorMapEntry color="#f2cca6" quantity="13" label="Moors and heathland"/>' +
    '<ColorMapEntry color="#f3e593" quantity="14" label="Sclerophyllous vegetation"/>' +
    '<ColorMapEntry color="#bf976f" quantity="15" label="Sparsely vegetated areas"/>' +
    '<ColorMapEntry color="#e6e6e6" quantity="16" label="Beaches, dunes, sands"/>' +
    '<ColorMapEntry color="#aaaaaa" quantity="17" label="Bare rocks, burnt areas, mines, dump, land without current use"/>' +
    '<ColorMapEntry color="#99ccff" quantity="18" label="Marshes"/>' +
    '<ColorMapEntry color="#403151" quantity="19" label="Peat bogs"/>' +
    '<ColorMapEntry color="#0066ff" quantity="20" label="Marine"/>' +
    '<ColorMapEntry color="#80f2e6" quantity="21" label="Rivers and Lakes"/>' +
  '</ColorMap>' +
'</RasterSymbolizer>';

Map.addLayer(classified.sldStyle(sld_intervals), {}, 'classified', false)

// ==========================================================================

// --------------------------- EXPORTS ------------------------

// Export classified image (classified)
var descr = nameSuffix + '_classified';
Export.image.toDrive({
  image: classified,
  description: descr,
  folder: 'GEE_output',
  scale: scale,
  region: AOI,
  crs: classified.projection(),
  maxPixels: 1e13 //max pixels allowed for download
});


print('end');