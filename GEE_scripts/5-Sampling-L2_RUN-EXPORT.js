/* =============================================================================

==== LIFE IP 4 NATURA - STEP 5/7 =====
Natalia Verde, AUTH/DUTH, 2020

BRIEF DESCRIPTION:
This GEE script exports samples from an image, based on the sample points uploaded by the user.

HOW TO USE:
1. set the appropriate training points and validation points paths in assets.
2. set the paths from the images for sampling from assets (all tiles).
3. set other properties you want in the "PARAMETER SETTINGS" section.
4. in the tasks tab, click RUN to all tasks in order to export the samples to Assets.

*/

//====================================================================================

//===============PARAMETER SETTINGS===================================================

// ---------- SET SAMPLE POINTS ----------
// computed in step 2

// TRAINING ----------
// manual
// var trainingFeatures = ee.FeatureCollection('users/n-verde/auth/LIFE_IP_4_NATURA/LIFE_GR_points_GM_NV_500');
var trainingFeatures = ee.FeatureCollection('users/n-verde/shared/LIFE-IP_4_NATURA/LIFE_points_all_500');

// VALIDATION ----------
// manual
// var validationFeatures = ee.FeatureCollection('users/n-verde/auth/LIFE_IP_4_NATURA/val_GM_LPIS_CORINE_32634_NV');
// automated validation samples are provided, only for a working example of the script:
var validationFeatures = ee.FeatureCollection('users/n-verde/shared/LIFE-IP_4_NATURA/LIFE_points_all_150');

// ---------- SET IMAGE COLLECTION TILES ----------

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


// ---------- EXPORT DETAILS ----------
var nameSuffix = 'LIFE_OB';

// ---------- OTHER SETTINGS ----------
var prop = 'class_L2';  // the property in the feature collection from which to train model (L1 or L2)
var scale = 10;
var exp = false; // whether to export to drive or not (otherwise will only export to assets)
var MP = 1e9; // maxPixels variable (according to AOI and pixel size you want)

//====================================================================================

//---------------MAIN PROGRAM--------------------------------

print('start');

// Mosaic tiles to create Image

// only one tile is provided, for a working example of the script: 
var tempCol = ee.ImageCollection.fromImages([I1]);
// var tempCol = ee.ImageCollection.fromImages([I1,I2,I3,I5,I6,I7,I8,I9,I10,I11,I12,I15,I16]);

var IC4classif = tempCol.mosaic();

print('IC4classif bands', IC4classif.bandNames())

// ---------- create training and validation data ----------

// sample the training image, to obtain spectral information for training points
var TrainingSamples = ee.Image(IC4classif).sampleRegions({
  collection: trainingFeatures,
  properties: [prop],
  scale: scale,
  tileScale: 16, //5,
  geometries: true
  });
  
// sample the training image, to obtain spectral information for validation points
var ValidationSamples = ee.Image(IC4classif).sampleRegions({
  collection: validationFeatures,
  properties: [prop],
  scale: scale,
  tileScale: 16, //5,
  geometries: true
  });
 
// ==========================================================================

// --------------------------- DISPLAYS ------------------------
Map.addLayer(IC4classif, {}, 'IC4classif', false)

// ==========================================================================

// --------------------------- EXPORTS ------------------------

// TO ASSETS ----------

// ---------- export training sample points to Assets ----------
var fileName = nameSuffix.concat('_').concat('TrainingSamples');
Export.table.toAsset({
  collection: TrainingSamples,
  description: fileName,
  assetId: fileName
});

// ---------- export validation sample points to Assets ----------
var fileName = nameSuffix.concat('_').concat('ValidationSamples');
Export.table.toAsset({
  collection: ValidationSamples,
  description: fileName,
  assetId: fileName
});

// TO DRIVE ----------

if (exp === true) {

  // ---------- export training sample points to Drive ----------
  var fileName = nameSuffix.concat('_').concat('TrainingSamples');
  Export.table.toDrive({
    collection: TrainingSamples,
    description: fileName,
    folder: 'GEE_output',
    fileFormat: 'SHP'
  });
  
  // ---------- export validation sample points to Drive ----------
  var fileName = nameSuffix.concat('_').concat('ValidationSamples');
  Export.table.toDrive({
    collection: ValidationSamples,
    description: fileName,
    folder: 'GEE_output',
    fileFormat: 'SHP'
  });

}

print('end');
