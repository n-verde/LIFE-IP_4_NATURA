/* =============================================================================

==== LIFE IP 4 NATURA - STEP 3/7 =====
Natalia Verde, AUTH/DUTH, 2020 

BRIEF DESCRIPTION:
This GEE script creates S2L2A, S2L1C and S1 seasonal mosaics and exports them in assets.
Mosaics created based on date ranges that are set for each function.

NOTES:
-

HOW TO USE:
1. set path from assets, for the AOI shapefile, to the "AOI" variable (also used in pervious step)
2. set path from assets, for the SNIC clusters, to the "clusters" variable (outputed to assets from previous step)
3. set the temporal mosaic interval as a asting ("month" or "season")
4. set other properties you want in the "PARAMETER SETTINGS" section.
5. in the tasks tab, click RUN to all tasks in order to export the mosaic images to Assets.

*/


//===============PARAMETER SETTINGS===================================================

// ---------- SET AOI ----------
var AOI = ee.FeatureCollection('users/n-verde/shared/LIFE-IP_4_NATURA/GR_boundingBox_buff1km'); 

// ---------- SET CLUSTERS ----------
// computed in step 1
var clusters = ee.Image('users/n-verde/shared/LIFE-IP_4_NATURA/LIFE_10m_OB_clusters'); 

// ---------- SET TEMPORAL MOSAIC INTERVAL ----------
// var interval = 'month';
var interval = 'season';

// ---------- IMAGE EXPORT DETAILS ----------
// var nameSuffix = 'LIFE_10m_M'; // for monthly composites
var nameSuffix = 'LIFE_10m_S'; // for seasonal composites

// ---------- SET EXTERNAL FUNCTIONS SCRIPTS ----------
var S2L2A = require('users/n-verde/auth-shared:LIFE-IP_4_NATURA/3.1-SeasonalMosaicsS2L2A_function'); 
var S2L1C = require('users/n-verde/auth-shared:LIFE-IP_4_NATURA/3.2-SeasonalMosaicsS2L1C_function');
var S1 = require('users/n-verde/auth-shared:LIFE-IP_4_NATURA/3.3-SeasonalMosaicsS1_function'); 

// ---------- OTHER SETTINGS ----------
var scale = 10;
var MP = 1e13; // maxPixels variable (according to AOI and pixel size you want)

//====================================================================================

// ---------- FUNCTIONS ----------

////////////////////////////////////////////////////////////
// This function creates a mask outside the clusters
function maskOut(image) {
    var mask = ee.Image.constant(0).float();
    mask = mask.where(clusters,1);
    var newimage = image.updateMask(mask);
  return newimage;
}

//====================================================================================

//---------------MAIN PROGRAM--------------------------------

print('start');

// ---------- CREATE SEASONAL MOSAICS COLLECTIONS ----------

// S2 L2A ----> 
var ICS2L2C = S2L2A.seasonalMosS2L2A(ee.Date.fromYMD(2018,5,1),ee.Date.fromYMD(2018,5,30),
                                        ee.Date.fromYMD(2019,5,1),ee.Date.fromYMD(2019,10,30),AOI,interval);
var IC_ImageS2 = ICS2L2C.toBands(); // turn collection to ee.Image

// add season property
var oldNames = IC_ImageS2.bandNames();
var slicingNum = ee.String(ICS2L2C.first().get('system:index')).length();
var newNames = oldNames.map(function slicing(x) {
                                      return ee.String(x).slice(slicingNum.add(1));
                                    });
var newNamesClient = newNames.map(function(x) {
                              var e = ee.String('S').cat(x); return e;
                              }); 
IC_ImageS2 = IC_ImageS2.select(oldNames, ee.List(newNamesClient));
print(IC_ImageS2)

var IC_ImageS2_masked = maskOut(IC_ImageS2);

IC_ImageS2_masked = IC_ImageS2_masked.divide(Math.pow(2,8)).float(); 


// S2 L1C ----> 
var ICS2L1C = S2L1C.seasonalMosS2L1C(ee.Date.fromYMD(2019,6,1),ee.Date.fromYMD(2019,8,30),AOI,interval);
var IC_ImageS2L1C = ICS2L1C.toBands(); // turn collection to ee.Image

// add season property
var oldNames = IC_ImageS2L1C.bandNames();
var slicingNum = ee.String(ICS2L1C.first().get('system:index')).length();
var newNames = oldNames.map(function slicing(x) {
                                      return ee.String(x).slice(slicingNum.add(1));
                                    });
var newNamesClient = newNames.map(function(x) {
                              var e = ee.String('S').cat(x); return e;
                              }); 
IC_ImageS2L1C = IC_ImageS2L1C.select(oldNames, ee.List(newNamesClient));
print(IC_ImageS2L1C)

var IC_ImageS2L1C_masked = maskOut(IC_ImageS2L1C);

IC_ImageS2L1C_masked = IC_ImageS2L1C_masked.divide(Math.pow(2,8)).float();


// S1 ----> 
var ICS1 = S1.seasonalMosS1(ee.Date.fromYMD(2019,5,1),ee.Date.fromYMD(2019,5,30),
                                        ee.Date.fromYMD(2019,6,1),ee.Date.fromYMD(2019,10,30),AOI,interval);
var IC_ImageS1 = ICS1.toBands(); // turn collection to ee.Image

// add season/month property
var oldNames = IC_ImageS1.bandNames();
var slicingNum = ee.String(ICS1.first().get('system:index')).length();
var newNames = oldNames.map(function slicing(x) {
                                      return ee.String(x).slice(slicingNum.add(1));
                                    });
var newNamesClient = newNames.map(function(x) {
                              var e = ee.String('S').cat(x); return e;
                              }); 
IC_ImageS1 = IC_ImageS1.select(oldNames, ee.List(newNamesClient));
print(IC_ImageS1)

var IC_ImageS1_masked = maskOut(IC_ImageS1).float();

// ==========================================================================

// --------------------------- DISPLAYS ------------------------
Map.addLayer(AOI, {}, 'AOI', false)
Map.addLayer(IC_ImageS2_masked, {}, 'S2L2A', false)
Map.addLayer(IC_ImageS2L1C_masked, {}, 'S2L1C', false)
Map.addLayer(IC_ImageS1_masked, {}, 'S1', false)

// ==========================================================================

// --------------------------- EXPORTS ------------------------

// TO ASSETS ----------

// Export S2L2A
// Name for your file export
var descr = nameSuffix + '_MosaicS2L2A_masked_float';
Export.image.toAsset({
  image: IC_ImageS2_masked,
  description: descr,
  scale: scale,
  crs: IC_ImageS2_masked.projection(),
  region: AOI,
  maxPixels: MP, //max pixels allowed for download
  pyramidingPolicy: {".default": "mean"}
});

// Export S2L1C
// Name for your file export
var descr = nameSuffix + '_MosaicS2L1C_masked_float';
Export.image.toAsset({
  image: IC_ImageS2L1C_masked,
  description: descr,
  scale: scale,
  crs: IC_ImageS2L1C_masked.projection(),
  region: AOI,
  maxPixels: MP, //max pixels allowed for download
  pyramidingPolicy: {".default": "mean"}
});

// Export S1
// Name for your file export
var descr = nameSuffix + '_MosaicS1_masked_float';
Export.image.toAsset({
  image: IC_ImageS1_masked,
  description: descr,
  scale: scale,
  crs: IC_ImageS1_masked.projection(),
  region: AOI,
  maxPixels: MP, //max pixels allowed for download
  pyramidingPolicy: {".default": "mean"}
});

print('end');
