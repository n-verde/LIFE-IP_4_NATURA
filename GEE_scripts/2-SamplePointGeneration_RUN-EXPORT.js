/*=============================================================================

==== LIFE IP 4 NATURA - STEP 2/7 =====
Natalia Verde, AUTH/DUTH, 2020 

BRIEF DESCRIPTION:
This GEE script creates functions for random polygon sampling from a segmented image and 
various vector layers.

HOW TO USE:
1. upload the reference data as thematic rasters, reclassified to the appropriate classes
  (classes should be in integers starting from "1" and "0" for no-data):
  a) HRL imperviousness layer, set path from assets to "HR" variable (contains only classes "1" and "2")
  b) CORINE/LPIS merged layer, set path from assets to "LP" variable (contains only classes "3" and "4")
  c) Natura2000 layer, set path from assets to "NA" variable
2. set path from assets, for the SNIC clusters, to the "clusters" variable (outputed to assets from previous step)
3. set path from assets, for the SNIC seeds, to the "seeds" variable (outputed to assets from previous step)
4. upload your validation points to assets, as a raster and import to the script. Rename it to "valRaster".
   If you don't have validation samples yet, set valRaster to a blank image (see included code line bellow)
   and run the code with nPoints = the desired number of validation points.
   Export and rasterize the validation points localy on your PC and then upload them to assets as a raster.
5. set path from assets, for the AOI shapefile, to the "aoi" variable (also used in pervious step)
6. set the property in your AOI shapefile, which uniquely defines your AOI polygon.
7. set other properties you want in the "PARAMETER SETTINGS" section.

 */

//====================================================================================

//===============PARAMETER SETTINGS===================================================

// ---------- SET REFERENCE DATASETS ----------
// Rasterized locally on PC and then uploaded
var HR = ee.Image('users/n-verde/shared/LIFE-IP_4_NATURA/HRL_IMD_2015_020m_EPSG3035_E50N10_E50N20');
var LP = ee.Image('users/n-verde/shared/LIFE-IP_4_NATURA/ILOTS-Corine_2018_LIFE_R');
var NA = ee.Image('users/n-verde/shared/LIFE-IP_4_NATURA/buffered_Natura_L2_R');

// ---------- SET CLUSTERS & SEEDS ----------
// computed in step 1
var clusters = ee.Image('users/n-verde/shared/LIFE-IP_4_NATURA/LIFE_10m_OB_clusters'); 
var seeds = ee.Image('users/n-verde/shared/LIFE-IP_4_NATURA/LIFE_10m_OB_seeds');

// ---------- SET AOI ----------
var aoi = ee.FeatureCollection('users/n-verde/shared/LIFE-IP_4_NATURA/GR_boundingBox_buff1km'); 
var idName = 'Id'; 

// ---------- SET RASTER WITH VALIDATION POINTS ----------
// to run the script for the first time if you haven't generated validation samples yet.
// var valRaster = ee.Image().byte().paint(aoi.geometry(), 0);

// binary raster with 1 where point exists and 0 otherwise. Rasterized locally on PC and then uploaded
var valRaster = ee.Image('users/n-verde/auth/LIFE_IP_4_NATURA/val_GM_LPIS_CORINE_32636_NV_R');

// ---------- OTHER SETTINGS ----------
//(it's not recomended to change these values)
var nPoints = 500; // number of points per class (150 points per class are validation, so training will be 500)
var proj = 'EPSG:4326'; // the projection of the broader area 
var res = 10; // spatial resolution of image, for stratified sampling
var exportFlag = true; // for exporting the samples to  Drive

// ---------- IMAGE EXPORT DETAILS ----------
var nameSuffix = 'LIFE';

//====================================================================================

// ---------- FUNCTIONS ----------

////////////////////////////////////////////////////////////
// This function creates a mask outside the geometry
function maskOut(image, geometry) {
  
  var mask = ee.Image.constant(0).byte();
  mask = mask.paint(geometry,1);
  var newimage = image.updateMask(mask);
  return newimage;
  
}

//====================================================================================

//---------------MAIN PROGRAM--------------------------------

print('start');

// make aoi a mask
aoi = aoi.map(function f(x) {
                          return x.set(idName,1);  
                            });
var aoiI = aoi.reduceToImage([idName],ee.Reducer.first()); 
// Map.addLayer(aoi, {}, 'aoiI')

var clusters = clusters.rename('clusters');
clusters = clusters.updateMask(aoiI);
// Map.addLayer(clusters, {}, 'clusters')

// find validation segments to mask and not take training samples there
var val_seg = clusters.addBands(valRaster).reduceConnectedComponents({
                                            reducer: ee.Reducer.max(), 
                                            labelBand: 'clusters', 
                                            maxSize: 100});
// Map.addLayer(val_seg, {min:0, max:1},'validation segments')

/////////////////////////////////////////////////////////////////////////////////////
// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ NATURA \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

var cl_NA_L2_I = NA.rename('class_L2');
cl_NA_L2_I = cl_NA_L2_I.unmask();  // turn masked areas to 0, otherwise mask affects stdDev
Map.addLayer(cl_NA_L2_I, {}, 'Natura', false)

// ---------- L2
// Compute per-cluster L2 class stdDev
var L2_segStdDev = clusters.addBands(cl_NA_L2_I).reduceConnectedComponents({
                                                  reducer: ee.Reducer.stdDev(), 
                                                  labelBand: 'clusters', 
                                                  maxSize: 100});

/*  !!!
You should specify the same neighborhood to both snic and reduceConnectedComponents, 
otherwise one of them might get truncated while the other one doesn't 
*/

// Compute per-cluster class
var L2_seg = cl_NA_L2_I.addBands(clusters).reduceConnectedComponents({
                                            reducer: ee.Reducer.mode(), 
                                            labelBand: 'clusters', 
                                            maxSize: 100});

// select only segments with one class (clear class segments)
// stdDev will be 0
var stdDevNot0 = L2_segStdDev.neq(0); // mask segments with stdDev <> 0
var L2_seg_cl = L2_seg.where(stdDevNot0, 0);
var gt28 = L2_seg_cl.gt(22); // mask segments with cl > number of Natura classes
L2_seg_cl = L2_seg_cl.where(gt28, 0);
var zero = L2_seg_cl.eq(0); // mask segments with cl = 0
L2_seg_cl = L2_seg_cl.where(zero, 0);

// SAMPLE SEGMENTS FROM objectPropertiesImage ----------

// create a mask containing "clean" L2 segments only
var samp = L2_seg_cl.gt(0);

// mask clusters and keep only the ones over clean segments
var clusters_cl = clusters.where(samp.not(), 0);
clusters_cl = clusters.updateMask(samp);
var clusters_cl_class = L2_seg_cl.updateMask(samp).toUint8();
Map.addLayer(clusters_cl_class, {}, 'clean segments Natura', false)

// mask validation segments
var clusters_cl_class2 = clusters_cl_class.updateMask(val_seg.not());
// Map.addLayer(clusters_cl_class2, {}, 'clean segments Natura (removed validation)')

// mask area out of seeds, to only take seeds as samples
var clusters_cl_class3 = clusters_cl_class2.updateMask(seeds);
// Map.addLayer(clusters_cl_class3, {}, 'clean Natura seeds(removed validation)')

// cerate random points within clean segments
var points_NA = clusters_cl_class3.stratifiedSample({
                               numPoints: nPoints,
                               region: aoi,
                               scale: res,
                               projection: clusters_cl_class3.projection(),
                               tileScale: 16,
                               geometries: true
                               });
Map.addLayer(points_NA, {}, 'Natura points', false)

// // ---------- export points SHPs to Assets ----------
// var fileName = nameSuffix.concat('_').concat('points_NA');
// Export.table.toAsset({
//   collection: points_NA,
//   description: fileName,
//   assetId: fileName
// });

// // ---------- export points SHPs to Drive ----------
// if (exportFlag === true) {

//   var fileName = nameSuffix.concat('_').concat('points_NA');
//   Export.table.toDrive({
//     collection: points_NA,
//     description: fileName,
//     folder: 'GEE_output',
//     fileFormat: 'SHP'
//   });
  
// }


// /////////////////////////////////////////////////////////////////////////////////////
// // \\\\\\\\\\\\\\\\\\\\\\\\\\\\\ LPIS + CORINE \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

var cl_LP_L2_I = LP.rename('class_L2');
cl_LP_L2_I = cl_LP_L2_I.unmask();  // turn masked areas to 0, otherwise mask affects stdDev
Map.addLayer(cl_LP_L2_I, {}, 'LPIS', false)

// ---------- L2
// Compute per-cluster L2 class stdDev
var L2_segStdDev = clusters.addBands(cl_LP_L2_I).reduceConnectedComponents({
                                                  reducer: ee.Reducer.stdDev(), 
                                                  labelBand: 'clusters', 
                                                  maxSize: 100});

/*  !!!
You should specify the same neighborhood to both snic and reduceConnectedComponents, 
otherwise one of them might get truncated while the other one doesn't 
*/

// Compute per-cluster class
var L2_seg = cl_LP_L2_I.addBands(clusters).reduceConnectedComponents({
                                            reducer: ee.Reducer.mode(), 
                                            labelBand: 'clusters', 
                                            maxSize: 100});

// select only segments with one class (clear class segments)
// stdDev will be 0
var stdDevNot0 = L2_segStdDev.neq(0); // mask segments with stdDev <> 0
var L2_seg_cl = L2_seg.where(stdDevNot0, 0);
var gt28 = L2_seg_cl.gt(5); // mask segments with cl > number of LPIS classes
L2_seg_cl = L2_seg_cl.where(gt28, 0);
var zero = L2_seg_cl.eq(0); // mask segments with cl = 0
L2_seg_cl = L2_seg_cl.where(zero, 0);

// SAMPLE SEGMENTS FROM objectPropertiesImage ----------

// create a mask containing "clean" L2 segments only
var samp = L2_seg_cl.gt(0);

// mask clusters and keep only the ones over clean segments
var clusters_cl = clusters.where(samp.not(), 0);
clusters_cl = clusters.updateMask(samp);
var clusters_cl_class = L2_seg_cl.updateMask(samp).toUint8();
Map.addLayer(clusters_cl_class, {}, 'clean segments LPIS', false)

// mask validation segments
var clusters_cl_class2 = clusters_cl_class.updateMask(val_seg.not());
// Map.addLayer(clusters_cl_class2, {}, 'clean segments Natura (removed validation)')

// mask area out of seeds, to only take seeds as samples
var clusters_cl_class3 = clusters_cl_class2.updateMask(seeds);
// Map.addLayer(clusters_cl_class3, {}, 'clean Natura seeds(removed validation)')

// cerate random points within clean segments
var points_LPCOR = clusters_cl_class3.stratifiedSample({
                               numPoints: nPoints,
                               region: aoi,
                               scale: res,
                               projection: clusters_cl_class3.projection(),
                               tileScale: 16,
                               geometries: true
                               });
Map.addLayer(points_LPCOR, {}, 'LPIS + CORINE points', false)

// // ---------- export points SHPs to Assets ----------
// var fileName = nameSuffix.concat('_').concat('points_LPCOR');
// Export.table.toAsset({
//   collection: points_LPCOR,
//   description: fileName,
//   assetId: fileName
// });

// // ---------- export points SHPs to Drive ----------
// if (exportFlag === true) {

//   var fileName = nameSuffix.concat('_').concat('points_LPCOR');
//   Export.table.toDrive({
//     collection: points_LPCOR,
//     description: fileName,
//     folder: 'GEE_output',
//     fileFormat: 'SHP'
//   });
  
// }

// /////////////////////////////////////////////////////////////////////////////////////
// // \\\\\\\\\\\\\\\\\\\\\\\\\\\\ HRL IMPERVIOUSNESS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

var hrl30 = HR.gt(1).and(HR.lt(30)); // imperviousness < = 30% low density
var hrl30cl = ee.Image(2).updateMask(hrl30).toUint8();
var hrl100 = HR.gte(30); // imperviousness 30% - 100% dense
var hrl100cl = ee.Image(1).updateMask(hrl100).toUint8(); 

var HRcl = (ee.ImageCollection.fromImages([hrl30cl, hrl100cl])).mosaic().toUint8();

var cl_HR_L2_I = HRcl.rename('class_L2');
cl_HR_L2_I = cl_HR_L2_I.unmask();  // turn masked areas to 0, otherwise mask affects stdDev
Map.addLayer(cl_HR_L2_I, {}, 'HRL impervioussness', false)

// ---------- L2
// Compute per-cluster L2 class stdDev
var L2_segStdDev = clusters.addBands(cl_HR_L2_I).reduceConnectedComponents({
                                                  reducer: ee.Reducer.stdDev(), 
                                                  labelBand: 'clusters', 
                                                  maxSize: 100});

/*  !!!
You should specify the same neighborhood to both snic and reduceConnectedComponents, 
otherwise one of them might get truncated while the other one doesn't 
*/

// Compute per-cluster class
var L2_seg = cl_HR_L2_I.addBands(clusters).reduceConnectedComponents({
                                            reducer: ee.Reducer.mode(), 
                                            labelBand: 'clusters', 
                                            maxSize: 100});

// select only segments with one class (clear class segments)
// stdDev will be 0
var stdDevNot0 = L2_segStdDev.neq(0); // mask segments with stdDev <> 0
var L2_seg_cl = L2_seg.where(stdDevNot0, 0);
var gt28 = L2_seg_cl.gt(3); // mask segments with cl > number of classes
L2_seg_cl = L2_seg_cl.where(gt28, 0);
var zero = L2_seg_cl.eq(0); // mask segments with cl = 0
L2_seg_cl = L2_seg_cl.where(zero, 0);

// SAMPLE SEGMENTS FROM objectPropertiesImage ----------

// create a mask containing "clean" L2 segments only
var samp = L2_seg_cl.gt(0);

// mask clusters and keep only the ones over clean segments
var clusters_cl = clusters.where(samp.not(), 0);
clusters_cl = clusters.updateMask(samp);
var clusters_cl_class = L2_seg_cl.updateMask(samp).toUint8();
Map.addLayer(clusters_cl_class, {}, 'clean segments HRL', false)

// mask validation segments
var clusters_cl_class2 = clusters_cl_class.updateMask(val_seg.not());
// Map.addLayer(clusters_cl_class2, {}, 'clean segments Natura (removed validation)')

// mask area out of seeds, to only take seeds as samples
var clusters_cl_class3 = clusters_cl_class2.updateMask(seeds);
// Map.addLayer(clusters_cl_class3, {}, 'clean Natura seeds(removed validation)')

// cerate random points within clean segments
var points_HR = clusters_cl_class3.stratifiedSample({
                               numPoints: nPoints,
                               region: aoi,
                               scale: res,
                               projection: clusters_cl_class3.projection(),
                               tileScale: 16,
                               geometries: true
                               });
Map.addLayer(points_HR, {}, 'HRL points', false)

// // ---------- export points SHPs to Assets ----------
// var fileName = nameSuffix.concat('_').concat('points_HR');
// Export.table.toAsset({
//   collection: points_HR,
//   description: fileName,
//   assetId: fileName
// });

// // ---------- export points SHPs to Drive ----------
// if (exportFlag === true) {

//   var fileName = nameSuffix.concat('_').concat('points_HR');
//   Export.table.toDrive({
//     collection: points_HR,
//     description: fileName,
//     folder: 'GEE_output',
//     fileFormat: 'SHP'
//   });
  
// }

// ---------- MERGE AND EXPORT ALL ---------

var all = points_NA.merge(points_LPCOR).merge(points_HR);

// ---------- export points SHPs to Assets ----------
var fileName = nameSuffix.concat('_').concat('points_all_').concat(ee.String(ee.Number(nPoints)).getInfo());
Export.table.toAsset({
  collection: all,
  description: fileName,
  assetId: fileName
});

// // ---------- export points SHPs to Drive ----------
// if (exportFlag === true) {

//   var fileName = nameSuffix.concat('_').concat('points_all');
//   Export.table.toDrive({
//     collection: all,
//     description: fileName,
//     folder: 'GEE_output',
//     fileFormat: 'SHP'
//   });
  
// }

print('end');
