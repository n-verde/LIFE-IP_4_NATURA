/* =============================================================================

==== LIFE IP 4 NATURA - STEP 4/7 =====
Natalia Verde, AUTH/DUTH, 2020 

BRIEF DESCRIPTION:
This GEE script creates an image collection for classification using S2L2A, S2L1C and S1 imagery, DEM and object properties.
Creates an image collection with bands and indices for seasonal composites, and an image collection with indices for monthly
composites.
Mosaics of bands are already created and exported in assets.

NOTES:
- For national-scale processing must run in tiles, otherwise memory errors occur. Exports images for classification in tiles and 
  then reconstructs the tiles to single image for classifying.

HOW TO USE:
1. upload the EU-DEM (available from Copernicus) of your AOI to assets, import to script by changing the assets path for the "DEM" variable.
2. import the SNIC clusters to the script by changing the assets path for the "OBJECTS" variable (outputed to assets from previous step)
3. split your AOI to tiles, of size approximately 45000.000 km² each. Create a shapefile with the tiles and create a field
   containing the number of the tile.
   Upload the shapefile to assets. Import the tiles shapefile to the script by changing the assets path for the "grid" variable.
4. set the temporal mosaic interval as a asting ("month" or "season") (set the same as in script of step 3)
5. set if you want to remove band features (in order to create a reduced feature set)
6. set the property ("idName" variable) in your grid shapefile, which uniquely defines each tile.
7. set the number of the tile you are processing each time ("field_id" variable).
8. run for all tiles.
9. set other properties you want in the "PARAMETER SETTINGS" section.
10. in the tasks tab, click RUN to all tasks in order to export the mosaic images to Assets.

*/

//====================================================================================

//===============PARAMETER SETTINGS===================================================

// ---------- SET DEM ----------
// downloaded from Copernicus website and then uploaded to assets
var DEM = ee.Image('users/n-verde/shared/LIFE-IP_4_NATURA/eudem_v11_GR_EPSG4326');

// ---------- SET OBJECT CLUSTERS ----------
// computed in step 1
var OBJECTS = ee.Image('users/n-verde/shared/LIFE-IP_4_NATURA/LIFE_10m_OB_clusters');

// ---------- GRID
// created locally on PC and then uploaded
var grid = ee.FeatureCollection('users/n-verde/shared/LIFE-IP_4_NATURA/GR_boundingBoxGrid_LIFE');
//set unique polygon id attribute name in shapefile
var idName = 'id';
// choose one polygon
var field_id = 1; //--> choose one polygon, CHANGE EACH TIME YOU RUN
var AOI = grid.filter(ee.Filter.equals(idName, field_id));
// (Single polygon instead of tiles only works for areas < 45000.000 km²)

// ---------- SET TEMPORAL MOSAIC INTERVAL ----------
// var interval = 'month';
var interval = 'season';

// ---------- IMAGE EXPORT DETAILS ----------
// change the name as you change the tile, EACH TIME YOU RUN
// var nameSuffix = 'LIFE_OB_M_1'; // for monthly composites
var nameSuffix = 'LIFE_OB_S_1'; // for seasonal composites

// ---------- REMOVE FEATURES OPTION ----------
// remove features (bands) (for reduced feature sets = true)
var remFeat = true;

// ---------- SET MOSAIC PATHS FROM ASSETS ----------
// computed in step 3
var S2L2A = ee.Image('users/n-verde/GR_10m_S_MosaicS2L2A_masked_float');
var S2L1C = ee.Image('users/n-verde/GR_10m_S_MosaicS2L1C_masked_float');
var S1 = ee.Image('users/n-verde/GR_10m_S_MosaicS1_masked_float');

// ---------- SET EXTERNAL FUNCTIONS SCRIPTS ----------
var indicesS2L2A = require('users/n-verde/auth-shared:LIFE-IP_4_NATURA/4.1-IndicesS2L2A_function'); 
var indicesS2L1C = require('users/n-verde/auth-shared:LIFE-IP_4_NATURA/4.2-IndicesS2L1C_function'); 
var indicesS1 = require('users/n-verde/auth-shared:LIFE-IP_4_NATURA/4.3-IndicesS1_function');

var dem = require('users/n-verde/auth-shared:LIFE-IP_4_NATURA/4.4-DEM_function');

var OBIA = require('users/n-verde/auth-shared:LIFE-IP_4_NATURA/4.5-ObjectStatistics_function');

// ---------- OTHER SETTINGS ----------
var scale = 10;
var exp = false; // whether to export to drive or not (otherwise will only export to assets)
var MP = 1e13; // maxPixels variable (according to AOI and pixel size you want)


//====================================================================================

//--------------- FUNCTIONS --------------------------------

////////////////////////////////////////////////////////////
// This function transforms a multiband image of x bands into an image collection (or List of images)
function bandsToCollection(image){
  var collection = ee.ImageCollection.fromImages(image.bandNames().map(function(bandName){
    return image.select(ee.String(bandName));
  }));
  return collection;
}

////////////////////////////////////////////////////////////
// This function adds a property with season info for an image
// based on band name
function addSeason(image) { 

  var band = ee.String(image.bandNames().get(0));
  var s = ee.Number.parse(band.slice(1,2));
  image = image.set({season: s});
  
  return image;
}

////////////////////////////////////////////////////////////
// This function stacks an Image Collection by a property
function stackCollectionByPropertyOptical (col, propertyName) {
  
  // get images of idividual season as Collection
  var distinctCol = col.distinct([propertyName]);
  
  // filter for matching property
  var filter = ee.Filter.equals({leftField: propertyName, rightField: propertyName});
  var join = ee.Join.saveAll('matches'); // create join
  
  // function to iterate through properties. Used for mapping over distinct collection
  function stackBySeason(d) {
    
    var distinctIm = ee.ImageCollection.fromImages([d]);
    var resultsCol = join.apply(col, distinctIm, filter); // apply join
    var stackedIm = ee.ImageCollection(resultsCol).toBands();
  
    // rename the bands (remove the numbers created by .toBands())
    var oldNames = stackedIm.bandNames();
    var newNames = oldNames.map(function slicing(x) {
                                          var splitted = ee.String(x).split('_');
                                          var splitted2 = ee.String(splitted.get(1)).split('B');
                                          return ee.String('B').cat(ee.String(splitted2.get(-1)));
                                        });
    stackedIm = stackedIm.select(oldNames, ee.List(newNames));
    
    // return stacked image WITH PROPERTY
    return stackedIm.copyProperties(d).set(propertyName, d.get(propertyName));
  
  }

  var stackedCol = distinctCol.map(stackBySeason);
  
  return stackedCol;

}

////////////////////////////////////////////////////////////
// This function stacks an Image Collection by a property
function stackCollectionByPropertySAR (col, propertyName) {
  
  // get images of idividual season as Collection
  var distinctCol = col.distinct([propertyName]);
  
  // filter for matching property
  var filter = ee.Filter.equals({leftField: propertyName, rightField: propertyName});
  var join = ee.Join.saveAll('matches'); // create join
  
  // function to iterate through properties. Used for mapping over distinct collection
  function stackBySeason(d) {
    
    var distinctIm = ee.ImageCollection.fromImages([d]);
    var resultsCol = join.apply(col, distinctIm, filter); // apply join
    var stackedIm = ee.ImageCollection(resultsCol).toBands();
  
    // rename the bands (remove the numbers created by .toBands())
    var oldNames = stackedIm.bandNames();
    var newNames = oldNames.map(function slicing(x) {
                                          var splitted = ee.String(x).split('_');
                                          return ee.String(splitted.get(1)).slice(-2);
                                        });
    stackedIm = stackedIm.select(oldNames, ee.List(newNames));
    
    // return stacked image WITH PROPERTY
    return stackedIm.copyProperties(d).set(propertyName, d.get(propertyName));
  
  }

  var stackedCol = distinctCol.map(stackBySeason);
  
  return stackedCol;

}

/////////////////////////////////////////////////////////////
// This function removes useless stdDev per object features from the image for classification.
function removeUselessStdDevs(image) {
  
  // remove features from image
  // texture indices (stdDevs of objects)
  var initialTextureList = ee.List(['CorB2_1','PANTEX_1','CorBCI_1','CorVV_1','CorVH_1']);
  
  var season = ee.String('');
  
  if (interval === 'month') {
    // MONTHLY MOSAICS
    season = ee.String('S5');
    var S5textureList = initialTextureList.map(function(x) { return season.cat(x)});
    season = ee.String('S6');
    var S6textureList = initialTextureList.map(function(x) { return season.cat(x)});
    season = ee.String('S7');
    var S7textureList = initialTextureList.map(function(x) { return season.cat(x)});
    season = ee.String('S8');
    var S8textureList = initialTextureList.map(function(x) { return season.cat(x)});
    season = ee.String('S9');
    var S9textureList = initialTextureList.map(function(x) { return season.cat(x)});
    season = ee.String('S1');
    var S1textureList = initialTextureList.map(function(x) { return season.cat(x)});
    
    var allSeasontextureList_final = S5textureList.cat(S6textureList).cat(S7textureList)
                                    .cat(S8textureList).cat(S9textureList).cat(S1textureList); // stdDevs
  
  } else if (interval === 'season') {
    
    // SEASONAL MOSAICS
    season = ee.String('S2');
    var S2textureList = initialTextureList.map(function(x) { return season.cat(x)});
    season = ee.String('S3');
    var S3textureList = initialTextureList.map(function(x) { return season.cat(x)});
    season = ee.String('S4');
    var S4textureList = initialTextureList.map(function(x) { return season.cat(x)});
    
    var allSeasontextureList_final = S2textureList.cat(S3textureList).cat(S4textureList); // stdDevs
    
  }
  
  // DEM features (stdDevs of objects)
  var demList_final = ee.List(['elev_1','slope_1','TRASP_1']);
  
  var allStdDevFeatures = allSeasontextureList_final.cat(demList_final);
  
  var imageWithoutStdDevs = image.select(image.bandNames().removeAll(allStdDevFeatures));
  
  return imageWithoutStdDevs;
}

function removeFeatures(image) {
  
  // remove features from image
  // bands
  var initialFeatureList = ee.List(['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12',
                                    'B2_1','B3_1','B4_1','B5_1','B6_1','B7_1','B8_1','B8A_1','B11_1','B12_1']);
  
  var season = ee.String('');
  
  
  if (interval === 'month') {
  
  // MONTHLY MOSAICS
  season = ee.String('S5');
  var S5featureList = initialFeatureList.map(function(x) { return season.cat(x)});
  season = ee.String('S6');
  var S6featureList = initialFeatureList.map(function(x) { return season.cat(x)});
  season = ee.String('S7');
  var S7featureList = initialFeatureList.map(function(x) { return season.cat(x)});
  season = ee.String('S8');
  var S8featureList = initialFeatureList.map(function(x) { return season.cat(x)});
  season = ee.String('S9');
  var S9featureList = initialFeatureList.map(function(x) { return season.cat(x)});
  season = ee.String('S1');
  var S1featureList = initialFeatureList.map(function(x) { return season.cat(x)});
  
  var allSeasonfeatureList = S5featureList.cat(S6featureList).cat(S7featureList)
                                  .cat(S8featureList).cat(S9featureList).cat(S1featureList); // bands
  } else if (interval === 'season') {
    
  // SEASONAL MOSAICS
  season = ee.String('S2');
  var S2featureList = initialFeatureList.map(function(x) { return season.cat(x)});
  season = ee.String('S3');
  var S3featureList = initialFeatureList.map(function(x) { return season.cat(x)});
  season = ee.String('S4');
  var S4featureList = initialFeatureList.map(function(x) { return season.cat(x)});
  
  var allSeasonfeatureList = S2featureList.cat(S3featureList).cat(S4featureList); // bands
  
  }
    
  var imageWithoutFeatures = image.select(image.bandNames().removeAll(allSeasonfeatureList));

  return imageWithoutFeatures;
}

//====================================================================================

//---------------MAIN PROGRAM--------------------------------

print('start');

// ---------- RECONSTRUCT IMAGE COLLECTIONS FROM IMAGES----------

var S2L2A_IC = bandsToCollection(S2L2A);
var S2L1C_IC = bandsToCollection(S2L1C);
var S1_IC = bandsToCollection(S1);

// Add month property to images based on band name
S2L2A_IC = S2L2A_IC.map(addSeason);
S2L1C_IC = S2L1C_IC.map(addSeason);
S1_IC = S1_IC.map(addSeason);

// Stack images by season
var propertyName = ee.String('season');

var S2L2A_IC_s = stackCollectionByPropertyOptical(S2L2A_IC, propertyName); // print('S2L2A_IC_s',S2L2A_IC_s)
var S2L1C_IC_s = stackCollectionByPropertyOptical(S2L1C_IC, propertyName); // print('S2L1C_IC_s',S2L1C_IC_s)
var S1_IC_s = stackCollectionByPropertySAR(S1_IC, propertyName); // print('S1_IC_s',S1_IC_s)

// ---------- ADD INDICES TO IMAGE COLLECTIONS ----------

S2L2A_IC_s = indicesS2L2A.IndicesS2L2A(S2L2A_IC_s,AOI);
S2L1C_IC_s = indicesS2L1C.IndicesS2L1C(S2L1C_IC_s,AOI);
S1_IC_s = indicesS1.IndicesS1(S1_IC_s,AOI);

// ---------- TURN TO IMAGES ----------
// S2L2A ----------
var S2L2A_I_s = S2L2A_IC_s.toBands(); // turn collection to ee.Image

// rename the bands (remove the numbers created by .toBands())
var oldNames = S2L2A_I_s.bandNames();
var newNames = oldNames.map(function slicing(x) {
                                      var splitted = ee.String(x).split('_');
                                      return ee.String(splitted.get(1));
                                    });
var newNamesClient = newNames.map(function(x) {
                              var e = ee.String('S').cat(x); return e;
                              }); 
S2L2A_I_s = S2L2A_I_s.select(oldNames, ee.List(newNamesClient));
// print('S2L2A_I_s', S2L2A_I_s)

var IC_ImageCat = S2L2A_I_s;

// S2L1C ----------
var S2L1C_I_s = S2L1C_IC_s.toBands(); // turn collection to ee.Image

// rename the bands (remove the numbers created by .toBands())
var oldNames = S2L1C_I_s.bandNames();
var newNames = oldNames.map(function slicing(x) {
                                      var splitted = ee.String(x).split('_');
                                      return ee.String(splitted.get(1));
                                    });
var newNamesClient = newNames.map(function(x) {
                              var e = ee.String('S').cat(x); return e;
                              }); 
S2L1C_I_s = S2L1C_I_s.select(oldNames, ee.List(newNamesClient));
// print('S2L1C_I_s', S2L1C_I_s)

IC_ImageCat = ee.Image.cat(IC_ImageCat, S2L1C_I_s);

// S1 ----------
var S1_I_s = S1_IC_s.toBands(); // turn collection to ee.Image

// rename the bands (remove the numbers created by .toBands())
var oldNames = S1_I_s.bandNames();
var newNames = oldNames.map(function slicing(x) {
                                      var splitted = ee.String(x).split('_');
                                      return ee.String(splitted.get(1));
                                    });
var newNamesClient = newNames.map(function(x) {
                              var e = ee.String('S').cat(x); return e;
                              }); 
S1_I_s = S1_I_s.select(oldNames, ee.List(newNamesClient));
// print('S1_I_s', S1_I_s)

IC_ImageCat = ee.Image.cat(IC_ImageCat, S1_I_s);

// Add DEM ---->
var ICfinal= dem.addDEM(IC_ImageCat,DEM,AOI); 
print('final image with indices', ICfinal)

// ---------- OBJECT STATISTICS ----------
  
var segOutputs = OBIA.objStats(OBJECTS, AOI, ICfinal);

var clusters = segOutputs.clusters; // Map.addLayer(clusters)
var objectPropertiesImage = segOutputs.objectPropertiesImage;
var ICfinalOBIA = objectPropertiesImage;

var IC4classif = ICfinalOBIA;

// remove useless stdDevs
IC4classif = removeUselessStdDevs(IC4classif);

// remove features (bands) (for reduced feature sets)
if (remFeat === true) {
  IC4classif = removeFeatures(IC4classif);
}

print('final image with objects statistics:',IC4classif) 
// bands with "_1" at the end are "mean"
// bands with original names are "stdDev"

// ==========================================================================

// --------------------------- DISPLAYS ------------------------
Map.addLayer(AOI, {}, 'AOI', false)
Map.addLayer(IC4classif, {}, 'IC4classif', false)

// ==========================================================================

// --------------------------- EXPORTS ------------------------

var firstBand = ee.String(IC4classif.bandNames().get(1));
var exportProjection = IC4classif.select(firstBand).projection();

// TO ASSETS ----------

// Export classification image (if obia = TRUE: objects with stdDev & Mean)
// Name for your file export
var descr = nameSuffix + '_IC4classif';
Export.image.toAsset({
  image:IC4classif,
  description: descr,
  scale: scale,
  crs: exportProjection,
  region: AOI,
  maxPixels: MP, //max pixels allowed for download
  pyramidingPolicy: {".default": "mean"}
});

// TO DRIVE ----------

if (exp === true) {
  
  //  ~~~~~~~~~~ IMAGES  ~~~~~~~~~~

  // make a list of the image dates from the original collection to export as csv
  var names = ee.List(ICfinal.bandNames());
  
  var namesfc = ee.FeatureCollection(names.map(function(nam){
      return ee.Feature(null, {name:nam});
    }));
    
  // Export band names
  
  // make var "nameSuffix" an GEE object from a client-side string
  nameSuffix = ee.String(nameSuffix).getInfo();
  
  // Name for description in Tasks
  var desc = nameSuffix + '_bandNames_' ;
  //Name for your file export
  var fileName = desc; 
      
  Export.table.toDrive({
    collection: namesfc,
    description: fileName,
    folder: 'GEE_Output',
    fileFormat: 'CSV'
  });
  
  // Export stack image
  //Name for your file export
  var descr = nameSuffix + '_stack';
  Export.image.toDrive({
    image:ICfinal,
    description: descr,
    folder: 'GEE_output',
    scale: scale,
    crs: exportProjection,
    region: AOI,
    maxPixels: MP //max pixels allowed for download
  });
  
  // export perimeters for testing reasons
    var descr = nameSuffix + 'perimeters';
    Export.image.toDrive({
      image:perimeterPixels,
      description: descr,
      folder: 'GEE_output',
      scale: scale,
      crs: perimeterPixels.projection(),
      region: AOI,
      maxPixels: MP //max pixels allowed for download
    });
    
  // Export classification image (if obia = TRUE: objects with stdDev & Mean)
  // Name for your file export
  var descr = nameSuffix + '_IC4classif';
  Export.image.toDrive({
    image:IC4classif,
    description: descr,
    folder: 'GEE_output',
    scale: scale,
    crs: exportProjection,
    region: AOI,
    maxPixels: MP //max pixels allowed for download
  });
  

}

print('end');


