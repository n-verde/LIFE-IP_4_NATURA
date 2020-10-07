/* =============================================================================

==== LIFE-IP 4 NATURA - STEP 1/7 =====
Natalia Verde, AUTH/DUTH, 2020

BRIEF DESCRIPTION:
This GEE script performs a SNIC segmentation to the july image of an image collection and
exports the segment perimeters and seeds as a raster to Drive and Assets.

NOTES:
- 

HOW TO USE:
1. upload your AOI shapefile to assets, and set the "AOI" variable with the path in your assets.
2. set other properties you want in the "PARAMETER SETTINGS" section.
3. in the tasks tab, click RUN to all tasks in order to export the segmentation images to Assets.

*/

//====================================================================================
 
//===============PARAMETER SETTINGS===================================================

// ---------- SET AOI ----------
// created locally on PC and then uploaded
var AOI = ee.FeatureCollection('users/n-verde/shared/LIFE-IP_4_NATURA/GR_boundingBox_buff1km')

// ---------- SET DATE RANGE FOR TRAINING IMAGE SEARCH ----------
// Date range (july month)
var startDate1 = ee.Date.fromYMD(2019,7,1); 
var endDate1 = ee.Date.fromYMD(2019,7,15);
var startDate2 = ee.Date.fromYMD(2019,7,15); 
var endDate2 = ee.Date.fromYMD(2019,7,30);

// ---------- SET CLOUD FREE IMAGE PERCENTAGE ----------
//(from 0 to 1) (recomended 1)
var NonCloudPerc = 0.9; //first perc of images in image coll., sorted by cloud coverage

// ---------- IMAGE EXPORT DETAILS ----------
var nameSuffix = 'LIFE_10m_OB_';

// ---------- OTHER SETTINGS ----------
//(it's not recomended to change these values)
var scale = 10;
var exp = false; // whether to export to drive or not (otherwise will only export to assets)
var obia = true; // whether to perform object-based image analysis
var MP = 1e13; // maxPixels variable (according to AOI and pixel size you want)
var seedGridSize = 10; // the size of the grid for superpixels (the biger, the bigger the segments)

//====================================================================================

//---------------FUNCTIONS------------------------------

function createCollection(startDate1,endDate1,startDate2,endDate2,aoi) {

  //--------------- FUNCTIONS ------------------------------
  
  /////////////////////////////////////////////////////////////
  /* 
   * This function searches a S2 collection according to date Range
   * and cloud coverage
   * variable "Icol" is the name of image collection (string)
   * "bounds" is the shp file
   */
  function filterCol(startdate, enddate, bounds) {
    //Asign proper cloud and collection filter
    var cloudfilter = 'CLOUDY_PIXEL_PERCENTAGE';
    var Icol = 'COPERNICUS/S2_SR';
  
    //Search Image Collection
    var IC = ee.ImageCollection(Icol)
      .filterBounds(bounds.geometry()) //reduce search to bounds
      .filterDate(startdate, enddate); //reduce search to date range
   
    // split image colleciton into less and more cloudy
    IC = IC
      .sort(ee.String(cloudfilter))
      .limit(IC.size().multiply(NonCloudPerc).int());
    
    IC = IC
      .sort('system:time_start');
  
    //Number of Images in Collection
    var Icount = IC.size();
    // print(Icol, ' collection size: ', Icount)
    
    //return a collection
    return ee.ImageCollection(IC)
            .select(['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11','B12',
            'QA60','SCL','MSK_CLDPRB']);
  }
  
  /////////////////////////////////////////////////////////////
  // This function masks clouds from a Sentinel-2 L2A image
  var maskCloudsS2L2A = function(image){
   
    var cloudProb = image.select('MSK_CLDPRB');
    var scl = image.select('SCL'); 
  
    var shadow = scl.eq(3); // 3 = cloud shadow
    var cirrus = scl.eq(10); // 10 = cirrus
    var mask = cloudProb.lt(5).and((cirrus).neq(1)).and((shadow).neq(1));
    
    return image.select(['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11','B12'])
                .updateMask(mask);
  
  };
  
  ////////////////////////////////////////////////////////////
  // This function adds a property with Month info for an image
  function addMonthForName(image) {
      image = ee.Image(image); //turn to ee.Image object so google knows what type it is
      //get month from image and turn it from string to number
      var month = ee.Number.parse(image.date().format('MM'));
    return image.set({'Month': month});
  }
  
  function addMonthForMos(image) { 
      image = ee.Image(image);
    return image.set({month: image.date().format('MM')});
  }
  
  ////////////////////////////////////////////////////////////
  // This function adds a property with season info for an image 
  // (winter, spring, summer, autumn)
  // function addSeasonForMos(image) { 
  //     image = ee.Image(image);
  //     if image.month
  //   return image.set({season: image.date().format('MM')});
  // }
  
  ////////////////////////////////////////////////////////////
  // This function mosaics a collection by month based on median values
  function mosaicByMonthTR(images) {
      var property = 'month';
      
      var distinct = images.distinct([property]); //Removes duplicates from a collection
  
      // ---------- Create a time filter to define a match as overlapping timestamps.
      var filter = ee.Filter.equals({leftField: property, rightField: property});
      
      // ---------- Define the join.
      var saveAllJoin = ee.Join.saveAll({
        matchesKey: 'matches',
        ordering: 'system:time_start',
        ascending: true
      });
      
      // ---------- Apply the join.
      var results = saveAllJoin.apply(distinct, images, filter);
      
      // mosaic
      var bandNames = ee.Image(images.first()).bandNames();
      results = results.map(function(i) {
        var mosaic = null;
        mosaic = ee.ImageCollection.fromImages(i.get('matches')).sort('Month')
                  .median().float();  // create a median mosaic
        
        return mosaic.copyProperties(i).set(property, i.get(property))
        .set('system:time_start', ee.Date(i.get(property)).millis());
      });
    return ee.ImageCollection(results);
  }
  
  ////////////////////////////////////////////////////////////
  // This function creates a mask outside the aoi
  function maskOut(image) {
      var mask = ee.Image.constant(0).float();
      mask = mask.paint(aoi,1);
      var newimage = image.updateMask(mask);
    return newimage;
  }
  
  //--------------- MAIN ALGORITHM --------------------------------

  // ---------- CREATE IMAGE COLLECTION ----------
  
  // filter by date range, cloud coverage and area of interest
  var Col1 = filterCol(startDate1, endDate1, aoi);
  var Col2 = filterCol(startDate2, endDate2, aoi);
  var Col = Col1.merge(Col2);
  
  // get the collection ID
  // usefull for  indices calculation later
  var colID = (ee.String((Col.first()).get("system:id")).replace("/[^/]*$", "")).getInfo();
  // print('found image collection for date range:', colID)
  
  // add month property to create monthly mosaics (ONLY MONTH, NO YEAR)
  Col = Col.map(addMonthForMos);
  // print('add month property to collection:', Col)
  
  // add month property to name bands (ONLY MONTH, NO YEAR)
  Col = Col.map(addMonthForName);
  // print('added months (name and for mos.):', Col)
  
  //Mask clouds
  Col = Col.map(maskCloudsS2L2A);
  // Map.addLayer(Col.first(),{},'MaskedCLouds')
  
  var ColMos = mosaicByMonthTR(Col); //create image collection with mosaics by month
  
  var imageVisParam = {"opacity":1,"bands":["B4","B3","B2"],"min":-111.92345461786806,"max":1994.995660354762,"gamma":1};
  Map.addLayer(ColMos,imageVisParam, 'ColMos', true)
  
  return ColMos;

}

function segmentation(segImage, ICfinal, seedgrid) {
  
  // ---------- SEGMENTATION ----------
  
  // Nick Clinton's code (modified by NV) ---------->
  
  // create superpixels grid: how many pixels away will the next superpixel be?
  var seeds = ee.Algorithms.Image.Segmentation.seedGrid(seedgrid); // superpixel every 10 pixels (100m)
  // Map.addLayer('seeds',seeds)
  
  /*  !!!
  You should specify the same neighborhood to both snic and reduceConnectedComponents, 
  otherwise one of them might get truncated while the other one doesn't 
  */
  
  // Run SNIC on the regular square grid, on img2348
  var clusters = ee.Algorithms.Image.Segmentation.SNIC({
    image: segImage, 
    // size: 10, // The superpixel seed location spacing, in pixels.
    compactness: 3,
    connectivity: 8,
    neighborhoodSize: 100,
    seeds: seeds
  }).select(['clusters']); // select only band with clusters (segments)
  // Map.addLayer(clusters.randomVisualizer(), {}, 'clusters')
  // print('clusters', clusters)
  
  // perimeter of clusters
  var minMax = clusters.reduceNeighborhood(ee.Reducer.minMax(), ee.Kernel.square(1));
  var perimeterPixels = minMax.select(0).neq(minMax.select(1)).rename('perimeter');
  
  // Compute per-cluster stdDev
  var stdDev = ICfinal.addBands(clusters).reduceConnectedComponents(ee.Reducer.stdDev(), 'clusters', 100);
  // Map.addLayer(stdDev, {min:80, max:400}, 'per cluster StdDev', false)
  
  // Compute per-cluster Mean
  var means = ICfinal.addBands(clusters).reduceConnectedComponents(ee.Reducer.mean(), 'clusters', 100);
  // Map.addLayer(stdDev, {min:80, max:400}, 'per cluster mean', false)
  // print('per cluster mean:', means)
  
  var objectPropertiesImage = ee.Image.cat([
    stdDev,
    means
  ]).float();
  // bands with _1 in their name are the MEANs, bands with origina names are STDDEVs
  // print('objectPropertiesImage:', objectPropertiesImage)
  // Map.addLayer(objectPropertiesImage, {}, 'objectPropertiesImage')
  
  // image for sampling and classification becomes objectPropertiesImage
  return {seeds:seeds, perimeterPixels:perimeterPixels, clusters:clusters, objectPropertiesImage:objectPropertiesImage};
  
}


//====================================================================================

//---------------MAIN PROGRAM--------------------------------

print('start');

// ---------- CREATE IMAGE COLLECTIONS ----------
var IC = createCollection(startDate1, endDate1, startDate2, endDate2, AOI);
// print(IC)
var IC_Image = IC.toBands(); // turn collection to ee.Image

var ICfinal = IC_Image;

// ---------- SEGMENTATION ----------
  
// select 10m bands of July month for segmentation
var julImg = IC.filterMetadata('month','equals','07').first();
// print('make sure this is a july image:', julImg)
  
var img2348 = julImg.select('B2','B3','B4','B8');
var segOutputs = segmentation(img2348, ICfinal, seedGridSize);

var seeds = segOutputs.seeds;
var clusters = segOutputs.clusters;
var perimeterPixels = segOutputs.perimeterPixels;


// ==========================================================================

// --------------------------- DISPLAYS ------------------------
Map.addLayer(AOI, {}, 'AOI', true)
Map.addLayer(perimeterPixels, {min: 0, max: 1}, 'perimeterPixels', false)
Map.addLayer(clusters, {}, 'clusters', false)

// ==========================================================================

// --------------------------- EXPORTS ------------------------

// TO ASSETS ----------

// export seeds
  var descr = nameSuffix + 'seeds';
  Export.image.toAsset({
    image:seeds,
    description: descr,
    scale: scale,
    crs: seeds.projection(),
    region: AOI,
    maxPixels: MP //max pixels allowed for download
  });
  
  // export clusters
  var descr = nameSuffix + 'clusters';
  Export.image.toAsset({
    image:clusters,
    description: descr,
    scale: scale,
    crs: clusters.projection(),
    region: AOI,
    maxPixels: MP //max pixels allowed for download
  });
  
  // export perimeters
  var descr = nameSuffix + 'perimeters';
  Export.image.toAsset({
    image:perimeterPixels,
    description: descr,
    scale: scale,
    crs: perimeterPixels.projection(),
    region: AOI,
    maxPixels: MP //max pixels allowed for download
  });

// TO DRIVE ----------

if (exp === true) {
  
  // export seeds
  var descr = nameSuffix + 'seeds';
  Export.image.toDrive({
    image:seeds,
    description: descr,
    folder: 'GEE_output',
    scale: scale,
    crs: seeds.projection(),
    region: AOI,
    maxPixels: MP //max pixels allowed for download
  });
  
  // export clusters
  var descr = nameSuffix + 'clusters';
  Export.image.toDrive({
    image:clusters,
    description: descr,
    folder: 'GEE_output',
    scale: scale,
    crs: clusters.projection(),
    region: AOI,
    maxPixels: MP //max pixels allowed for download
  });
  
  // export perimeters
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

}

print('end');
