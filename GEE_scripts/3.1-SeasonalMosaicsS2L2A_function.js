/*=============================================================================

==== LIFE IP 4 NATURA ===== 
Natalia Verde, AUTH/DUTH, 2020

BRIEF DESCRIPTION:
This GEE function creates seasonal mosaics of S2L2A images for a specific date range
for an area of interest (AOI) to later use in creating indices.

Bands:
-B2 
-B3
-B4
-B5
-B6
-B7
-B8
-B8A
-B11
-B12

 
*/

//====================================================================================

//===============PARAMETER SETTINGS===================================================

// ---------- SET CLOUD FREE IMAGE PERCENTAGE ----------
//(from 0 to 1) (recomended 1)
var NonCloudPerc = 0.9; //first perc of images in image coll., sorted by cloud coverage

//====================================================================================

//--------------- MAIN FUNCTION --------------------------

function seasonalMosS2L2A(startDate1,endDate1,startDate2,endDate2,aoi,interval) {

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
  // This function adds a property with month and season info for an image
  function addMonthSeason(image) { 
    image = ee.Image(image);
    
    // set month
    var month = ee.Number.parse(image.date().format('M'));
    image = image.set({month: month});
    
    // set season (1=winter, 2=spring, 3=summer, 4=autumn)
    var s = ee.Algorithms.If(month.eq(ee.Number(1)),ee.Number(1),
            ee.Algorithms.If(month.eq(ee.Number(2)),ee.Number(1),
            ee.Algorithms.If(month.eq(ee.Number(3)),ee.Number(2),
            ee.Algorithms.If(month.eq(ee.Number(4)),ee.Number(2),
            ee.Algorithms.If(month.eq(ee.Number(5)),ee.Number(2),
            ee.Algorithms.If(month.eq(ee.Number(6)),ee.Number(3),
            ee.Algorithms.If(month.eq(ee.Number(7)),ee.Number(3),
            ee.Algorithms.If(month.eq(ee.Number(8)),ee.Number(3),
            ee.Algorithms.If(month.eq(ee.Number(9)),ee.Number(4),
            ee.Algorithms.If(month.eq(ee.Number(10)),ee.Number(4),
            ee.Algorithms.If(month.eq(ee.Number(11)),ee.Number(4),
            ee.Algorithms.If(month.eq(ee.Number(12)),ee.Number(1)))))))))))));
    
    image = image.set({season: s});
    
    return image;
  }
  
  ////////////////////////////////////////////////////////////
  // This function mosaics a S2 collection by season based on median values
  function mosaicBySeasonS2(images) {
      var property = interval;
      
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
        mosaic = ee.ImageCollection.fromImages(i.get('matches')).sort(interval)
                  .median().uint16();  // create a median mosaic
        
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
  
  //====================================================================================
  
  //--------------- MAIN ALGORITHM --------------------------------

  // ---------- CREATE IMAGE COLLECTION ----------
  
  // ---------- S2 ----------
  
  // filter S2 by date range, cloud coverage and area of interest
  var ColS21 = filterCol(startDate1, endDate1, aoi);
  var ColS22 = filterCol(startDate2, endDate2, aoi);
  var ColS2 = ColS21.merge(ColS22);
  print(ColS2)
  
  // get the collection ID
  // usefull for  indices calculation later
  var colID = (ee.String((ColS2.first()).get("system:id")).replace("/[^/]*$", "")).getInfo();
  // print('found image collection for date range:', colID)
  
  // add month and season property
  ColS2 = ColS2.map(addMonthSeason);
  // print('add month and season property to S2 collection:', ColS2)
  
  //Mask clouds
  ColS2 = ColS2.map(maskCloudsS2L2A);
  // Map.addLayer(ColS2.first(),{},'MaskedCLouds')
  
  var ColMosS2 = mosaicBySeasonS2(ColS2); //create image collection with mosaics by season
  // print('number of seasonal mosaics for training (S2):', ColMosS2.size())
  // print('mosaic by season S2 IC:', ColMosS2)
  
  // mask outside AOI
  ColMosS2 = ColMosS2.map(maskOut);
  // Map.addLayer(ColMosS2.first(),{},'S2 Mosaic of first season')
  
  var bands_listS2 = ['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12'];  

  // select only the specific bands from the mosaic collection
  ColMosS2 = ColMosS2.select(bands_listS2); 
  
  // first rename all the bands by date
  function renameBands(image){
    var int = ee.String(image.get(interval));
    var new_bands = [ (int.cat('B2')),     // 1
                      (int.cat('B3')),     // 2
                      (int.cat('B4')),     // 3
                      (int.cat('B5')),     // 4
                      (int.cat('B6')),     // 5
                      (int.cat('B7')),     // 6
                      (int.cat('B8')),     // 7
                      (int.cat('B8A')),    // 8
                      (int.cat('B11')),    // 9
                      (int.cat('B12'))];   // 10
    return image.select(bands_listS2,new_bands);
  }
  ColMosS2 = ColMosS2.map(renameBands);
  // print('rename bands:', ColMosS2)
  
  // set images ids to null, so ".toBands()" will not assign image id to band names
  ColMosS2 = ColMosS2.map(function x(im){
                        im = im.set({'system:index':''});
                        return im;
                      });
  
return ColMosS2;

}

exports.seasonalMosS2L2A = seasonalMosS2L2A;
