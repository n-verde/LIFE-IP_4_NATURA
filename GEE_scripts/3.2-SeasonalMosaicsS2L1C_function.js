/*=============================================================================

==== LIFE IP 4 NATURA ===== 
Natalia Verde, AUTH/DUTH, 2020

BRIEF DESCRIPTION:
This GEE function creates seasonal mosaics of Sentinel-2 L1C data. Used for later
creation of indices.

*/

//====================================================================================
  
// ---------- SET CLOUD FREE IMAGE PERCENTAGE ----------
//(from 0 to 1) (recomended 1)
var NonCloudPerc = 0.9; //first perc of images in image coll., sorted by cloud coverage
 

function seasonalMosS2L1C(startDateB, endDateB, aoi, interval){

  //====================================================================================
  
  //---------------FUNCTIONS------------------------------
  
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
    var Icol = 'COPERNICUS/S2'; // searches in S2 L1C because there are no coefficients for L2A data!
  
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
              .select(['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9', 'B10','B11','B12','QA60']);
  }                   //  0,   1,   2,   3,   4,   5,   6,   7,    8,   9,    10,   11,   12,   13
  
  /////////////////////////////////////////////////////////////
  // This function sets the correct parametes in a S2 image, to later apply cloud and shadow
  // masking functions
  function sentinel2toa(img) {
    var toa = img .addBands(img.select(['QA60']))
                  .set('solar_azimuth',img.get('MEAN_SOLAR_AZIMUTH_ANGLE'))
                  .set('solar_zenith',img.get('MEAN_SOLAR_ZENITH_ANGLE'));
      return toa;
  }
  /////////////////////////////////////////////////////////////
  // This function finds clouds in S2 image
  function ESAcloud(toa) {
    // author: Nick Clinton
    
    var qa = toa.select('QA60');
    
    // Bits 10 and 11 are clouds and cirrus, respectively.
    var cloudBitMask = Math.pow(2, 10);
    var cirrusBitMask = Math.pow(2, 11);
    
    // clear if both flags set to zero.
    var clear = qa.bitwiseAnd(cloudBitMask).eq(0).and(
               qa.bitwiseAnd(cirrusBitMask).eq(0));
    
    var cloud = clear.eq(0);
    return cloud;
  }
  /////////////////////////////////////////////////////////////
  // This function finds shadows from clouds in S2 image
  function shadowMask(toa,cloud){
    // Author: Gennadii Donchyts
    // License: Apache 2.0
    
    // solar geometry (radians)
    var azimuth =ee.Number(toa.get('solar_azimuth')).multiply(Math.PI).divide(180.0)
                          .add(ee.Number(0.5).multiply(Math.PI));
    var zenith  =ee.Number(0.5).multiply(Math.PI ).subtract(ee.Number(toa.get('solar_zenith'))
                          .multiply(Math.PI).divide(180.0));
  
    // find where cloud shadows should be based on solar geometry
    var nominalScale = cloud.projection().nominalScale();
    var cloudHeights = ee.List.sequence(200,10000,500);
    var shadows = cloudHeights.map(function(cloudHeight){
      cloudHeight = ee.Number(cloudHeight);
      var shadowVector = zenith.tan().multiply(cloudHeight);
      var x = azimuth.cos().multiply(shadowVector).divide(nominalScale).round();
      var y = azimuth.sin().multiply(shadowVector).divide(nominalScale).round();
      return cloud.changeProj(cloud.projection(), cloud.projection().translate(x, y));
    });
    var potentialShadow = ee.ImageCollection.fromImages(shadows).max();
    
    // shadows are not clouds
    potentialShadow = potentialShadow.and(cloud.not());
    
    // (modified by Sam Murphy) dark pixel detection 
    var darkPixels = toa.normalizedDifference(['B3','B12']).gt(0.25).rename(['dark_pixels']);
    
    // shadows are dark
    var shadow = potentialShadow.and(darkPixels).rename('shadows');
    return shadow;
  }
  /////////////////////////////////////////////////////////////
  // This function, combines the previous functions
  // returns a S2 image with mask at clouds and shadows
  function cloud_and_shadow_maskS2(img) {
    var toa = sentinel2toa(img);
    var cloud = ESAcloud(toa);
    var shadow = shadowMask(toa,cloud);
    var mask = cloud.or(shadow).eq(0);
    return toa.updateMask(mask);
  }
  
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
  // This function mosaics a collection by season based on median values
  function mosaicBySeason(images) {
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
  
  
  //--------------- MAIN FUNCTION --------------------------

  // filter by date range, cloud coverage and area of interest
  var Col = filterCol(startDateB, endDateB, aoi);
  print(Col)
  
  // get the collection ID
  // usefull for  indices calculation later
  var colID = (ee.String((Col.first()).get("system:id")).replace("/[^/]*$", "")).getInfo();
  // print('found image collection for date range:', colID)
  
  // add month and season property
  Col = Col.map(addMonthSeason);
  // print('add month and season property to collection:', Col)
  
  //Mask clouds and shadows from them
  Col = Col.map(cloud_and_shadow_maskS2);
  // print('masked clouds collection:', Col)
  // Map.addLayer(Col.first(), {}, 'Masked Collection')
  
  var ColMos = mosaicBySeason(Col); //create image collection with mosaics by season
  // print('number of seasonal mosaics:', ColMos.size())
  // print('mosaic by season ic:', ColMos)
  
  // mask outside AOI
  ColMos = ColMos.map(maskOut);
  
  // bands to rename
  var bands_array = ['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9', 'B10','B11','B12'];
  
  // select only the specific bands from the mosaic collection
  ColMos = ColMos.select(bands_array);
  
  // first rename all the bands by date
  function renameBands(image){
    var int = ee.String(image.get(interval));  
    var new_bands = [ (int.cat('B1')),         // 1
                      (int.cat('B2')),         // 2
                      (int.cat('B3')),         // 3
                      (int.cat('B4')),         // 4
                      (int.cat('B5')),         // 5
                      (int.cat('B6')),         // 6
                      (int.cat('B7')),         // 7
                      (int.cat('B8')),         // 8
                      (int.cat('B8A')),        // 9
                      (int.cat('B9')),         // 10
                      (int.cat('B10')),        // 11
                      (int.cat('B11')),        // 12
                      (int.cat('B12')),        // 13
                      ];    
    return image.select(bands_array,new_bands);
  }
  ColMos = ColMos.map(renameBands);
  
  // returns a collection
  return ColMos;
  
}

exports.seasonalMosS2L1C = seasonalMosS2L1C;