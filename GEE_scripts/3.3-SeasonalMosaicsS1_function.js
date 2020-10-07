/*=============================================================================

==== LIFE IP 4 NATURA ===== 
Natalia Verde, AUTH/DUTH, 2020

BRIEF DESCRIPTION:
This GEE function creates seasonal mosaics of Sentinel-1 IW GRD images based on a date range.
for an area of interest (AOI) to later use for indices creation.

Bands:
-VV
-VH
 
*/

//====================================================================================

//--------------- MAIN FUNCTION --------------------------

function seasonalMosS1(startDate1,endDate1,startDate2,endDate2,aoi,interval) {

  //--------------- FUNCTIONS ------------------------------
  
  ////////////////////////////////////////////////////////////
  // This function removes edge pixels in S1 using the connected Componets
  function maskEdge(img) {
          return img.clip(img.geometry().buffer(-2000));
  }
  
  /////////////////////////////////////////////////////////////
  /* 
   * This function searches a S1 collection according to date Range
   * returns IW mode VV and VH polarization
   * variable "Icol" is the name of image collection (string)
   * "bounds" is the shp file
   */
  function filterColS1(startdate, enddate, bounds) {
    //Asign proper cloud and collection filter
    var Icol = 'COPERNICUS/S1_GRD_FLOAT';
  
    //Search Image Collection
    var IC = ee.ImageCollection(Icol)
      .filterBounds(bounds.geometry()) //reduce search to bounds
      .filterDate(startdate, enddate) //reduce search to date range
      .filterMetadata('instrumentMode', 'equals', 'IW')
      .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
      .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'));
    
    IC = IC.map(maskEdge);  // remove black edges by buffering
    
    IC = IC
      .sort('system:time_start');
  
    //Number of Images in Collection
    var Icount = IC.size();
    // print(Icol, ' collection size: ', Icount)
    
    //return a collection
    return ee.ImageCollection(IC);
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
  // This function mosaics a S1 collection by season based on mean values
  function mosaicBySeasonS1(images) {
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
                  .mean();  // create a mean mosaic
        
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
  
    /////////////////////////////////////////////////////////////
  // This function  takes an image and then scales the output to [0,1].
  // ex. normalize(img, 'NDVI', img.geometry(), imgMin, imgMax)
  // ATTENTION! THIS IS RESCALE, NOT STRETCH!!!
  function normalize(img, bandName, min, max) {
    
    var Min = ee.Number(min);
    var Max = ee.Number(max);
    
    var imgBand = img.select(bandName);
      
    if (Min.lt(0)) {
      imgBand = imgBand.add(Min.abs()); // turn negative values to possitive
      Max = Max.add(Min.abs()); // change new max 
    }
      
    imgBand = imgBand.float();
    
    var imgBandNorm = imgBand.divide(Max); // normalize the data to 0 - 1
    
    return imgBandNorm;
    
  }

  /////////////////////////////////////////////////////////////
  // This function does cross-normalization for S1 images for use in multitemporal analysis
  // First, truncates histogram to the 98-percentile, then does a linear rescale to [0-1] (normalization)
  // histogram equalization between multitemporalimages is mapped through another function (crossNormalizeS1b)
  // (Dellepiane et al., 2012)
  function crossNormalizeS1a(image) {
    // var precVV = image.select('VV').reduceRegion({
    //                                               reducer: ee.Reducer.percentile([0,98]),
    //                                               geometry: aoi,
    //                                               scale: 10,
    //                                               tileScale: 16,
    //                                               maxPixels: 1e13
    //                                               });
    // var precVH = image.select('VH').reduceRegion({
    //                                               reducer: ee.Reducer.percentile([0,98]),
    //                                               geometry: aoi,
    //                                               scale: 10,
    //                                               tileScale: 16,
    //                                               maxPixels: 1e13
    //                                               });
  
    // var imageVV = normalize(image, 'VV', precVV.get('VV_p0'), precVV.get('VV_p98'));
    // var imageVH = normalize(image, 'VH', precVH.get('VH_p0'), precVH.get('VH_p98'));
    
    var imageVV = normalize(image, 'VV', 0.0000187, 0.1643852); // p0 for whole Greece for S2 and S3
    var imageVH = normalize(image, 'VH', 0.0000066, 0.0301957); // p98 for whole Greece for S2 and S3
    
    var rescaled = imageVV.addBands(imageVH);
    
    return rescaled.copyProperties(image).set(interval, image.get(interval));
  }
  
  // won't run for whole Greece for scale 10
  function crossNormalizeS1b(image) {
    var meansStdDevs = image.reduceRegion({
                                          reducer: ee.Reducer.mean().combine(ee.Reducer.stdDev(),null,true),
                                          geometry: aoi,
                                          scale: 500, // 10,
                                          tileScale: 16,
                                          maxPixels: 1e13
                                          });  
    var means = ee.Image.constant(meansStdDevs.values(["VV_mean", "VH_mean"]));
    var stdDev = ee.Image.constant(meansStdDevs.values(["VV_stdDev", "VH_stdDev"]));
    var normalized = image.select(['VV', 'VH']).subtract(means).divide(stdDev);
    return normalized.copyProperties(image).set(interval, image.get(interval));
  }


  //====================================================================================
  
  //--------------- MAIN ALGORITHM --------------------------------

  // ---------- CREATE IMAGE COLLECTION
  
  // filter S1 by date range, and area of interest
  var ColS11 = filterColS1(startDate1, endDate1, aoi);
  var ColS12 = filterColS1(startDate2, endDate2, aoi);
  var ColS1 = ColS11.merge(ColS12);
  print(ColS1)
  // Map.addLayer(ColS1.first(),{},'S1 collection')
  
  // add month and season property
  ColS1 = ColS1.map(addMonthSeason);
  // print('add month and season property to S1 collection:', ColS1)
  
  var ColMosS1 = mosaicBySeasonS1(ColS1); //create image collection with mosaics by season
  // print('number of seasonal mosaics for training (S1):', ColMosS1.size())
  // print('mosaic by season S1 IC:', ColMosS1)
  
  // mask outside AOI
  ColMosS1 = ColMosS1.map(maskOut);
  
  // cross-Normalize
  ColMosS1 = ColMosS1.map(crossNormalizeS1a);
  ColMosS1 = ColMosS1.map(crossNormalizeS1b);

  var bands_listS1 = ['VV',
                      'VH',
                      ];

  // select only the specific bands from the mosaic collection
  ColMosS1 = ColMosS1.select(bands_listS1); 
  // Map.addLayer(ColMosS1, {}, 'ColMosS1')
  
  // first rename all the bands by date
  function renameBands(image){
    var int = ee.String(image.get(interval));
    var new_bands = [ (int.cat('VV')),         // 1
                      (int.cat('VH')),         // 2
                      ];     
    return image.select(bands_listS1,new_bands);
  }
  ColMosS1 = ColMosS1.map(renameBands);
  // print('rename S1 bands:', ColMosS1)
  
  // set images ids to null, so ".toBands()" will not assign image id to band names
  ColMosS1 = ColMosS1.map(function x(im){
                        im = im.set({'system:index':''});
                        return im;
                      });
  
return ColMosS1;

}

exports.seasonalMosS1 = seasonalMosS1;
