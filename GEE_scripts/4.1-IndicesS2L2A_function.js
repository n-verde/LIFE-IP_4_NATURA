/*=============================================================================

==== LIFE IP 4 NATURA =====
Natalia Verde, AUTH/DUTH, 2020

BRIEF DESCRIPTION:
This GEE function creates an image collection of Sentinel-2 L2A images and indices 
for an area of interest (AOI) to later use in training and classification.

Indices calculated and added to the collection are:
-reNDVI
-NDWI
-GRVI
-NDRBI
-MSI
-EVI
-Correlation GLCM texture metric for B2, 7x7
-PANTEX
 
*/

//====================================================================================

//--------------- MAIN FUNCTION --------------------------

function IndicesS2L2A(IC,aoi) {

  //--------------- FUNCTIONS ------------------------------
  
  ////////////////////////////////////////////////////////////
  // This function creates a mask outside the aoi
  function maskOut(image) {
      var mask = ee.Image.constant(0).float();
      mask = mask.paint(aoi,1);
      var newimage = image.updateMask(mask);
    return newimage;
  }
  
  /////////////////////////////////////////////////////////////
  // This function calculates red edge NDVI from a S2 image, and adds it to the image stack.
  // bands 8 (NIR) and 5 (RE1) (Forkuor, 2018)
  function addreNDVI(image) {
    var rendvi = (image.normalizedDifference(['B8', 'B5'])).rename('reNDVI');
    return image.addBands(rendvi.rename('reNDVI')).float(); //change new band name
  }
  
  /////////////////////////////////////////////////////////////
  // This function calculates NDWI from a S2 image, and adds it to the image stack.
  // NDWI (Normalized Difference Water Index) index is most appropriate for water body mapping.
  function addNDWI(image) {
    var ndwi = (image.normalizedDifference(['B3', 'B8'])).rename('NDWI');
    return image.addBands(ndwi.rename('NDWI')).float(); //change new band name
  }
  
  /////////////////////////////////////////////////////////////
  // This function calculates GRVI from a S2 image, and adds it to the image stack.
  // GRVI (Green-Red Vegetation Index) is for plant phenology
  function addGRVI(image) {
    var grvi = (image.normalizedDifference(['B3', 'B4'])).rename('GRVI');
    return image.addBands(grvi.rename('GRVI')).float(); //change new band name
  }
  
  /////////////////////////////////////////////////////////////
  // This function calculates NDRBI from a S2 image, and adds it to the image stack.
  // Normalized Difference Red-Blue Index (NDRBI), enhances red roofs
  // (Mallinis et. al, 2014)
  function addNDRBI(image) {
    var ndrbi = (image.normalizedDifference(['B4', 'B2'])).rename('NDRBI');
    return image.addBands(ndrbi.rename('NDRBI')).float(); //change new band name
  }
  
  /////////////////////////////////////////////////////////////
  // This function calculates MSI from a S2 image, and adds it to the image stack.
  // MSI, Moisture Stress Index (Hunt, 1989), (Dotzler, 2015)
  function addMSI(image) {
    var swir = image.select('B11');
    var nir = image.select('B8');
    
    var msi = (swir.divide(nir)).rename('MSI');
    return image.addBands(msi.rename('MSI')).float(); //change new band name
  }
  
  ///////////////////////////////////////////////////////////// 
  // This function calculates EVI from a S2 image, and adds it to the image stack.
  // EVI, (Enhanced Vegetation Index) 
  // (Huete et al., 2002) 
  function addEVI(image) {
    var evi = image.expression('2.5 * ((NIR - R) / (NIR + 6.0 * R - 7.5 * B + 1.0))',{
              'B': (image.select('B2')).divide(ee.Image(10000)),
              'R': (image.select('B4')).divide(ee.Image(10000)),
              'NIR': (image.select('B8')).divide(ee.Image(10000))
              }).rename('EVI');
    
    // turn to float and rescale to same range of other indices
    evi = evi.expression("1000*b(0)").float();
    return image.addBands(evi.rename('EVI')).float();
  }
  
  /////////////////////////////////////////////////////////////
  // This function computes correlation texture metric from the Gray Level Co-occurrence Matrix.
  // neighborhood size: 7x7, band: B2 (blue)
  function addB2Corr7(image) {
    
    // var band = image.select('B2');
    var band = image.select('B2').expression("10000*b(0)").int32(); // turn to int32
    var glcm = band.glcmTexture({size:7});
    var correlation = glcm.select(3).rename('CorB2');
    
    // turn to float and rescale to same range of other indices
    correlation = correlation.expression("b(0)/1e9").float();
    return image.addBands(correlation.rename('CorB2'));
  }
  
  /////////////////////////////////////////////////////////////
  // This function applies the PANTEX method, for built-up detection with GLCM.
  // PANTEX uses the blue band or the max of the RGB bands. Then calculates contrast in the 
  // GLCM and finds the intersection image of the different angles in which it is calculated.
  // FOR SENTINEL 2 IMAGES
  // FOR MAPPING OVER A COLLECTION!
  // (Pesaresi et al., 2008)
    
  function addPANTEX(image) {
    var angles = null; var pantex = null; var glcm = null;

    // option #1: blue band (Pesaresi et al., 2008)
    // var band = image.select('B2');
  
    // option #2: max of RBG bands (MASADA)
    var b = image.select('B2').expression("10000*b(0)").int32().rename('B2'); // turn to int32
    var g = image.select('B3').expression("10000*b(0)").int32().rename('B2'); // turn to int32
    var r = image.select('B4').expression("10000*b(0)").int32().rename('B2'); // turn to int32
    var bands = ee.ImageCollection.fromImages([b,g,r]);
    var band = bands.max().expression("10000*b(0)").int32().rename('B2'); // turn to int32
 
    var k = ee.Kernel.fixed(5, 5, [
                                [0, 0, 0, 1, 0],
                                [0, 0, 0, 1, 1],
                                [0, 0, 0, 1, 1],
                                [0, 0, 1, 1, 1],
                                [0, 0, 1, 1, 0]
                                ]);
    
    glcm = band.glcmTexture({size:5, kernel:k, average:false});
    var d1 = glcm.select('B2_contrast_1_-2').rename('angles');
    var d2 = glcm.select('B2_contrast_1_-1').rename('angles');
    var d3 = glcm.select('B2_contrast_2_-1').rename('angles');
    var d4 = glcm.select('B2_contrast_1_0').rename('angles');
    var d5 = glcm.select('B2_contrast_2_0').rename('angles');
    var d6 = glcm.select('B2_contrast_0_1').rename('angles');
    var d7 = glcm.select('B2_contrast_1_1').rename('angles');
    var d8 = glcm.select('B2_contrast_2_1').rename('angles');
    var d9 = glcm.select('B2_contrast_0_2').rename('angles');
    var d10 = glcm.select('B2_contrast_1_2').rename('angles');
    
    angles = ee.ImageCollection.fromImages([d1,d2,d3,d4,d5,d6,d7,d8,d9,d10]);
    // pantex = angles.min().rename('PANTEX');
    
    // turn to float and rescale to same range of other indices
    pantex = angles.min().expression("b(0)/1e17").float();
    return image.addBands(pantex.rename('PANTEX'));
  }

  //====================================================================================
  
  //--------------- MAIN ALGORITHM --------------------------------

  var ColMosS2 = ee.ImageCollection(IC);
  
  // ---------- INDICES S2 ---------
  
  // calculate reNDVI for each image and add it to stack    
  ColMosS2 = ColMosS2.map(addreNDVI);
  // print('add reNDVI band to collection:', ColMosS2)
  
  // calculate NDWI index and add as band to collection
  ColMosS2 = ColMosS2.map(addNDWI);
  // print('add NDWI band to collection:', ColMosS2)
  
  // calculate GRVI
  ColMosS2 = ColMosS2.map(addGRVI);
  // print('add GRVI band to collection:', ColMosS2)
  
  // calculate NDRBI 
  ColMosS2 = ColMosS2.map(addNDRBI);
  // print('add NDRBI band to collection:', ColMosS2)
  
  // calculate MSI
  ColMosS2 = ColMosS2.map(addMSI);
  // print('add MSI band to collection:', ColMosS2)
  
  // calculate EVI
  ColMosS2 = ColMosS2.map(addEVI);
  // print('add EVI band to collection:', ColMosS2)
  
  // calculate correlation texture
  ColMosS2 = ColMosS2.map(addB2Corr7);
  // print('add correlation band to collection:', ColMosS2)
  
  // calculate PANTEX                                                                 
  ColMosS2 = ColMosS2.map(addPANTEX);
  // print('add PANTEX band to collection:', ColMosS2)
  
  var bands_listS2 = [
                       'reNDVI',
                       'NDWI',
                       'GRVI',
                       'NDRBI',
                       'MSI',
                       'EVI',
                       'CorB2', 
                       'PANTEX',
                       'B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12'];  
                       // features from which to train model

  // select only the specific bands from the mosaic collection
  var ColIndS2 = ColMosS2.select(bands_listS2); 
  // print('select S2 indices:', ColIndS2)
  
  // first rename all the bands by date
  function renameBands(image){
    var season = ee.String(image.get('season'));
    var new_bands = [ 
                      (season.cat('reNDVI')), // 0 
                      (season.cat('NDWI')),   // 1 
                      (season.cat('GRVI')),   // 2 
                      (season.cat('NDRBI')),  // 3 
                      (season.cat('MSI')),    // 4
                      (season.cat('EVI')),    // 5
                      (season.cat('CorB2')),  // 6
                      (season.cat('PANTEX')), // 7
                      (season.cat('B2')),     // 8
                      (season.cat('B3')),     // 9
                      (season.cat('B4')),     // 10
                      (season.cat('B5')),     // 11
                      (season.cat('B6')),     // 12
                      (season.cat('B7')),     // 13
                      (season.cat('B8')),     // 14
                      (season.cat('B8A')),    // 15
                      (season.cat('B11')),    // 16
                      (season.cat('B12'))];   // 17
    return image.select(bands_listS2,new_bands);
  }
  ColIndS2 = ColIndS2.map(renameBands);
  
  // print('rename bands:', ColIndS2)
  
  // set images ids to null, so ".toBands()" will not assign image id to band names
  ColIndS2 = ColIndS2.map(function x(im){
                        im = im.set({'system:index':''});
                        return im;
                      });
  
return ColIndS2;

}

exports.IndicesS2L2A = IndicesS2L2A;
