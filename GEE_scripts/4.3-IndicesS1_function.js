/*=============================================================================

==== LIFE IP 4 NATURA ===== 
Natalia Verde, AUTH/DUTH, 2020

BRIEF DESCRIPTION:
This GEE function creates an image collection of Sentinel-1 IW GRD images and indices 
for an area of interest (AOI) to later use in training and classification.

Indices calculated and added to the collection are:
-VV/VH ratio
-Correlation GLCM texture metric for VV, 7x7
-Correlation GLCM texture metric for VH, 7x7

*/


//====================================================================================

//--------------- MAIN FUNCTION --------------------------

function IndicesS1(IC,aoi) {

  //--------------- FUNCTIONS ------------------------------
  
  /////////////////////////////////////////////////////////////
  // This function calculates VV-VH ratio, and adds it to the image stack.
  function addRATIO(image) {
    var ratio = image.select('VV').subtract(image.select('VH'))
                      .rename('RATIO'); //change new band name
  
    return image.addBands(ratio.rename('RATIO')); //change new band name
  
  }
  
  /////////////////////////////////////////////////////////////
  // This function computes correlation and homogeneity texture metric from the Gray Level Co-occurrence Matrix.
  // neighborhood size: 7x7, band: VV
  function addVVCorr7(image) {
    
    // var band = image.select('VV');
    var band = image.select('VV').expression("10000*b(0)").int32(); // turn to int32
    var glcm = band.glcmTexture({size:7});
    var correlation = glcm.select(3).rename('CorVV');

    // turn to float and rescale to same range of other indices
    correlation = correlation.expression("b(0)/1e9").float();
    
    return image.addBands(correlation.rename('CorVV'));
  }
  
    /////////////////////////////////////////////////////////////
  // This function computes correlation and homogeneity texture metric from the Gray Level Co-occurrence Matrix.
  // neighborhood size: 7x7, band: VH
  function addVHCorr7(image) {
    
    // var band = image.select('VH');
    var band = image.select('VH').expression("10000*b(0)").int32(); // turn to int32
    var glcm = band.glcmTexture({size:7});
    var correlation = glcm.select(3).rename('CorVH');
    
    // turn to float and rescale to same range of other indices
    correlation = correlation.expression("b(0)/1e9").float();
    
    return image.addBands(correlation.rename('CorVH'));
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

  var ColMosS1 = ee.ImageCollection(IC);
  
  // mask outside AOI
  ColMosS1 = ColMosS1.map(maskOut);
  
  // ---------- INDICES S1 ---------
  
  // calculate VV/VH (ratio)
  ColMosS1 = ColMosS1.map(addRATIO);
  // print('add VV/VH RATIO to collection:', ColMosS1)
  
  // add Correletion texture metrics for VV and VH
  ColMosS1 = ColMosS1.map(addVVCorr7);
  ColMosS1 = ColMosS1.map(addVHCorr7);
  
  var bands_listS1 = ['VV',
                      'VH',
                      'RATIO',
                      'CorVV',
                      'CorVH'
                      ];  
                      // bands from which to train model

  // select only the specific bands from the mosaic collection
  var ColIndS1 = ColMosS1.select(bands_listS1); 
  // print('select S1 indices:', ColIndS1)
  // Map.addLayer(ColIndS1, {}, 'ColIndS1')
  
  // first rename all the bands by date
  function renameBands(image){
    var season = ee.String(image.get('season'));
    var new_bands = [ (season.cat('VV')),         // 1
                      (season.cat('VH')),         // 2
                      (season.cat('RATIO')),      // 3
                      (season.cat('CorVV')),      // 4
                      (season.cat('CorVH'))       // 5
                      ];     
    return image.select(bands_listS1,new_bands);
  }
  ColIndS1 = ColIndS1.map(renameBands);
  // print('rename S1 bands:', ColIndS1)
  
  // set images ids to null, so ".toBands()" will not assign image id to band names
  ColIndS1 = ColIndS1.map(function x(im){
                        im = im.set({'system:index':''});
                        return im;
                      });
  
return ColIndS1;

}

exports.IndicesS1 = IndicesS1;
