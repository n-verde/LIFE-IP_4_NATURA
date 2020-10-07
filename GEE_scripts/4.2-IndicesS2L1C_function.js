/*=============================================================================

==== LIFE IP 4 NATURA =====
Natalia Verde, AUTH/DUTH, 2020

BRIEF DESCRIPTION:
This GEE function creates an image collection of Tasseled Cap Transformation Brightness, Greenness and Wetness bands of Sentinel-2 L1C data, 
as well as the BCI index.

Indices calculated and added to the collection are:
-TC Brightness
-TC Greenness
-TC Wetness
-BCI
-Correlation GLCM texture metric for BCI, 3x3

*/

//====================================================================================
  
// ---------- MAIN FUNCTION ----------

function IndicesS2L1C(IC, aoi){

  //====================================================================================
  
  //---------------FUNCTIONS------------------------------
  
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
  // This function calculates Tasseled Cap Components from a S2 L1C image, (Brightness, Greenness, and Wetness)
  // IMAGES MUST BE TOA AND NOT REFLECTANCE!
  // and adds it to the image stack.
  // (Shi, 2019)
  
  function calculateTasseledCap(image){
    var b = null; var brightness_ = null; var greenness_= null; var wetness_= null; var sum = null;
    var brightness = null; var greenness = null; var wetness = null; var tct = null;
    
    b = image.select([0,1,2,3,4,5,6,7,8,9,10,11,12]); // ATTENTION!!! MUST MATCH BANDS IN PAPER (Shi, 2019)!
    
    b = b.multiply(Math.pow(2,16)); // turn to uint16 from float
    
    brightness_= ee.Image([ 0.2381, 0.2569, 0.2934, 0.3020, 0.3099, 0.3740, 0.4180, 0.3580, 0.3834, 0.0103, 0.0020, 0.0896, 0.0780]);
     greenness_= ee.Image([-0.2266,-0.2818,-0.3020,-0.4283,-0.2959, 0.1602, 0.3127, 0.3138, 0.4261, 0.1454,-0.0017,-0.1341,-0.2538]);
       wetness_= ee.Image([ 0.1825, 0.1763, 0.1615, 0.0486, 0.0170, 0.0223, 0.0219,-0.0755,-0.0910,-0.1369, 0.0003,-0.7701,-0.5293]);
    sum = ee.call('Reducer.sum');
    
    brightness = b.multiply(brightness_).reduce('sum');
    greenness = b.multiply(greenness_).reduce('sum');
    wetness = b.multiply(wetness_).reduce('sum');
    
    tct = ee.Image(brightness).addBands(greenness).addBands(wetness);
    tct = tct.select([0,1,2],['Bri','Gree','Wet']);
     
    return image.addBands(tct);
  }
  
  /////////////////////////////////////////////////////////////
  // This function is used for rescaling TCT double bands to float.
  // Used for mapping over collection
  function mappingFunRescaleTC(image) {
    
    var Bres = image.select('Bri').divide(Math.pow(2,16)).float();
    var Gres = image.select('Gree').divide(Math.pow(2,16)).float();
    var Wres = image.select('Wet').divide(Math.pow(2,16)).float();
    
    var tctres = image.select('BCI');
    tctres = tctres.addBands(Bres.rename('Bri')).addBands(Gres.rename('Gree')).addBands(Wres.rename('Wet'));
      
    return tctres;
    
  }
  
  /////////////////////////////////////////////////////////////
  // This function calculates BCI from a S2 or L8 image, and adds it to the image stack.
  // BCI, (Biophysical Composition Index), based on the Tasseled Cap transformation.
  // BCI is able to discriminate correctly between urban and bare land covers. In addition, all other 
  // land cover types appear with zero or negative values.
  // (Deng & Wu, 2012)
  
  // input image must already contain bands named 'Brightness','Greenness','Wetness' !!!!
  // if needed, define the scale and the geometry of the reducers (lines 19-24), for faster processing
  
  function addBCI(image){
    
    var bci = null;
    
    // normalize to [0,1] as indicated in (Deng & Wu, 2012)
    var tc1 = image.select('Bri');
    var TC1 = normalize(tc1,'Bri',1004.80,54465.84);
    var tc2 = image.select('Gree');
    var TC2 = normalize(tc2,'Gree',-10420.02,6874.25);
    var tc3 = image.select('Wet');
    var TC3 = normalize(tc3,'Wet',-15365.60,1955.01);
    
    var TC1min = ee.Image(ee.Number(4.7)); var TC1max = ee.Image(ee.Number(38));
    var TC2min = ee.Image(ee.Number(-11.5)); var TC2max = ee.Image(ee.Number(21));
    var TC3min = ee.Image(ee.Number(-22)); var TC3max = ee.Image(ee.Number(4));
    
    var H = (TC1.subtract(TC1min)).divide(TC1max.subtract(TC1min));
    var V = (TC2.subtract(TC2min)).divide(TC2max.subtract(TC2min));
    var L = (TC3.subtract(TC3min)).divide(TC3max.subtract(TC3min));
    
    image = image.addBands(H.rename('H')); 
    image = image.addBands(V.rename('V')); 
    image = image.addBands(L.rename('L'));
    
    bci = image.expression('((H+L)/2-V)/((H+L)/2+V)',{
              'H': image.select('H'),
              'V': image.select('V'),
              'L': image.select('L')
    }).rename('BCI');
    
    image = image.select(image.bandNames().removeAll(['H', 'V', 'L']));
    var bcires = bci.float();
  
    return image.addBands(bcires.rename('BCI'));
  }
  
  /////////////////////////////////////////////////////////////
  // This function computes correlation texture metric from the Gray Level Co-occurrence Matrix.
  // neighborhood size: 3x3, band: BCI
  function addBCICorr3(image) {
    
    var band = image.select('BCI').expression("10000*b(0)").int32(); // turn to int32
    var glcm = band.glcmTexture({size:3});
    var correlation = glcm.select(3).rename('CorBCI');

    // turn to float and rescale to same range of other indices
    correlation = correlation.expression("b(0)/1e6").float();
    
    // return image.addBands(corrres.rename('CorBCI'));
    return image.addBands(correlation.rename('CorBCI')).float();
  }
  
  //--------------- MAIN FUNCTION --------------------------

  var ColMos = ee.ImageCollection(IC);
  
  // mask outside AOI
  ColMos = ColMos.map(maskOut);
  
  // ---------- TCT ---------
  ColMos = ColMos.map(calculateTasseledCap);
  // print('add TCT bands to collection:', ColMos)
  
  // ---------- BCI ----------
  ColMos = ColMos.map(addBCI);
  // print('add BCI bands to collection:', ColMos)
  
  // // ---------- rescale TC to float
  ColMos = ColMos.map(mappingFunRescaleTC);
  
  // ----------- TEXTURE ---------
  ColMos = ColMos.map(addBCICorr3);
  
  // bands to rename
  var bands_array = ['Bri',
                     'Gree',
                     'Wet',
                     'BCI',
                     'CorBCI'
                     ];
  
  // select only the specific bands from the mosaic collection
  var ColTCBCI = ColMos.select(bands_array);
  
  // first rename all the bands by date
  function renameBands(image){
    var season = ee.String(image.get('season'));
    var new_bands = [ (season.cat('Bri')),        // 0
                      (season.cat('Gree')),       // 1
                      (season.cat('Wet')),        // 2
                      (season.cat('BCI')),        // 3
                      (season.cat('CorBCI'))      // 4
                      ];    
    return image.select(bands_array,new_bands);
  }
  ColTCBCI = ColTCBCI.map(renameBands);
  
  // returns a collection
  return ColTCBCI;
  
}

exports.IndicesS2L1C = IndicesS2L1C;