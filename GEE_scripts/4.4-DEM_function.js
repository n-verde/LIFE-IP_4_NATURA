/*=============================================================================

==== LIFE IP 4 NATURA ===== 
Natalia Verde, AUTH/DUTH, 2020

BRIEF DESCRIPTION:
This GEE function adds elevation slope and TRASP to a given image collection.

*/

//====================================================================================

//---------------FUNCTIONS------------------------------

  /////////////////////////////////////////////////////////////
  // This function  takes an image and then scales the output to uint16.
  // ex. rescale(img, 'NDVI', img.geometry())
  // ATTENTION! THIS IS RESCALE, NOT STRETCH!!!
  
  function rescale(img, bandName, aoi, min, max) {
    
    var Min = ee.Number(min);
    var Max = ee.Number(max);
    
    var imgBand = img.select(bandName);
      
    if (Min.lt(0)) {
      imgBand = imgBand.add(Min.abs()); // turn negative values to possitive
      Max = Max.add(Min.abs()); // change new max 
    }
      
    imgBand = imgBand.float();
    
    var imgBandNorm = imgBand.divide(Max); // normalize the data to 0 - 1
    var imgBandScale = imgBandNorm.multiply(65535); // Now scale to 65535
    var rescaledBand = imgBandScale.uint16();
    
    return rescaledBand;
    
  } 

/////////////////////////////////////////////////////////////
// This function calculates TRASP from an aspect image. The result is 
// a continuous variable between 0 - 1
// TOPOGRAPHIC SOLAR-RADIATION INDEX (TRASP)
// (Roberts and Cooper, 1989)

function TRASP(aspectImage) {
  var trasp =  null;
  var temp = aspectImage.expression('(pi/180)*(a-30)',{
            'a': aspectImage,
            'pi': 3.14159265359
            });
  trasp = (ee.Image(1).subtract(temp.cos())).divide(ee.Image(2));
  return trasp.rename('TRASP');
}

//====================================================================================

//--------------- MAIN FUNCTION --------------------------

function addDEM(inputIC,DEM,AOI) {

  // ---------- Add the DEM to the IC

  // match projection of IC to DEM projection
  var proj = DEM.projection().getInfo();
  var transform = proj.transform;
  var crs = proj.crs;

  // clip DEM to AOI
  var dem = DEM.clip(AOI);
  
  // elevation -->
  var demRepr = dem.reproject(crs, null, 25);
  var demres = dem;
  // demres = rescale(dem,'b1',AOI,-11.993,2830.168).toUint16(); // turn to unsigned 16-bit
  
  inputIC = inputIC.addBands(demres.rename('elev'));
  
  // slope-->
  var slope = ee.Terrain.slope(dem);
  var slopeRepr = slope.reproject(crs, null, 40); // reproject to 40m to remove artifacts
  var slopeRes = slopeRepr;
  // slopeRes = rescale(slopeRepr,'slope',AOI,0,90).toUint16(); // turn to unsigned 16-bit
  
  inputIC = inputIC.addBands(slopeRes);

  // aspect-->
  var aspect = ee.Terrain.aspect(dem);
  
  // trasp -->
  var TR = TRASP(aspect);
  
  var TRRepr = TR.reproject(crs, null, 40); // reproject to 40m to remove artifacts
  var TRRes = TRRepr;
  // TRRes = rescale(TRRepr,'TRASP',AOI,0,1).toUint16(); // turn to unsigned 16-bit
  TRRes = TRRepr.float(); // turn to float from double
  
  // final image collection with indices, tc, bci, dem etc...  
  var ICfinal = inputIC.addBands(TRRes);
  
  return ICfinal;

}

exports.addDEM = addDEM;
