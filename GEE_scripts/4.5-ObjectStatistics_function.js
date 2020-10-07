/*=============================================================================

==== LIFE IP 4 NATURA ===== 
Natalia Verde, AUTH/DUTH, 2020

BRIEF DESCRIPTION:
This GEE function computes statistcs of an image fora each object.
Inputs:
1. the image segmented image (created by GEE SNIC)
2. the image to calculate stdDev and Mean based on the clusters

Returns:
1. the seeds for the segmentation
2. the perimeters of the objects
3. the clusters
4. an object image with mean and stdDev properties for each band (so returns an image with twice the bands)

*/

//====================================================================================

//--------------- MAIN FUNCTION --------------------------

function objStats(objectImage, aoi, ICfinal) {
  
  // create superpixels grid: how many pixels away will the next superpixel be?
  // must be the same for a specific AOI
  var seedgrid = 10; // superpixel every 10 pixels (100m in 10m S2 image)
  var seeds = ee.Algorithms.Image.Segmentation.seedGrid(seedgrid); 
  
  // object image resulting from SNIC
  var clusters = objectImage.rename('clusters').clipToCollection(aoi);
  
  // perimeter of clusters
  var minMax = clusters.reduceNeighborhood(ee.Reducer.minMax(), ee.Kernel.square(1));
  var perimeterPixels = minMax.select(0).neq(minMax.select(1)).rename('perimeter');
  
  // Compute the perimeter of each cluster.
  var perimeter = perimeterPixels
    .addBands(clusters).reduceConnectedComponents({
      reducer: ee.Reducer.sum(),
      labelBand: 'clusters',
      maxSize: 100
    }).rename('perimeter');
  
  // Compute the area of each cluster.
  var area = ee.Image.pixelArea()
    .addBands(clusters).reduceConnectedComponents({
      reducer: ee.Reducer.sum(),
      labelBand: 'clusters',
      maxSize: 100
    }).rename('area');
  
  var tempI = area.addBands(perimeter);
  
  // Compute the form factor (Jiao & Liu, 2012)
  var form = tempI.expression('(4*pi*A)/(P**2)', {
                    'A': tempI.select('area'),
                    'P': tempI.select('perimeter'),
                    'pi': ee.Number(3.14159265359)
                    }).rename('form').multiply(100);
  
  // Compute the square pixel metric (Jiao & Liu, 2012)
  var square = tempI.expression('1-(4*(A**0.5))/P', {
                    'A': tempI.select('area'),
                    'P': tempI.select('perimeter')
                    }).rename('square').multiply(1000);
  
  // Compute the fractal dimension (Jiao & Liu, 2012)
  var fractal = tempI.expression('2*log(P/4)/log(A)', {
                    'A': tempI.select('area'),
                    'P': tempI.select('perimeter')
                    }).rename('fractal').multiply(1000);
  
  // Compute the shape index (Jiao & Liu, 2012)
  var shape = tempI.expression('(P/(4*(A**0.5)))', {
                    'A': tempI.select('area'),
                    'P': tempI.select('perimeter')
                    }).rename('shape').multiply(1000);
  // Compute per-cluster stdDev
  var stdDev = (ICfinal.addBands(clusters)).reduceConnectedComponents(ee.Reducer.stdDev(), 'clusters', 100);
  // Map.addLayer(stdDev, {}, 'per cluster StdDev', false)
  
  // Compute per-cluster Mean
  var means = (ICfinal.addBands(clusters)).reduceConnectedComponents(ee.Reducer.mean(), 'clusters', 100);
  // Map.addLayer(stdDev, {}, 'per cluster mean', false)
  
  var objectPropertiesImage = ee.Image.cat([
    stdDev,
    means,
    perimeter.multiply(10),
    area.divide(10),
    form,
    square,
    fractal,
    shape
  ]).float();
  // bands with _1 in their name are the MEANs, bands with origina names are STDDEVs
  
  // image for sampling and classification becomes objectPropertiesImage
  return {seeds:seeds, clusters:clusters, perimeterPixels:perimeterPixels, objectPropertiesImage:objectPropertiesImage};
  
}
exports.objStats = objStats;