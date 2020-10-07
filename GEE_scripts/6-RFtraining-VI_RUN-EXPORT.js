/*=============================================================================

==== LIFE IP 4 NATURA - STEP 6/7 =====
Natalia Verde, AUTH/DUTH, 2020

BRIEF DESCRIPTION:
This GEE script trains a Random Forest classifier based on the training and validation samples imported by the used
and finds the most important variables.
mtry and ntree for each model have derrived from R 'randomForest' functions.

HOW TO USE: 
1. set the paths for training and validation samples from assets.
2. set other properties you want in the "PARAMETER SETTINGS" section.
3. in the tasks tab, click RUN to all tasks in order to export the error matrix and the variable importance to Google Drive.


*/

//====================================================================================

//===============PARAMETER SETTINGS===================================================

// ---------- SAMPLES ----------
var auto = true; // if samples were obtained in automatic manner (auto=true) or namual (auto=false)

// // Seasonal all (mtry = 94)
// var tr = ee.FeatureCollection('users/n-verde/auth/LIFE_IP_4_NATURA/GR_OB_TrainingSamples_aut_S_All_500');
// var va = ee.FeatureCollection('users/n-verde/auth/LIFE_IP_4_NATURA/GR_OB_ValidationSamples_man_S_All');

// Seasonal indices (only indices) (mtry = 42) 
var tr = ee.FeatureCollection('users/n-verde/shared/LIFE-IP_4_NATURA/LIFE_OB_TrainingSamples');
var va = ee.FeatureCollection('users/n-verde/shared/LIFE-IP_4_NATURA/LIFE_OB_ValidationSamples');

// // Monthly indices (monthly indices) (mtry = 94) 
// var tr = ee.FeatureCollection('users/n-verde/auth/LIFE_IP_4_NATURA/GR_OB_TrainingSamples_aut_S_Indices_350');
// var va = ee.FeatureCollection('users/n-verde/auth/LIFE_IP_4_NATURA/GR_OB_ValidationSamples_man_S_Indices');

// ---------- RANDOM FOREST PARAMETERS ----------
var ntree = 400; // value derrived from analysis in R 'randomForest'
// var mtry = 94; // value derrived from analysis in R 'tuneRF','randomForest'
var mtry = 42; // value derrived from analysis in R 'tuneRF','randomForest'

var randomSeed = 12344; // random seed for RF
var IVnum = 20; // number of most important variables to show in console

// ---------- ERROR MATRIX AND VARIABLE IMPORTANCE EXPORT DETAILS ----------
var nameSuffix = 'LIFE';

// ---------- OTHER SETTINGS ----------
var prop = 'class_L2';  // the property in the feature collection from which to train model

//====================================================================================

//---------------MAIN PROGRAM--------------------------------

print('start');

// ---------- SAMPLES ----------

var allTraining = tr;
var allValidation = va;

var bands = allTraining.first().propertyNames();
bands = bands.filter(ee.Filter.neq('item','system:index'));  // remove 'system:index' property
bands = bands.filter(ee.Filter.neq('item', prop));  // remove 'class_L2' property
var bandNumber = bands.length().getInfo();
print('Bands:', bands)

// ---------- remove NAs from samples ----------
allTraining = allTraining.filter(ee.Filter.notNull(bands));
allValidation = allValidation.filter(ee.Filter.notNull(bands));

// ---------- CLASSIFIER ----------

// // if you don't know the optimal mtry
// // mtry is usually the square root of the number of bands you use for classification
// var mtry =  Math.round(Math.sqrt(bandNumber));
// print('mtry value:',mtry)

// create the RF classifier
var rf = ee.Classifier.smileRandomForest({
      numberOfTrees: ntree,  // (ntree) most used value 500 
      variablesPerSplit: mtry,  // (mtry) is usually the square root of the number of bands you use for classification
      bagFraction: 0.632,  // 1 = take no OOB samples, 0.632 = default in R
      seed: randomSeed  // random seed
    });

// train the created classifier with the training image and training samples
var classifier = rf.train({
  features: allTraining,
  classProperty: prop,
  inputProperties: bands
  });
var classifierDetails = classifier.explain();
print('RF model and importance:', classifierDetails)
print('OOB error estimate:', classifierDetails.get('outOfBagErrorEstimate'))

// --------------------------- VARIABLE IMPORTANCE -------------------------------

var varImp = ee.Feature(null, ee.Dictionary(classifierDetails).get('importance'));
var varImpDict = varImp.toDictionary();
// print('Variable importance:', varImpDict)

// find most important
var variables = varImpDict;

var features = ee.FeatureCollection(variables.map(function(k,v) { return ee.Feature(null, {'band': k, 'importancy': v})}).values());
var mostImportant = features.sort('importancy', false).limit(IVnum);
// print('Most important variables for RF model:', mostImportant)

// mostImportant = mostImportant.distinct('band') // The keys in a dict cannot contain duplicates
var dic = ee.Dictionary.fromLists({
  keys:mostImportant.aggregate_array('band'), values: mostImportant.aggregate_array('importancy')});
print('Most important variables for RF model:', dic)

// ---------- CHART of variable importance
var chart =
ui.Chart.feature.byProperty(varImp)
.setChartType('ColumnChart')
.setOptions({
title: 'GEE Variable Importance',
legend: {position: 'none'},
hAxis: {title: 'Bands'},
vAxis: {title: 'Importance'}
});
print(chart)


// ==========================================================================

// --------------------------- ACCURACY ------------------------------------- 

var validated = allValidation.classify(classifier);

var list = ee.List.sequence(1,21);

var errorMatrix = validated.errorMatrix(prop, 'classification',list);
print('Error Matrix:', errorMatrix)

var OA = errorMatrix.accuracy();
print('Model OA:', OA)
var PA = errorMatrix.producersAccuracy();
print('Model PA (Omission accuracy):', PA)
var UA = errorMatrix.consumersAccuracy();
print('Model UA (Commission Accuracy):', UA)
var Kappa = errorMatrix.kappa();
print('Model Kappa:', Kappa)

// --------------------------- EXPORTS ------------------------

var exportErrorMatrix = ee.Feature(null, {matrix: errorMatrix.array()});

// Export Confusion Matrix
var descr = nameSuffix + '_errorMatrix';
Export.table.toDrive({
  collection: ee.FeatureCollection(exportErrorMatrix),
  description: descr,
  folder : 'GEE_output',
  fileFormat: 'CSV'
});

// Export Variable Importance
var descr = nameSuffix + '_variableImportance';
Export.table.toDrive({
  collection: ee.FeatureCollection(varImp),
  description: descr,
  folder : 'GEE_output',
  fileFormat: 'CSV'
});

print('end');

exports.classifier = classifier;