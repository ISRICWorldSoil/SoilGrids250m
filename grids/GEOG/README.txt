Preliminary predictions SoilGrids at 250m (the same were aggregated to 1km). Produced using machine learning algorithms and fully documented via https://github.com/ISRICWorldSoil/SoilGrids250m/ 

Model fitting and predictions by: T. Hengl (tom.hengl@isric.org)
Server installation (Geonode) and maintaince: J. Mendes de Jesus (jorge.mendesdejesus@wur.nl)

Inputs: various national soil profile DBs which are either publicly available or have been given to ISRIC (for a complete list see: http://www.isric.org/data/wosis). As covariates we use global SoilGrids250m covariates (https://github.com/ISRICWorldSoil/SoilGrids250m/blob/master/grids/covs1t/SoilGrids250m_COVS250m.csv) ca 160 layers;
Outputs: Predicted values of soil properties at 7 standard depths (0, 5, 15, 30, 60, 100, 200 cm) based on 'ranger' and 'xgboost', predicted probabilities per class based on randomForest predictions (fitted using the 'ranger' and 'nnet' softwares)
Lineage: Preparation of point data is documented at: https://github.com/ISRICWorldSoil/SoilGrids250m/ and model fitting at: https://github.com/ISRICWorldSoil/SoilGrids250m/tree/master/grids/model. SoilGrids documentation can be also found at: http://www.isric.org/content/soilgrids

Disclaimer: These are preliminary (unvalidated) results subject to constant updates. Please, do NOT distribute these maps without consulting the authors (use only for testing purposes). See also: http://www.isric.org/content/disclaimer-soilgrids. Final maps will (most likely) be made publicly available under the Open Data Base License (http://opendatacommons.org/licenses/odbl/summary/) in mid 2016. To report a bug or artifact in the map please use: https://github.com/ISRICWorldSoil/SoilGrids250m/issues.