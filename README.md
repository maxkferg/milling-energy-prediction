# Energy Consumption Prediction with GPR
An energy consumption model for a milling machine, based on Gaussian Process Regression (GPR). We use the MATLAB-PMML package to save the eergy consumption models in the PMML file format.

## About
The main purpose of the codein this repository is to:

* Train a Gaussian Process Regression Model to predict the energy consumption of a milling machine
* Save the trained models to PMML using the [MATLAB PMML Package](https://github.com/maxkferg/matlab-pmml)
* Make energy predictions using the PMML model files

## PMML
The energy consumption models are saved to the models directory. Each of the models uses a different number of features to predict the energy consumption of the milling machine. For this reason, some of the models will perform better than others.

## Usage
Clone the repository and setup the dependancies
```
git clone https://github.com/maxkferg/milling-energy-prediction
git submodule init
git submodule update
```

From the MATLAB editor:

1. Open the root folder with MATLAB
2. Run the `train_energy_model.m` file from MATLAB to generate PMML models
3. Run the `predict_energy_usage.m` file from MATLAB to make predictions

## Credits
The energy prediction model was first written by [@JinkyooPark](https://www.researchgate.net/profile/Jinkyoo_Park) and used as the basis for several publications. The PMML support was added by [@maxkferg](https://www.researchgate.net/profile/Max_Ferguson) for "Gaussian Process Regression (GPR) Representation in Predictive Model Markup Language (PMML)" [[link](https://compass.astm.org/DIGITAL_LIBRARY/JOURNALS/SSMS/PAGES/SSMS20160008.htm)].

## License
MIT
