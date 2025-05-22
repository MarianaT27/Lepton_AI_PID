# Lepton AI PID

The first step is to generate and reconstruct the events used for training. These should follow flat distributions for both signal and background.

The next steps assume that each signal and background sample is stored in ROOT files containing all the relevant information for training.

## Training

Make sure the file paths in lines 35–36 of `TMVAClassification.C` are correct. This script trains models using the MLP and BDT methods with all 9 variables. To customize, simply remove or add variables as needed. Run the training with:

```
root -l TMVAClassification.C
```

## Testing

Use `TMVAClassificationApplication.C` to test the trained models. The output is a ROOT file that contains the input variables and the scores assigned by each model:

```
root -l TMVAClassificationApplication.C
```

## Validation

Validation is performed both on the signal (signal preservation) and the background (background suppression).

`Rad_electron.C` is an older version of the code and can be used as a reference for basic testing. For the most up-to-date validation, use `Rad_analysis.C` and `Background_updated.C`. These scripts use `clas12root`, but they may also run under standard `root` if the CLAS12-specific libraries are removed.

### Signal Validation

This method relies on the fact that leptons radiate photons as they travel through the detector, allowing association between a lepton and its corresponding photon. A peak at zero in the Δθ (DeltaTheta) distribution indicates a match.

To validate the signal, use `Rad_analysis.C`:

- **Parameter 1**: Name of the file  
- **Parameter 2**: Configuration (`Fall18inbending = -18`, `Fall18outbending = 18`, `Spring19 = -19`)  
- **Parameter 3**: BDT model used (6 or 9 variables)  
- **Parameter 4**: Lepton type (`electron = 11`, `positron = -11`)  

```
clas12root 'Rad_analysis.C("pos_S19_BDT", -19, 6, -11)' -q -b
```

### Background Validation

For background validation, the process uses the reaction `ep → e⁻ π⁺ (X)`, where the pion is misidentified as an electron (`PID = -11`). The parameters are:

- **Parameter 1**: Name of the file  
- **Parameter 2**: Configuration  
- **Parameter 3**: BDT model used  
- **Parameter 4**: Training tag  
- **Parameter 5**: Lepton type  

```
clas12root 'Background_updated.C("pos_S19_BDT", -19, 6, "jpsitcs", -11)' -q -b
```

## Visualization of Data

Modify the main function inside `plot.C` to select which plots to generate:

- **`Variable_Plots_TMVA`**: Plots the distributions of the 9 variables, showing combinations such as:
  - True Positives vs. True Negatives  
  - True Positives vs. False Positives  
  - Negative Sample vs. False Positives  
  - Positive Sample vs. True Positives

- **`True_False`**: Plots:
  - True Positive vs. False Positive  
  - True Negative vs. False Negative

- **`ROC`**: Plots the ROC curve of the dataset (simulations only).

- **`Print_Table`**: Displays counts of True Positives, False Positives, True Negatives, and False Negatives in the terminal.
