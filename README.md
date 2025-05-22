# Lepton AI PID
The first step is to generate and reconstruct the events used for training. Flat distributions of the signal and the background.
The next steps are assuming each sample for signal and for background are stored in root files that contains all the relevant information for training.
## Training
Make sure that the location of the files in lines 35-36 is correct. `TMVAClassification.C` will train using the MLP and BDT methods using all 9 variables. To edit this just remove or add the variables as needed. You can then run it using:
```
root -l TMVAClassification.C
```
## Testing
Use `TMVAClassificationApplication` to test the model. The result is a root file that contains the variables and the score assigned from each model.
```
root -l TMVAClassificationApplication.C
```

## Visualization of data
Modify the main funtion inside `plot.C`.

-`Variable_Plots_TMVA`: It plots the distributions of the 9 variables overlapping: True Positives and True Negatives, True Positives and False Positives, Negative Sample and False Positives, and Positive Sample and True Positives.

-`True_False`: Plots only True Positive vs False Positive, True Negative vs False Negative.

-`ROC`: Plots the ROC of the dataset. This is for simulations only.

-`Print_Table`: Displays True Positives, False Positives, True Negatives and False Negative values on screen. 
