# Notes of the evaluation workflow

## Predict with the model
Run from the root directory
`python c12a_predictor.py Eval_on_HT1-2/input_ht1_2.txt Eval_on_HT1-2/output_ht1_2.txt`

## Compare the prediction with the "truth" (experimentally measured indel frequencies)
Put the predicted results together with the "truth" from the publication Kim2018 in one file `Evalulation_Brien2023.xlsx`

Run the R script `Eval_on_HT1-2/draw_scatter_plot.R`
