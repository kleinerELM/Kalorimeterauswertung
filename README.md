# Jupyter notebook to automatically preprocess the measurements of a TAM Air Calorimeter

This script loads exported xls files from the TAM Air software, crops the dataset, creates PNG/SVG graphs and creates a summary XLSX file containing all measurements in a folder.

The script searches for the first value for "Normalized heat", where the values are constantly larger than 0 J/g and drops all values before.

## requirements
- Python >= 3.11
- Jupyter
- pandas
- matplotlib
- tk
- xlswriter
- xlrd
