# constitexprImportdata

constitexprImportdata is a repository to keep my data import routines.

## List of files

- loadFunctionsAndPackages 

+ You can start with compareData.Rmd
- It is calling other scripts which are in the same git repo you will see:
  - loadFunctionsAndPackages
  - readMMData.R
  - transformMMData.R
-These contain the necessary codes, to import de MM data.

+The output of that are two dataframes:
-myframes_to_mycells : one row per observation
-mycells : one row per cell

+And you should enter the location of your data into the dataList csv file.
