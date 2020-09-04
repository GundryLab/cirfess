# cirfess
An online tool for aiding in proteomic strategies of the cell surface and extracellular proteome.  Data compilation by Matt Waas, R Shiny interface taken from SurfaceGenie (https://github.com/GundryLab/SurfaceGenie)

Everything needed to replicate the R Shiny app, CIRFESS, is found in this git repository.  These are the instructions to recreate it.  Other instructions on how to update it are in the readme in the data directory.   

Directory Structure:

```
root - R scripts to run the site
  |
  | 
  + - data - python scripts to create the sqlite3 database used by CIRFESS
  |    |       
  |    + - sources - data from phobius, tmhmm, etc to populate the db
  |
  + - www - images, stylesheets, etc for the display of the website
       | 
       + - data - json files for the protvista style display of peptides
```
Instructions -

- clone the git repo from github
- run data/makeLiteSingle.py to create the sqlite3 db, cirfess.db 
  - This will take a while
- run data/createJSON.py to create the json files in www/data for the protvista display.
  - This will take a long while
  - If it errors out, just run again.  It will pick up from where it left off.

