Heterogeneity in the effect of short-term exposure to fine particulate matter on cardiovascular diseases
================
Emma Thomas
30/01/2020

This repository contains code for processing data and performing analysis for the paper "Heterogeneity in the effect of short-term exposure to fine particulate matter on cardiovascular diseases".

# Directories

  - `code`
      - The code directory should be used to store all code used as part
        of the project. If there are a large number of code files
        additional directories within the `code` directory are
        recommended, especially if there are multiple workflows. There
        should also be a description of the order code should be run in,
        or some other description of the contents and a method of
        indicating the workflow.
  - `data`
      - This directory should contain all raw data and processed data
        used by a given project. All csv and rds files within this
        directory are assumed to be large and by default are not tracked
        by git. If there are a large number of data files additional
        internal structure is recommended. The `.gitignore` file will
        need to be updated if more structure is used. One common
        paradign is to have a `data/raw_data` directory storing unedited
        files as receinved from the source and a `data/analysis_data`
        folder containing the assembled and cleaned data ready for use
        with models.
  - `figures`
      - This directory should contain all of the figures generated
        through the course of reaserch. The readme file in this
        directory should list the files, have a brief discription of the
        figures, and list the file that creates the figure. This
        directory should be tracked by git.
  - `reports`
      - This directory should contain Rmarkdown files used to summarize
        and describe reserach processes and data features. The HTML or
        PDF outputs of the markdown files should also be stored in this
        directory. This directory should be tracked by git.
  - `results`
      - This directory should be used to store objects such as models,
        tables, or other similar products of analysis. Whether or not
        files in this directory should be tracked by github is a
        question of their size and sensitivity and likely varies by
        project. By default, .RDS files in this directory will not be
        tracked, but all other files will.
