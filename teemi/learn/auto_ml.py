#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

""" This module is focused on Machine Learning purposes. And simplifying these workflows"""

import h2o
from h2o.automl import H2OAutoML
import pandas as pd


def autoML_on_partitioned_data(
    feature_cols: list,
    training_column: str,
    df_input_for_ml,
    path="",
    partitions=5,
    training_time=5,
    nfold=10,
) -> None:
    """Runs over a pandas dataframe and trains MLs according to specified partition lenght.

    Parameters
    ----------
    feature_cols: list
        the column you want to train on fx : feature_cols = ['0', '1', '2', '3']
    training_column : str
        the column you want the models to be able to predict
    df_input_for_ml : pd.DataFrame
        A pandas dataframe with your data
    training_time : int
        the amount of time you want to train the models.
        If you set it to 0, it will not have a time limit and will train untill it reaches saturation i.e. best models

    Returns
    -------
    csv files in the specified path
    """

    all_mae = []

    # partitioning
    step = int(len(df_input_for_ml) / partitions) + 1
    partitions = [i for i in range(0, len(df_input_for_ml), step)]

    # Partion columns  - used for getting the right output
    partitions_col = partitions[1:]
    partitions_col.append(len(df_input_for_ml))

    # INCREASING THE SIZE OF THE DATASET
    partitions_list = [
        df_input_for_ml[partitions[0] : partitions[i]]
        for i in range(1, len(partitions))
    ]
    # add the last_full partition
    partitions_list.append(df_input_for_ml[partitions[0] :])

    ### Making the dataframes into h2o dfs
    list_of_df_test_frames = []

    for df in partitions_list:
        # initialize a h20 dataframe
        df_test = h2o.H2OFrame(pd.concat([df], axis="columns"))

        # changing columns to strings
        for col in df_input_for_ml.columns:
            if col != training_column:
                col = str(col)

        # making the dataframes categorical except the training column
        for column in df_test.columns:
            if col != training_column:
                df_test[column] = df_test[column].asfactor()
        list_of_df_test_frames.append(df_test)

    ##### setting up ML
    autoML_dataclasses_list = []

    # Initialize 5 - H2O autoML class
    for i in range(len(list_of_df_test_frames)):
        AutoML = H2OAutoML(
            max_runtime_secs=training_time,  # 1 hour =int(3600 * 1) , if unlimited time is wanted then set this to zero = 0
            max_models=None,  # None =  no limit
            nfolds=nfold,  # number of folds for k-fold cross-validation (nfolds=0 disables cross-validation)
            seed=1,  # Reproducibility
            sort_metric="MAE",
            keep_cross_validation_predictions=True,
        )
        autoML_dataclasses_list.append(AutoML)

    ##### Training the models on partitioned data
    for i in range(len(autoML_dataclasses_list)):
        autoML_dataclasses_list[i].train(
            x=feature_cols, y=training_column, training_frame=list_of_df_test_frames[i]
        )

        print(
            "len of dataframes that are being trained on :",
            len(list_of_df_test_frames[i]),
        )

    ### getting the mae for each model
    model_name = []
    cv_sd_mae = []
    cv_mean_mae = []
    best_models_mae = []
    for model in autoML_dataclasses_list:
        # Mae for each model train
        best_model = model.get_best_model()
        best_models_mae.append(best_model.mae())

        # CV metrics
        best_model_cv_summary = (
            best_model.cross_validation_metrics_summary().as_data_frame()
        )
        mean = float(best_model_cv_summary.iloc[0:1, 0:3]["mean"])
        sd = float(best_model_cv_summary.iloc[0:1, 0:3]["sd"])
        # save ot
        cv_mean_mae.append(mean)
        cv_sd_mae.append(sd)

        ## save names
        model_name.append(best_model.model_id)

    # saving ALL maes
    all_mae.append(best_models_mae)
    df = pd.DataFrame(all_mae, columns=partitions_col, dtype=float)
    df = df.T

    # add cv mean mae and sd
    df["CV_mean_MAE"] = cv_mean_mae
    df["CV_SD_MAE"] = cv_sd_mae
    df["Model_name"] = model_name

    # getting a unique name
    from datetime import datetime

    now = datetime.now()  # current date and time
    time = now.strftime("%Y_%m_%d_%H:%M")

    df.to_csv(path + time + "_ml_models_running_over_partioned_data.csv")
