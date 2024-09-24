import numpy as np
import xgboost as xgb
import shap
import scipy.io
import pandas as pd
import os
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler

    # Read the train label
    y_train = pd.read_excel('path_to_training_HC_data.xlsx')
    # Read the test label
    y_test = pd.read_excel('path_to_diagnostic_group_data.xlsx')

    # Load features for brain age prediction
    X_train = scipy.io.loadmat('training_features_path.mat')
    X_test = scipy.io.loadmat('testing_features_path.mat')

    # Normalize the features
    scaler = MinMaxScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.1, random_state=42)
    # Create training data object
    dtrain = xgb.DMatrix(data=X_train, label=y_train)
    dval = xgb.DMatrix(X_val, label=y_val)
    dtest = xgb.DMatrix(X_test)
    # Define hyperparameters for the model
    params = {
        'objective': 'reg:squarederror',  # loss function
        'eta': 0.05,  # learning rate
        'max_depth': 4,  # maximum depth of a tree
        'num_boost_round': 2000  # number of iterations
    }

    # training model
    model = xgb.train(params, dtrain, num_boost_round=params['num_boost_round'],
                      evals=[(dtrain, 'train'), (dval, 'validation')],
                      early_stopping_rounds=20)
    y_pred = model.predict(dtest)


    # Calculate correlation
    df = pd.DataFrame({'data1': y_test, 'data2': y_pred})
    correlation = df['data1'].corr(df['data2'])
    print("Correlation coefficient:", correlation)

    # Calculate MAE
    mae = np.mean(np.abs(y_test - y_pred))
    # Print MAE
    print("Mean Absolute Error (MAE):", mae)

    # Save model
    model.save_model('xgboost_model.model')
    # losd model
    loaded_model = xgb.Booster()
    loaded_model.load_model('xgboost_model.model')

    # Create SHAP interpreter
    explainer = shap.Explainer(loaded_model)

    # Explain model predictions
    shap_values = explainer.shap_values(X_test)

    # Print SHAP value
    print(shap_values)

    # Save SHAP value
    result_df = pd.DataFrame(shap_values)
    output_file_path = f':/output_path.xlsx'
    result_df.to_excel(output_file_path, index=False)


