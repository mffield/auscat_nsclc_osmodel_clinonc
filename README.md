# Federated learning survival model and potential radiotherapy decision support impact assessment for non-small cell lung cancer using real-world data. #

This is the code repository for the study:
"Federated learning survival model and potential radiotherapy decision support impact assessment for non-small cell lung cancer using real-world data" in Clinical Oncology 

The code developed in MATLAB utilised in this study analysis involving federated learning (FL) is contained in the folder "MATLAB_FL". It comprises a central algorithm script ("CentralLungDSS.m") which coordinates the FL and a local algorithm ("LocalLungModel.m") which runs on the infrastructure storing the data. This code is dependent on a functioning FL infrastructure utlilising the function in "Matlab_FL\@MPIservice" as described in Field et al. [2], [3], with code found here on Github: [AusCAT Federated Learning](https://github.com/mffield/auscat_federated_learning).

An addiitonal Python script is provided under the Python folder for validating the model. An example of a data set with dummy values is provided to demonstrate the script.

To cite this FL model development and decision support work please use this:

Field M, Vinod S, Aherne N, Carolan M, Dekker A, Delaney G, Greenham S, Hau E, Lehmann J, Ludbrook J, Miller A, Rezo A, Selvaraj J, Sykes J, Thwaites D, Holloway L, Federated learning survival model and potential radiotherapy decision support impact assessment for non-small cell lung cancer using real-world data, Clinical Oncology 2024.


## Contact ##

matthew.field@unsw.edu.au

matthew.field@aimedeng.com


## References ##

[1] Field M, Vinod S, Aherne N, Carolan M, Dekker A, Delaney G, Greenham S, Hau E, Lehmann J, Ludbrook J, Miller A, Rezo A, Selvaraj J, Sykes J, Thwaites D, Holloway L, Federated learning survival model and potential radiotherapy decision support impact assessment for non-small cell lung cancer using real-world data, Clinical Oncology 2024 (in press).

[2] Field M, Thwaites DI, Carolan M, Delaney GP, Lehmann J, Sykes J, Vinod S, Holloway L. Infrastructure platform for privacy-preserving distributed machine learning development of computer-assisted theragnostics in cancer. Journal of Biomedical Informatics 2022 Oct;134:104181. doi: https://10.1016/j.jbi.2022.104181.

[3] Field M, Vinod S, Delaney G, Aherne N, Bailey M, Carolan M, Dekker A, Greenham S, Hau E, Lehmann J, Ludbrook J, Miller A, Rezo A, Selvaraj J, Sykes J, Holloway L, Thwaites D, Implementation of the Australian computer-assisted theragnostics (AusCAT) network for radiation oncology data extraction, reporting and distributed learning, Journal of Medical Imaging and Radiation Oncology 65 (5) (2021) 627–636. doi:https://doi.org/10.1111/1754-9485.13287.
