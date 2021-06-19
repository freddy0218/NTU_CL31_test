import numpy as np
import warnings
warnings.filterwarnings("ignore",category=RuntimeWarning)
from netCDF4 import Dataset
from scipy.ndimage.filters import gaussian_filter
import scipy
import glob,os
import xarray
from numpy import dtype
from tqdm import tqdm

import ceilometer_help
import peak_finding


s = input("file path: ")
if os.path.isdir(s):
    filelist = sorted(glob.glob(s+'*'))
else:
    print("Directory does not exist")
    
CAMPAIGN = 'DYNAMO'#'CPS08'

if CAMPAIGN=='CPS08':
    datename = ['0901','0902','0903','0904','0905','0906',\
                '0909','0910','0911','0912','0914','0915','0916','0917','0918','0919',\
                '0920','0921','0922','0923','0924','0926','0927','0928','0929']
    for i in tqdm(range(len(datename))):
        tmp_file = glob.glob(s+'*'+datename[i]+'*')
        
        backscatterCPS = []
        cbh_1CPS = []     
        for j in range(len(tmp_file)):
            tmp_dataset = xarray.open_dataset(tmp_file[j])
            backscatterCPS.append(tmp_dataset['backscatter'])
            cbh_1CPS.append(tmp_dataset['cbh_1'])
        
        cbh1CPS_smpl = np.array([])
        backscatterCPS_smpl = np.concatenate(backscatterCPS,0)       
        for j2 in range(len(tmp_file)):
            cbh1CPS_smpl = np.concatenate([cbh1CPS_smpl,np.asarray(cbh_1CPS[j2])],axis=0)
        
        # open a netCDF file to write
        ncout = Dataset('/scra6/ft21894/py_script/TAHOPE/obs/CPS08_ceilometer/processed/CPS08_CL31_'+\
                        str(datename[i])+'.nc', 'w', format='NETCDF4')
        # define axis size
        ncout.createDimension('time', None)  # unlimited
        ncout.createDimension('vert_lv', 385)

        # create variable array
        backscatter_out = ncout.createVariable('backscatter_out', dtype('double').char, ('time', 'vert_lv'))
        cbh1_out = ncout.createVariable('cbh1_out', dtype('double').char, ('time'))
        
        # copy axis from original dataset
        backscatter_out[:] = np.asarray(backscatterCPS_smpl[:])
        cbh1_out[:] = np.asarray(cbh1CPS_smpl[:])
        
        del tmp_file, tmp_dataset, backscatterCPS, cbh_1CPS, cbh1CPS_smpl, backscatterCPS_smpl
        # close files
        ncout.close()
    
elif CAMPAIGN=='DYNAMO':
    datename = ['1002','1003','1004','1005','1006',\
                '1009','1010','1011','1012','1014','1015','1016','1017','1018','1019',\
                '1020','1021','1022','1023','1024','1026','1027','1028','1029','1030','1031',\
                '1101','1102','1103','1104','1105','1106',\
                '1109','1110','1111','1112','1114','1115','1116','1117','1118','1119',\
                '1120','1121','1122','1123','1124','1126','1127','1128','1129','1130',\
                '1201','1202','1203','1204','1205','1206',\
                '1209','1210','1211','1212','1214','1215','1216','1217','1218','1219',\
                '1220','1221','1222','1223','1224','1226','1227','1228','1229','1230','1231']
    for i in tqdm(range(len(datename))):
        tmp_file = glob.glob(s+'A1'+datename[i]+'*')
        print(tmp_file)
        backscatterCPS = []
        cbh_1CPS = []
        for j in range(len(tmp_file)):
            tmp_dataset = xarray.open_dataset(tmp_file[j])
            backscatterCPS.append(tmp_dataset['backscatter'])
            cbh_1CPS.append(tmp_dataset['cbh_1'])
        
        cbh1CPS_smpl = np.array([])
        backscatterCPS_smpl = np.concatenate(backscatterCPS,0)
        for j2 in range(len(tmp_file)):
            cbh1CPS_smpl = np.concatenate([cbh1CPS_smpl,np.asarray(cbh_1CPS[j2])],axis=0)
            
        # open a netCDF file to write
        ncout = Dataset('/scra6/ft21894/py_script/TAHOPE/Data/DYNAMO_DiegoGarcia_Ceilometer/processed/DYNAMO_DG_CL31_'+\
                        str(datename[i])+'.nc', 'w', format='NETCDF4')
        # define axis size
        ncout.createDimension('time', None)  # unlimited
        ncout.createDimension('vert_lv', 770)
        # create variable array
        backscatter_out = ncout.createVariable('backscatter_out', dtype('double').char, ('time', 'vert_lv'))
        cbh1_out = ncout.createVariable('cbh1_out', dtype('double').char, ('time'))
        # copy axis from original dataset
        backscatter_out[:] = np.asarray(backscatterCPS_smpl[:])
        cbh1_out[:] = np.asarray(cbh1CPS_smpl[:])
        del tmp_file, tmp_dataset, backscatterCPS, cbh_1CPS, cbh1CPS_smpl, backscatterCPS_smpl
        # close files
        ncout.close()

        
