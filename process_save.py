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

#from tqdm import tqdm
s = input("file path: ")
if os.path.isdir(s):
    filelist = sorted(glob.glob(s+'*.nc'))
else:
    print("Directory does not exist")
CAMPAIGN = 'DYNAMO_DiegoGarcia'

if CAMPAIGN=='TAHOPE':
    data = []
    for i in range(len(filelist)):
        data.append(xarray.open_dataset(filelist[i]))
        
    time_utc = []
    backscatter = []
    cbh_1 = []    
    level = data[0]['level']
    for i in range(len(filelist)):
        time_utc.append(data[i]['time_utc'])
        backscatter.append(data[i]['backscatter'])
        cbh_1.append(data[i]['cbh_1'])
    cbh1_smpl = np.array([])
    time_smpl = np.array([])
    backscatter_smpl = np.concatenate(backscatter,0)
    
    for i in range(len(filelist)):
        cbh1_smpl = np.concatenate([cbh1_smpl,np.asarray(cbh_1[i])],axis=0)
        time_smpl = np.concatenate([time_smpl,np.asarray(time_utc[i])],axis=0)
    
    def process_data(scatterer=None):
        import scipy.signal
        CWT0_5_0610 = np.zeros((720*8,4))
        CWT18_24_0610 = np.zeros((720*6,4))
        hs_cl31_smoo0610 = scipy.signal.wiener(np.asarray(scatterer[:,0:600]),9)
        hs_cl31_smoo0610n = hs_cl31_smoo0610[720*18:720*24,0:600]
        gradient0610 = peak_finding.gradient_method(hs_cl31_smoo0610[0:720*8,0:600],np.linspace(0,599,600)).gradient_peak()
        gradient0610n = peak_finding.gradient_method(hs_cl31_smoo0610n[:,0:600],np.linspace(0,599,600)).gradient_peak()
        
        for i in tqdm(range(720*8)):
            haar_CWT,temp = peak_finding.CWT_transform(hs_cl31_smoo0610[i,0:600]).CWT_peaks(a=40)
            if (temp.shape[0]<4) and (temp.shape[0]>0):
                CWT0_5_0610[i,0:temp.shape[0]] = temp[0:temp.shape[0]]
                CWT0_5_0610[i,temp.shape[0]:4] = temp[-1]+10
            elif (temp.shape[0]==0):
                CWT0_5_0610[i,0:4] = np.nan
            else:    
                CWT0_5_0610[i,:] = temp[0:4]
            del haar_CWT,temp
            
        for i in tqdm(range(720*6-1)):
            haar_CWT,temp = peak_finding.CWT_transform(hs_cl31_smoo0610n[i,0:600]).CWT_peaks(a=40)
            if (temp.shape[0]<4) and (temp.shape[0]>0):
                CWT18_24_0610[i,0:temp.shape[0]] = temp[0:temp.shape[0]]
                CWT18_24_0610[i,temp.shape[0]:4] = temp[-1]+10
            elif (temp.shape[0]==0):
                CWT18_24_0610[i,0:4] = np.nan
            else:    
                CWT18_24_0610[i,:] = temp[0:4]
            del haar_CWT,temp
        return CWT0_5_0610,CWT18_24_0610,gradient0610,gradient0610n,hs_cl31_smoo0610,hs_cl31_smoo0610n
    
    CWTpeak_morning = []
    CWTpeak_night = []
    gpeak_morning = []
    gpeak_night = []
    scatterer_morning = []
    scatterer_night = []
    
    for i in range(len(filelist)):
        temp1,temp2,temp3,temp4,temp5,temp6 = process_data(scatterer=backscatter[i])
        CWTpeak_morning.append(temp1)
        CWTpeak_night.append(temp2)
        gpeak_morning.append(temp3)
        gpeak_night.append(temp4)
        scatterer_morning.append(temp5)
        scatterer_night.append(temp6)
    
    #------------------
    # write netCDF file
    #------------------
    fileNAME = ['0609','0610','0611']
    
    for i in range(len(filelist)):
        # open a netCDF file to write
        ncout = Dataset('/scra6/ft21894/py_script/TAHOPE/obs/NTU_CL31/dryrun/processed/TAHOPE_CL31_'+\
                        str(fileNAME[i])+'.nc', 'w', format='NETCDF4')
        # define axis size
        ncout.createDimension('time_morning', None)  # unlimited
        ncout.createDimension('time_night', None)  # unlimited
        ncout.createDimension('scat_lvl', 600)
        ncout.createDimension('CWTp_lvl', 4)
        ncout.createDimension('gp_lvl', 8)

        # create variable array
        CWTpeak_morningout = ncout.createVariable('CWTpeak_morning', dtype('double').char, ('time_morning', 'CWTp_lvl'))
        CWTpeak_nightout = ncout.createVariable('CWTpeak_night', dtype('double').char, ('time_night', 'CWTp_lvl'))
        gpeak_morningout = ncout.createVariable('gpeak_morning', dtype('double').char, ('time_morning', 'gp_lvl'))
        gpeak_nightout = ncout.createVariable('gpeak_night', dtype('double').char, ('time_night', 'gp_lvl'))
        scatterer_morningout = ncout.createVariable('scatterer_morning', dtype('double').char, ('time_morning', 'scat_lvl'))
        scatterer_nightout = ncout.createVariable('scatterer_night', dtype('double').char, ('time_night', 'scat_lvl'))
        
        # copy axis from original dataset
        CWTpeak_morningout[:] = np.asarray(CWTpeak_morning[i][:])
        CWTpeak_nightout[:] = np.asarray(CWTpeak_night[i][:])
        gpeak_morningout[:] = np.asarray(gpeak_morning[i][:])
        gpeak_nightout[:] = np.asarray(gpeak_night[i][:])
        scatterer_morningout[:] = np.asarray(scatterer_morning[i][:])
        scatterer_nightout[:] = np.asarray(scatterer_night[i][:])
        
        # close files
        ncout.close()

if CAMPAIGN=='CPS08':
    data = []
    for i in range(len(filelist)):
        data.append(xarray.open_dataset(filelist[i]))
        
    backscatter = []
    cbh_1 = []    
    for i in range(len(filelist)):
        backscatter.append(data[i]['backscatter_out'])
        cbh_1.append(data[i]['cbh1_out'])
    cbh1_smpl = np.array([])
    backscatter_smpl = np.concatenate(backscatter,0)
    
    for i in range(len(filelist)):
        cbh1_smpl = np.concatenate([cbh1_smpl,np.asarray(cbh_1[i])],axis=0)
    
    def process_data(scatterer=None):
        import scipy.signal
        CWT0_5_0610 = np.zeros((225*18,4))
        hs_cl31_smoo0610 = scipy.signal.wiener(np.asarray(scatterer[0:225*18,0:200]),9)
        gradient0610 = peak_finding.gradient_method(hs_cl31_smoo0610[0:225*18,0:200],np.linspace(0,199,200)).gradient_peak()
        
        for i in tqdm(range(225*18)):
            haar_CWT,temp = peak_finding.CWT_transform(hs_cl31_smoo0610[i,0:200]).CWT_peaks(a=40)
            if (temp.shape[0]<4) and (temp.shape[0]>0):
                CWT0_5_0610[i,0:temp.shape[0]] = temp[0:temp.shape[0]]
                CWT0_5_0610[i,temp.shape[0]:4] = temp[-1]+10
            elif (temp.shape[0]==0):
                CWT0_5_0610[i,0:4] = np.nan
            else:    
                CWT0_5_0610[i,:] = temp[0:4]
            del haar_CWT,temp
        return CWT0_5_0610,gradient0610,hs_cl31_smoo0610
    
    CWTpeak = []
    gpeak = []
    scatterer = []
    
    for i in range(len(filelist)):
        temp1,temp2,temp3 = process_data(scatterer=backscatter[i])
        CWTpeak.append(temp1)
        gpeak.append(temp2)
        scatterer.append(temp3)
    
    #datename = ['0801','0802','0803','0806','0807','0808','0809','0810','0811','0812','0813','0814','0816','0817','0818','0819',\
    #            '0820','0821','0822','0823','0824','0825','0826','0827','0828','0829','0830']
    datename = ['0901','0902','0903','0904','0905','0906',\
                '0909','0910','0911','0912','0914','0915','0916','0917','0918','0919',\
                '0920','0921','0922','0923','0924','0926','0927','0928','0929']
    
    for i in range(len(filelist)):
        # open a netCDF file to write
        ncout = Dataset('/scra6/ft21894/py_script/TAHOPE/obs/CPS08_ceilometer/processed/CPS08_CL31_proc_'+\
                        str(datename[i])+'.nc', 'w', format='NETCDF4')
        # define axis size
        ncout.createDimension('time', None)  # unlimited
        ncout.createDimension('scat_lvl', 200)
        ncout.createDimension('CWTp_lvl', 4)
        ncout.createDimension('gp_lvl', 8)

        # create variable array
        CWTpeak_out = ncout.createVariable('CWTpeak_out', dtype('double').char, ('time', 'CWTp_lvl'))
        gpeak_out = ncout.createVariable('gpeak_out', dtype('double').char, ('time', 'gp_lvl'))
        scatterer_out = ncout.createVariable('scatterer_out', dtype('double').char, ('time', 'scat_lvl'))
        
        # copy axis from original dataset
        CWTpeak_out[:] = np.asarray(CWTpeak[i][:])
        gpeak_out[:] = np.asarray(gpeak[i][:])
        scatterer_out[:] = np.asarray(scatterer[i][:])
        
        # close files
        ncout.close()

if CAMPAIGN=='PECAN_SPol':
    data = []
    for i in range(len(filelist)):
        data.append(xarray.open_dataset(filelist[i],decode_times=False))
        
    backscatter = []
    cbh_1 = []    
    for i in range(len(filelist)):
        backscatter.append(data[i]['backscatter'])
        cbh_1.append(data[i]['cbh'])
    
    def process_data(scatterer=None):
        import scipy.signal
        CWT0_5_0610 = np.zeros((240*19,4))
        hs_cl31_smoo0610 = scipy.signal.wiener(np.asarray(scatterer[0:240*19,0:300]),9)
        gradient0610 = peak_finding.gradient_method(hs_cl31_smoo0610[0:240*19,0:300],np.linspace(0,299,300)).gradient_peak()
        
        for i in tqdm(range(240*19)):
            haar_CWT,temp = peak_finding.CWT_transform(hs_cl31_smoo0610[i,0:300]).CWT_peaks(a=40)
            if (temp.shape[0]<4) and (temp.shape[0]>0):
                CWT0_5_0610[i,0:temp.shape[0]] = temp[0:temp.shape[0]]
                CWT0_5_0610[i,temp.shape[0]:4] = temp[-1]+10
            elif (temp.shape[0]==0):
                CWT0_5_0610[i,0:4] = np.nan
            else:    
                CWT0_5_0610[i,:] = temp[0:4]
            del haar_CWT,temp
        return CWT0_5_0610,gradient0610,hs_cl31_smoo0610
    
    CWTpeak = []
    gpeak = []
    scatterer = []
    
    for i in range(len(filelist)):
        temp1,temp2,temp3 = process_data(scatterer=backscatter[i])
        CWTpeak.append(temp1)
        gpeak.append(temp2)
        scatterer.append(temp3)
    
    #datename = ['0602','0603','0604','0605','0606','0608','0609','0610','0611','0612','0613','0614','0616','0617','0618','0619',\
    #            '0620','0621','0622','0624','0625','0626','0627','0628','0629','0630']
    datename = ['0701','0702','0703','0704','0705','0706','0707','0708','0709','0710','0711']
    
    for i in range(len(filelist)):
        # open a netCDF file to write
        ncout = Dataset('/scra6/ft21894/py_script/TAHOPE/obs/PECAN15_SPolCL31/processed/PECAN_SPolCL31_proc_'+\
                        str(datename[i])+'.nc', 'w', format='NETCDF4')
        # define axis size
        ncout.createDimension('time', None)  # unlimited
        ncout.createDimension('scat_lvl', 300)
        ncout.createDimension('CWTp_lvl', 4)
        ncout.createDimension('gp_lvl', 8)

        # create variable array
        CWTpeak_out = ncout.createVariable('CWTpeak_out', dtype('double').char, ('time', 'CWTp_lvl'))
        gpeak_out = ncout.createVariable('gpeak_out', dtype('double').char, ('time', 'gp_lvl'))
        scatterer_out = ncout.createVariable('scatterer_out', dtype('double').char, ('time', 'scat_lvl'))
        
        # copy axis from original dataset
        CWTpeak_out[:] = np.asarray(CWTpeak[i][:])
        gpeak_out[:] = np.asarray(gpeak[i][:])
        scatterer_out[:] = np.asarray(scatterer[i][:])
        
        # close files
        ncout.close()

def process_data(scatterer=None,fileshape=None,smooth=None,vert_res=None):
    import scipy.signal
    CWT0_5_0610 = np.zeros((fileshape[0],4))
    hs_cl31_smoo0610 = scipy.signal.wiener(np.asarray(scatterer[0:fileshape[0],0:fileshape[1]]),smooth)
    gradient0610= peak_finding.gradient_method(hs_cl31_smoo0610[0:fileshape[0],\
                                                                0:fileshape[1]],\
                                               np.linspace(0,fileshape[1]-1,fileshape[1]),vert_res).gradient_peak()
        
    for i in tqdm(range(fileshape[0])):
        haar_CWT,temp = peak_finding.CWT_transform(hs_cl31_smoo0610[i,0:fileshape[1]],vert_res).CWT_peaks(a=40)
        if (temp.shape[0]<4) and (temp.shape[0]>0):
            CWT0_5_0610[i,0:temp.shape[0]] = temp[0:temp.shape[0]]
            CWT0_5_0610[i,temp.shape[0]:4] = temp[-1]+10
        elif (temp.shape[0]==0):
            CWT0_5_0610[i,0:4] = np.nan
        else:    
            CWT0_5_0610[i,:] = temp[0:4]
        del haar_CWT,temp
    return CWT0_5_0610,gradient0610,hs_cl31_smoo0610
    
if CAMPAIGN=='DYNAMO_Gan':
    data = []
    for i in range(len(filelist)):
        data.append(xarray.open_dataset(filelist[i],decode_times=False))
    backscatter = []
    cbh_1 = []    
    for i in range(len(filelist)):
        backscatter.append(data[i]['backscatter'])
        cbh_1.append(data[i]['first_cbh'])
    CWTpeak = []
    gpeak = []
    scatterer = []
    for i in range(len(filelist)):
        temp1,temp2,temp3 = process_data(scatterer=backscatter[i],fileshape=[225*24,250],smooth=10,vert_res=30)
        CWTpeak.append(temp1)
        gpeak.append(temp2)
        scatterer.append(temp3)
    
    datename = ['1002','1003','1004','1005','1006','1007','1008','1009','1010','1011','1012','1013','1014',\
                '1103','1104','1105','1106','1107','1108','1109','1110','1111','1112','1113','1114','1115','1116','1117',\
                '1204','1205','1206','1207',\
                '1225','1226','1227','1228','1229','1230']
    
    for i in range(len(filelist)):
        # open a netCDF file to write
        ncout = Dataset('/scra6/ft21894/py_script/TAHOPE/Data/DYNAMO/DYNAMO_GanCeliometer/processed/DYNAMO_Gan_proc_'+\
                        str(datename[i])+'.nc', 'w', format='NETCDF4')
        # define axis size
        ncout.createDimension('time', None)  # unlimited
        ncout.createDimension('scat_lvl', 250)
        ncout.createDimension('CWTp_lvl', 4)
        ncout.createDimension('gp_lvl', 8)

        # create variable array
        CWTpeak_out = ncout.createVariable('CWTpeak_out', dtype('double').char, ('time', 'CWTp_lvl'))
        gpeak_out = ncout.createVariable('gpeak_out', dtype('double').char, ('time', 'gp_lvl'))
        scatterer_out = ncout.createVariable('scatterer_out', dtype('double').char, ('time', 'scat_lvl'))
        
        # copy axis from original dataset
        CWTpeak_out[:] = np.asarray(CWTpeak[i][:])
        gpeak_out[:] = np.asarray(gpeak[i][:])
        scatterer_out[:] = np.asarray(scatterer[i][:])
        
        # close files
        ncout.close()

if CAMPAIGN=='DYNAMO_DiegoGarcia':
    data = []
    for i in range(len(filelist)):
        data.append(xarray.open_dataset(filelist[i],decode_times=False))
    backscatter = []
    cbh_1 = []    
    for i in range(len(filelist)):
        backscatter.append(data[i]['backscatter_out'])
        cbh_1.append(data[i]['cbh1_out'])
    CWTpeak = []
    gpeak = []
    scatterer = []
    for i in range(len(filelist)):
        temp1,temp2,temp3 = process_data(scatterer=backscatter[i],fileshape=[223*24,300],smooth=7,vert_res=10)
        CWTpeak.append(temp1)
        gpeak.append(temp2)
        scatterer.append(temp3)
    
    datename = ['1003','1005','1006','1017','1018','1019','1020','1021','1022','1023','1024','1028','1029','1030',\
                '1105','1106','1110','1111','1112','1114','1126','1127','1130',\
                '1205','1206','1209','1210','1211','1212','1214','1215','1216',\
                '1223','1224','1227','1228','1229','1230','1231']
    
    for i in range(len(filelist)):
        # open a netCDF file to write
        ncout = Dataset('/scra6/ft21894/py_script/TAHOPE/Data/DYNAMO/DYNAMO_DiegoGarcia_Ceilometer/processed/DYNAMO_DG_proc_'+\
                        str(datename[i])+'.nc', 'w', format='NETCDF4')
        # define axis size
        ncout.createDimension('time', None)  # unlimited
        ncout.createDimension('scat_lvl', 300)
        ncout.createDimension('CWTp_lvl', 4)
        ncout.createDimension('gp_lvl', 8)

        # create variable array
        CWTpeak_out = ncout.createVariable('CWTpeak_out', dtype('double').char, ('time', 'CWTp_lvl'))
        gpeak_out = ncout.createVariable('gpeak_out', dtype('double').char, ('time', 'gp_lvl'))
        scatterer_out = ncout.createVariable('scatterer_out', dtype('double').char, ('time', 'scat_lvl'))
        
        # copy axis from original dataset
        CWTpeak_out[:] = np.asarray(CWTpeak[i][:])
        gpeak_out[:] = np.asarray(gpeak[i][:])
        scatterer_out[:] = np.asarray(scatterer[i][:])
        
        # close files
        ncout.close()
