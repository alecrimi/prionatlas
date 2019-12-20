#!/usr/bin/env python3
import os
import SimpleITK as sitk
import numpy as np
from scipy.stats import t


def average_volumes(fnames, out_path ): #, spacing=[0.025, 0.025, 0.025]
    reader = sitk.ImageFileReader()
    images = []
    for f in fnames:
        print f
        reader.SetFileName(f)
        volume = reader.Execute()
        volume_flat = sitk.GetArrayFromImage(volume)
        dims = volume_flat.shape
        volume_flat = volume_flat.flatten()
        images.append(volume_flat.astype(np.float))
    reference_flat = np.mean(images,axis=0)
    reference_flat = reference_flat.astype(np.uint16)
    avg_vol = sitk.GetImageFromArray(reference_flat.reshape(dims))
    avg_vol.SetSpacing(spacing)
    sitk.WriteImage(avg_vol, out_path)

 
def iterative_affine_reg(fnames):
    niter = 1
    print(niter)
    for i in xrange(niter):
       # print "-->\tAFFINE ITERATION", i+1, "/", niter
        if i == 0:
            print("init")
            ref_name = "bootstrap_ref.tif"
        else:
            ref_name = "affine_ref.tif" 
        ref = os.path.join(final_out_dir, ref_name)
        k = 1
        for f in fnames:
          if f.endswith(".tif"): 
            #print "-->\t[",k,"/",len(fnames),"]:", f
            fname = os.path.join(data_dir, f)
            out_dir = os.path.join(reg_dir, "reg_" + f)
            if not os.path.exists(out_dir): 
                os.makedirs(out_dir)
            print "-->\t Changing to", out_dir
            cmd = "elastix -f %s -m %s -out %s -p params/rigid.txt -p params/affine.txt > /tmp/out.txt 2>&1" % (ref, fname, out_dir)
            os.system(cmd)
            k += 1
        next_ref_name = "affine_ref.tif" 
        affine_fnames = [os.path.join(reg_dir, "reg_"+f, "result.1.tif") for f in fnames]
        average_volumes(affine_fnames, os.path.join(final_out_dir, next_ref_name))


def iterative_warp_reg(fnames, restart=True):
    niter = 2
    print("Warp")
    for i in xrange(niter):
        print "-->\tWARP ITERATION", i+1, "/", niter
        if i == 0 and restart:
            ref_name = "affine_ref.tif"
        else: 
            ref_name = "warp_ref.tif" 
        #print "Using as ref:", ref_name
        ref = os.path.join(final_out_dir, ref_name)
        k = 1
        for f in fnames:
          if f.endswith(".tif"): 
            print "-->\t[",k,"/",len(fnames),"]:", f
            #fname = os.path.join(data_dir, f)
            fname = os.path.join(reg_dir,"reg_" + f,"result.1.tif")

            out_dir = os.path.join(reg_dir, "reg_" + f)
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            print "-->\t Changing to", out_dir
            cmd = "elastix -f %s -m %s -out %s -p   params/Parameters_BSpline.txt > /tmp/out.txt 2>&1" % (ref, fname, out_dir)
            os.system(cmd) 
            fname_old = os.path.join(reg_dir,"reg_" + f,"result.0.tif")
            fname2 =    os.path.join(reg_dir,"reg_" + f,"result.2.tif")
            rename = "mv "+fname_old +" " +fname2
            os.system(rename)

            k += 1
        next_ref_name = "warp_ref.tif" 
        warp_fnames = [os.path.join(reg_dir, "reg_"+f, "result.2.tif") for f in fnames]
        average_volumes(warp_fnames, os.path.join(final_out_dir, next_ref_name))


################################### MAIN BLOCK ###############################

ref = "template_25_half_set_right_border.tif"
data_dir = "downsamp_tiff"
reg_dir = "reg"
out_dir = "output"
final_out_dir = "final_output"
pixel_size = 8*(3.3/1000.) #4*(5.16/1000.) # in mm
pixel_size_z = 8*(3/1000.) # in mm
spacing=[pixel_size, pixel_size, pixel_size_z]
reader = sitk.ImageFileReader()

'''
# Scale volumes
 
# Convert to nifti
k = 1
for file in os.listdir(data_dir):
    if file.endswith(".tif"):

        print(os.path.join(data_dir, file)) 
        fname = os.path.join(data_dir, file)

        out_dir = os.path.join(reg_dir, "reg_" + file)

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        cmd = "elastix -f %s -m %s -out %s -p params/rigid.txt > /tmp/out.txt 2>&1" % (ref, fname, out_dir)
        os.system(cmd)
        k += 1

# average together
#rigid_fnames = [os.path.join(reg_dir, "reg_"+f, "result.0.nrrd") for f in data_dir]
#average_volumes_old(rigid_fnames, os.path.join(final_out_dir, "bootstrap_ref.nrrd"))


images = []
for file in os.listdir(data_dir):
   if file.endswith(".tif"):
        f = os.path.join(reg_dir , "reg_"+file, "result.0.tif")
        #f = os.path.join(reg_dir , "reg_"+file)

        reader.SetFileName(f)
        volume = reader.Execute()
        volume_flat = sitk.GetArrayFromImage(volume)
        dims = volume_flat.shape
        print file
        print dims
        volume_flat = volume_flat.flatten()
        images.append(volume_flat.astype(np.float))
        print np.shape(images)
reference_flat = np.mean(images,axis=0)
reference_flat = reference_flat.astype(np.uint16)
avg_vol = sitk.GetImageFromArray(reference_flat.reshape(dims))
avg_vol.SetSpacing(spacing)
next_ref_name = "bootstrap_ref.tif" 
out_path = os.path.join(final_out_dir, next_ref_name)
sitk.WriteImage(avg_vol, out_path)
 
fnames = os.listdir(data_dir) # filenames to be included

# perform iterative affine registration
#iterative_affine_reg(fnames)

# perform iterative warping registration
iterative_warp_reg(fnames, restart=True)
''' 
#########################################
# Create other mean

#print os.listdir("temp")
temp1 = "res/nbs" 
temp2 = "res/RML"

images_pr = []
for file in os.listdir(data_dir):
   if file.endswith("RML_sc.tif"):
        f = os.path.join(reg_dir , "reg_"+file, "result.2.tif")
        print f
        reader.SetFileName(f)
        volume = reader.Execute()
        volume_flat = sitk.GetArrayFromImage(volume)
        dims = volume_flat.shape
        volume_flat = volume_flat.flatten()
        images_pr.append(volume_flat.astype(np.float))  
reference_flat = np.mean(images_pr,axis=0)
mean_prion = reference_flat.astype(np.uint16)  
row,col =  np.shape(images_pr)

images_nc = []
for file in os.listdir(data_dir):
   if file.endswith("NBH.tif"):
        f = os.path.join(reg_dir , "reg_"+file, "result.2.tif")
        print f
        reader.SetFileName(f)
        volume = reader.Execute()
        volume_flat = sitk.GetArrayFromImage(volume)
        dims = volume_flat.shape
        volume_flat = volume_flat.flatten()
        images_nc.append(volume_flat.astype(np.float))  
reference_flat = np.mean(images_nc,axis=0)
mean_control = reference_flat.astype(np.uint16)  
row2,col =  np.shape(images_nc)


sigma = np.sqrt( np.var(images_pr)/row +  np.var(images_nc)/col   ) 
 
ttest =   (mean_prion - mean_control ) / sigma
n = row + row2
pval =  t.sf(np.abs(ttest), n-1)*2 
reference_flat = pval.astype(np.uint16)
pval_vol = sitk.GetImageFromArray(reference_flat.reshape(dims))
pval_vol.SetSpacing(spacing)
 
sitk.WriteImage(pval_vol, "result_pval.tif")
