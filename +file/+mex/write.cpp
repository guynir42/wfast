#include "mex.h"
#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <thread>
#include <chrono>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "SaveData.h"
#include "SaveDataHDF5.h"
#include "MyMatrix.h"
#include "MyAttribute.h"


const char *parse_extension(const char *filename);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){
					  
	if(nrhs<1){ // usage is detailed when function is called without arguments!
		SaveData::print_help();
		return;
	}
	
	// check minimal inputs
	if(nrhs<2) mexErrMsgIdAndTxt("MATLAB:file:mex:mexWrite:invalidNumInputs", "At least 2 arguments are required!");
	if(mxIsChar(prhs[0])!=1) mexErrMsgIdAndTxt( "MATLAB:file:mex:mexWrite:inputNotString", "filename must be a string.");
	if(mxGetM(prhs[1])!=1) mexErrMsgIdAndTxt( "MATLAB:file:mex:mexWrite:inputNotVector", "mex_flag must be a row vector.");
	if(mxGetN(prhs[1])!=2) mexErrMsgIdAndTxt( "MATLAB:file:mex:mexWrite:inputSizeMismatch", "mex_flag must have exactly 2 elements!");
	
	SaveData *s=0;
	const char *filename=mxArrayToString(prhs[0]);
	const char *extension=parse_extension(filename);
	
	if (strlen(extension)==0 || SaveData::cs(extension, "h5", "h5z", "hdf5")) s=new SaveDataHDF5(filename);
	else mexErrMsgIdAndTxt( "MATLAB:file:mex:mexWrite:unknownExtension", "filename extension is unknown. Try 'h5' or 'hdf5'...");
	
	s->mex_flag=mxGetPr(prhs[1]);
	s->parseVararginPairs(nrhs, prhs);
		
	if(s->debug_bit>2){// get back some feedback on inputs...

		mexPrintf("\nfilename: '%s'\n", s->filename);
		mexPrintf("mex_flag_write= %g %g\n", s->mex_flag[0], s->mex_flag[1]);
		
		s->images.printout();
		s->images_raw.printout();
		s->images_cal.printout();
		s->cutouts_raw.printout();
		s->cutouts_cal.printout();
		s->positions.printout();
		s->full_sum.printout();
		s->timestamps.printout();
		s->psfs.printout();
		s->lightcurves.printout();
		
		mexPrintf("debug= %d | deflate= %d | chunk= %d | async= %d \n\n", s->debug_bit, s->deflate, s->chunk_size, s->async_write);
		
	}
	
	if(s->debug_bit>5){ // print all the parameters in the cell array too
		
		for(int i=0;i<s->parameter_addresses_vector.size();i++){
			
			printf("METADATA OBJECT IN: %s\n", s->parameter_addresses_vector[i].c_str());
			
			for(int j=0;j<s->parameter_attributes_2D_vector[i].size(); j++){
				
				s->parameter_attributes_2D_vector[i][j].printout();
				
			}
			
		}
		
	}
	
	if(s->async_write==0) s->write();
	else{
		std::thread mythread(&SaveData::write, s);
		mythread.detach();
	}
	
	
}

const char *parse_extension(const char *filename){
	
	const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
	
}
	
/*
std::vector<Attribute> getAttributes(const mxArray *pars){
	
	std::vector<Attribute> vec=std::vector<Attribute>();
	
	int N=mxGetNumberOfFields(pars);
	
	for(int i=0;i<N;i++){
		
		const char *name=mxGetFieldNameByNumber(pars, i);
		mxArray *value=mxGetFieldByNumber(pars, 0, i);
		
		int is_empty=mxIsEmpty(value);
		int is_numeric=mxIsNumeric(value);
		int is_char=mxIsChar(value);
				
		if(!is_numeric && !is_char) continue; // skip objects and other classes
		if(mxGetN(value)>1 && mxGetM(value)>1) continue; // skip matrix attributes...
		if(is_empty) continue; // skip empty attributes
		
		hsize_t rows=mxGetN(value);
		hsize_t cols=mxGetM(value);
		int is_vec=!mxIsScalar(value);
		
		Attribute a;
		a.name=std::string(name);
		
		if(is_char){
			a.str=std::string(mxArrayToString(value));
			a.is_str=1;	
			vec.push_back(a);
		}
		else if(is_vec){
			size_t len=mxGetN(value);
			if(len==1) len=mxGetM(value);
			a.vec.assign(mxGetPr(value), mxGetPr(value)+len);
			a.is_vec=1;
			vec.push_back(a);
		}
		else if(is_numeric){
			a.scalar=mxGetScalar(value);
			a.is_scalar=1;
			vec.push_back(a);
		}
				
	}
	
	return vec;
	
}

void writeHDF5(InputData data){
	
	data.mex_flag[0]=1;
	data.mex_flag[1]=0;	
	data.mex_flag[2]=0;
	
	if(data.debug_bit>3) mexPrintf("began writing\n");
	
	hid_t file_id, dataset_id, dataspace_id, plist_id, att_dataspace_id;
	file_id = H5Fcreate(data.filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
	int status;
	
	///////////////// write images //////////////////////////////
	if(data.images_length){ 
	
		if(data.debug_bit>9) mexPrintf("writing images...\n");
	
		hsize_t img_dims_c[4]={data.num_cutouts, data.num_frames, data.num_cols, data.num_rows}; // c style array with a flipped order of dims (c style!)
		
		dataspace_id = H5Screate_simple(4,img_dims_c, NULL);
		if(dataspace_id<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:create_dataspace", "Something wrong when creating dataspace for images!"); }
		
		plist_id  = H5Pcreate(H5P_DATASET_CREATE);
		
		if(data.deflate){
			
			if(data.debug_bit>2) mexPrintf("deflating data at level %d\n", data.deflate);
			
			hsize_t chunk_dims[4]={1, 1, (hsize_t) data.chunk_size, (hsize_t) data.chunk_size};
			
			for(int i=0;i<4;i++) if(chunk_dims[i]>img_dims_c[i]) chunk_dims[i]=img_dims_c[i]; // make sure chunk size is no bigger than actual size of image input...
			
			status = H5Pset_chunk (plist_id, 4, chunk_dims);
			if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:set_chunk", "Something wrong when setting chunk for images!"); }
			
			status = H5Pset_deflate (plist_id, data.deflate);
			if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:set_deflate", "Something wrong when setting deflate! for images"); }
			
		}
		if(data.images_int){// write uint16 data
			dataset_id = H5Dcreate(file_id, "/images", H5T_NATIVE_USHORT, dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT); // can use H5T_NATIVE_INT / H5T_NATIVE_USHORT
			if(dataset_id<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:create_dataset", "Something wrong when creating dataset for images (int16)!"); }
			
			status = H5Dwrite(dataset_id, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.images_int);
			if(data.debug_bit>2) mexPrintf("IMAGES: file_id= %d | dataspace_id= %d | dataset_id= %d | status= %d\n", file_id, dataspace_id, dataset_id, status);
			if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:write_dataset", "Something wrong when writing dataset for images (int16)!"); }
		}
		else if(data.images_double){// write double data
			dataset_id = H5Dcreate(file_id, "/images", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT); // can use H5T_NATIVE_INT / H5T_STD_U16BE
			if(dataset_id<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:create_dataset", "Something wrong when creating dataset for images (double)!"); }
			
			status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.images_double);
			if(data.debug_bit>2) mexPrintf("IMAGES: file_id= %d | dataspace_id= %d | dataset_id= %d | status= %d\n", file_id, dataspace_id, dataset_id, status);
			if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:write_dataset", "Something wrong when writing dataset for images (double)!"); }
		}
		
		H5Pclose(plist_id);
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		
    }
	
	///////////////// write cut pos /////////////////////////////////
	if(data.cut_pos_length){
		
		if(data.debug_bit>9) mexPrintf("writing cut_pos...\n");
		
		hsize_t cut_pos_dims[2] = {2, data.cut_pos_num_rows}; // c style dimensions are reversed! 
		dataspace_id = H5Screate_simple(2, cut_pos_dims, NULL);
		if(dataspace_id<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:create_dataspace", "Something wrong when creating dataspace for cut_pos!"); }
		
		dataset_id = H5Dcreate(file_id, "/cut_pos", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); // are you sure we don't need 64 bit == long int?
		if(dataset_id<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:create_dataset", "Something wrong when creating dataset for cut_pos!"); }
				
		status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.cut_pos);
		if(data.debug_bit>2) mexPrintf("CUT_POS: file_id= %d | dataspace_id= %d | dataset_id= %d | status= %d\n", file_id, dataspace_id, dataset_id, status);
		if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:write_dataset", "Something wrong when writing dataset for cut_pos!"); } 
		
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
				
	}
	
	///////////////// write timestamps //////////////////////////////
	if(data.timestamps_length){
	
		if(data.debug_bit>9) mexPrintf("writing timestamps...\n");
		
		hsize_t timestamps_dims_c[1]={data.timestamps_length};
		dataspace_id = H5Screate_simple(1, timestamps_dims_c, NULL);
		if(dataspace_id<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:create_dataspace", "Something wrong when creating dataspace for timestamps!"); }
			
		dataset_id = H5Dcreate(file_id, "/times", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); // can use H5T_NATIVE_INT / H5T_STD_U16BE 
		if(dataset_id<0){ data.mex_flag[2]=-1; 	mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:create_dataset", "Something wrong when creating dataset for timestamps!"); }
		
		status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.timestamps);
		if(data.debug_bit>2) mexPrintf("TIMESTAMPS: file_id= %d | dataspace_id= %d | dataset_id= %d | status= %d\n", file_id, dataspace_id, dataset_id, status);
		if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:write_dataset", "Something wrong when writing dataset for timestamps!");	} 
		
		hsize_t dims=1;
		hid_t attribute_id; 
		
		if(data.file_write_timestamp>=0){// ********* write the timestamps file_write attributes *********
		
			att_dataspace_id=H5Screate_simple(1, &dims, NULL);

			attribute_id=H5Acreate2(dataset_id, "file_write_timestamp", H5T_NATIVE_DOUBLE, att_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
			status=H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &data.file_write_timestamp);
			if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:write_timestamp_attributes", "Something wrong when writing file_write_timestamp attribute!"); }
			
			status = H5Aclose(attribute_id);
			if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:close_att", "Something wrong when closing file_write_timestamp attribute!"); }				
			
			status = H5Sclose(att_dataspace_id);
			if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:close_att_dataspace", "Something wrong when closing file_write_timestamp attribute dataspace!"); }
		
		}
		
		if(data.file_write_datestr.length()){// *********** file_write_datestr atribute ************
				
			att_dataspace_id=H5Screate_simple(1, &dims, NULL);

			hid_t atype = H5Tcopy(H5T_C_S1);
			H5Tset_size(atype, data.file_write_datestr.length());
			H5Tset_strpad(atype,H5T_STR_NULLTERM);
			
			attribute_id=H5Acreate2(dataset_id, "file_write_datestr", atype, att_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
			
			status=H5Awrite(attribute_id, atype, data.file_write_datestr.c_str());
			if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:write_timestamp_attributes", "Something wrong when writing file_write_datestr attribute!"); }
			
			status = H5Aclose(attribute_id);
			if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:close_att", "Something wrong when closing file_write_datestr attribute!"); }				
			
			status = H5Sclose(att_dataspace_id);
			if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:close_att_dataspace", "Something wrong when closing file_write_datestr attribute dataspace!"); }
			
		}
		
		if(data.file_start_datestr.length()){// *********** file_start_datestr atribute ************
				
			att_dataspace_id=H5Screate_simple(1, &dims, NULL);

			hid_t atype = H5Tcopy(H5T_C_S1);
			H5Tset_size(atype, data.file_start_datestr.length());
			H5Tset_strpad(atype,H5T_STR_NULLTERM);
			
			attribute_id=H5Acreate2(dataset_id, "file_start_datestr", atype, att_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
			
			status=H5Awrite(attribute_id, atype, data.file_start_datestr.c_str());
			if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:write_timestamp_attributes", "Something wrong when writing file_start_datestr attribute!"); }
			
			status = H5Aclose(attribute_id);
			if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:close_att", "Something wrong when closing file_start_datestr attribute!"); }				
			
			status = H5Sclose(att_dataspace_id);
			if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:close_att_dataspace", "Something wrong when closing file_start_datestr attribute dataspace!"); }
			
		}
		
		H5Pclose(plist_id);	
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		
	}
	
	/////////////// write PSFs //////////////////////
	if(data.psfs){ 
	
		if(data.debug_bit>9) mexPrintf("writing psfs...\n");
	
		hsize_t psfs_dims_c[4] = {data.psfs_num_cutouts, data.psfs_num_frames, data.psfs_num_cols, data.psfs_num_rows}; // c style array with a flipped order of dims (c style!)
				
		dataspace_id = H5Screate_simple(4,psfs_dims_c, NULL);
		if(dataspace_id<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:create_dataspace", "Something wrong when creating dataspace for PSFs!"); }
		
		plist_id  = H5Pcreate (H5P_DATASET_CREATE);
		
		if(data.deflate){ 
			
			if(data.debug_bit>2) mexPrintf("deflating data at level %d\n", data.deflate);
			if(data.debug_bit>2) mexPrintf("deflating data at level %d\n", 1);
			
			const hsize_t cdims[4]={1,1,data.psfs_num_cols, data.psfs_num_rows};
			
			status = H5Pset_chunk(plist_id, 4, cdims);
			if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:set_chunk", "Something wrong when setting chunk for PSFs!"); }
			
			status = H5Pset_deflate(plist_id, 1); // replaced "deflate" with 1 because having high deflate values for PSFs causes a crash...
			if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:set_deflate", "Something wrong when setting deflate for PSFs!"); }
			
		}
		
		dataset_id = H5Dcreate(file_id, "/psfs", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT); 
		if(dataset_id<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:create_dataset", "Something wrong when creating dataset for PSFs!"); }
		
		status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.psfs);
		if(data.debug_bit>2) mexPrintf("PSFS: file_id= %d | dataspace_id= %d | dataset_id= %d | status= %d\n", file_id, dataspace_id, dataset_id, status);
		if(status<0){ mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:create_dataset", "Something wrong when creating dataset for PSFs!"); }
				
		// ********* write the psf_sampling attribute *********
		hsize_t dims=1;
		att_dataspace_id=H5Screate_simple(1, &dims, NULL);

		hid_t attribute_id=H5Acreate2(dataset_id, "psf_sampling", H5T_NATIVE_DOUBLE, att_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
		status=H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &data.psf_sampling);
		if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:write_psf_sampling", "Something wrong when writing psf_sampling attribute!"); }
		
		status = H5Aclose(attribute_id);
		if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:close_att", "Something wrong when closing psf_sampling attribute!"); }				
		
		status = H5Sclose(att_dataspace_id);
		if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:close_att_dataspace", "Something wrong when closing psf_sampling attribute dataspace!"); }
		
		H5Pclose(plist_id);	
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		
    }	
	
	/////////////// write lightcurves ///////////////
	if(data.lightcurves_length){ 
			
		if(data.debug_bit>9) mexPrintf("writing lightcurves...\n");
	
		hsize_t lcs_dims_c[2] = {data.lc_num_cols, data.lc_num_rows}; // c style array with a flipped order of dims (c style!)
					
		dataspace_id = H5Screate_simple(2,lcs_dims_c, NULL);
		if(dataspace_id<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:create_dataspace", "Something wrong when creating dataspace for LCs!"); }
		
		plist_id  = H5Pcreate (H5P_DATASET_CREATE);
			
		dataset_id = H5Dcreate(file_id, "/lightcurves", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT); 
		if(dataset_id<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:create_dataset", "Something wrong when creating dataset for LCs!"); }
		
		status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.lightcurves);
		if(data.debug_bit>2) mexPrintf("LIGHTCURVES: file_id= %d | dataspace_id= %d | dataset_id= %d | status= %d\n", file_id, dataspace_id, dataset_id, status);
		if(status<0){ data.mex_flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:write_dataset", "Something wrong when writing dataset for LCs!"); }
		
		H5Pclose(plist_id);	
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
	}
	
	writeHDF5metadata(file_id, "/pars", data.parameters, data.mex_flag, data.debug_bit);
	writeHDF5metadata(file_id, "/sim_pars", data.sim_pars, data.mex_flag, data.debug_bit);
	
	H5Fclose(file_id);
	
	data.mex_flag[1]=1;
	
	
}

void writeHDF5metadata(hid_t file_id, const char *par_name, std::vector<Attribute> attributes, double *flag, int debug_bit){
	
	if(attributes.empty()) return; // short circuit in case no pars are given...
	
	if(debug_bit>6) printf("writing metadata...\n");
	
	// hsize_t size=1;
	// hid_t dataspace_id = H5Screate_simple(1, &size, NULL);
	// hid_t dataset_id = H5Dcreate(file_id, par_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); // can use H5T_NATIVE_INT / H5T_STD_U16BE
	hid_t group_id = H5Gcreate2(file_id, par_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	herr_t status;

	// start going over the parameter fields
	size_t N=attributes.size();
	hsize_t dims = 1;
	
	for(int i=0;i<N;i++){
		
		if(attributes[i].is_str){ 
			
			hsize_t num_strings=1; 
			hid_t att_dataspace_id=H5Screate_simple(1, &num_strings, NULL);
			
			if(debug_bit>9) mexPrintf("name: %s dims= %d\n", attributes[i].name.c_str(), attributes[i].str.size());
						
			hid_t atype = H5Tcopy(H5T_C_S1);
			H5Tset_size(atype, attributes[i].str.size());
			H5Tset_strpad(atype,H5T_STR_NULLTERM);
			
			hid_t attribute_id=H5Acreate2(group_id, attributes[i].name.c_str(), atype, att_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
			
			status=H5Awrite(attribute_id, atype, attributes[i].str.c_str());
			if(status<0){ flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:write_att_char", "Something wrong when writing string attribute!"); }
			
			status = H5Aclose(attribute_id);
			if(status<0){ flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:close_att", "Something wrong when closing string attribute!"); }
			
			status = H5Sclose(att_dataspace_id);
			if(status<0){ mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:close_att_dataspace", "Something wrong when closing string attribute dataspace!"); }
			
		}
		else if(attributes[i].is_scalar){ 
			
			dims=1;
			if(debug_bit>9) mexPrintf("name: %s (scalar)\n", attributes[i].name.c_str());
			hid_t att_dataspace_id=H5Screate_simple(1, &dims, NULL);
						
			hid_t attribute_id=H5Acreate2(group_id, attributes[i].name.c_str(), H5T_NATIVE_DOUBLE, att_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
			
			status=H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &attributes[i].scalar);
			if(status<0){ flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:write_att_numeric", "Something wrong when writing scalar attribute!"); }
			
			status = H5Aclose(attribute_id);
			if(status<0){flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:close_att", "Something wrong when closing attribute!"); }
			
			status = H5Sclose(att_dataspace_id);
			if(status<0){ flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:close_att_dataspace", "Something wrong when closing scalar attribute dataspace!"); }
			
		}
		else if(attributes[i].is_vec){
			
			dims=attributes[i].vec.size();	
			if(debug_bit>9) mexPrintf("name: %s dims=", attributes[i].name.c_str(), dims);
			hid_t att_dataspace_id=H5Screate_simple(1, &dims, NULL);
			
			hid_t attribute_id=H5Acreate2(group_id, attributes[i].name.c_str(), H5T_NATIVE_DOUBLE, att_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
			
			status=H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, attributes[i].vec.data());
			if(status<0){ flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:write_att_numeric", "Something wrong when writing vector attribute!"); }
			
			status = H5Aclose(attribute_id);
			if(status<0){flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:close_att", "Something wrong when closing attribute!"); }
			
			status = H5Sclose(att_dataspace_id);
			if(status<0){ flag[2]=-1; mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:writeHDF5:close_att_dataspace", "Something wrong when closing vector attribute dataspace!"); }
						
		}
		
		
	}

	// finishup

    // status = H5Sclose(dataspace_id);
	// status = H5Dclose(dataset_id);
	status = H5Gclose(group_id);

}

*/